// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDINTERSECTION_HH
#define DUNE_GRID_YASPGRIDINTERSECTION_HH

/** \file
 * \brief The YaspIntersection class
 *
   YaspIntersection provides data about intersection with
   neighboring codim 0 entities.
 */

namespace Dune {

  /** \brief  YaspIntersection provides data about intersection with
     neighboring codim 0 entities.
   */
  template<class GridImp>
  class YaspIntersection
  {
    constexpr static int dim = GridImp::dimension;
    constexpr static int dimworld = GridImp::dimensionworld;
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

    friend class YaspIntersectionIterator<GridImp>;

  public:
    // types used from grids
    typedef typename GridImp::YGridLevelIterator YGLI;
    typedef typename GridImp::YGrid::Iterator I;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    void update() {

      // vector with per-direction movements
      std::array<int,dim> dist{{0}};

      // first move: back to center
      dist[_dir] = 1 - 2*_face;

      // update face info
      _dir = _count / 2;
      _face = _count % 2;

      // second move: to new neighbor
      dist[_dir] += -1 + 2*_face;

      // move transforming iterator
      _outside.transformingsubiterator().move(dist);
    }

    /*! return true if we are on the boundary of the domain
        unless we are periodic in that direction
     */
    bool boundary () const
    {
      // Coordinate of intersection in its direction
      int coord = _inside.transformingsubiterator().coord(_dir) + _face;
      if (_inside.gridlevel()->mg->isPeriodic(_dir))
        return false;
      else
        return coord == 0
               ||
               coord == _inside.gridlevel()->mg->levelSize(_inside.gridlevel()->level(),_dir);
    }

    //! return true if neighbor across intersection exists in this processor
    bool neighbor () const
    {
      // Coordinate of intersection in its direction
      int coord = _inside.transformingsubiterator().coord(_dir) + _face;
      return coord > _inside.gridlevel()->overlap[0].dataBegin()->min(_dir)
             &&
             coord <= _inside.gridlevel()->overlap[0].dataBegin()->max(_dir);
    }

    //! Yasp is always conform
    bool conforming () const
    {
      return true;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const
    {
      return Entity(_inside);
    }

    //! return Entity on the outside of this intersection
    Entity outside() const
    {
      return Entity(_outside);
    }

    //! identifier for boundary segment from macro grid
    //! (attach your boundary condition as needed)
    int boundarySegmentIndex() const
    {
      if(! boundary())
        DUNE_THROW(GridError, "called boundarySegmentIndex while boundary() == false");
      // size of local macro grid
      const std::array<int, dim> & size = _inside.gridlevel()->mg->begin()->overlap[0].dataBegin()->size();
      const std::array<int, dim> & origin = _inside.gridlevel()->mg->begin()->overlap[0].dataBegin()->origin();
      std::array<int, dim> sides;
      {
        for (int i=0; i<dim; i++)
        {
          sides[i] =
            ((_inside.gridlevel()->mg->begin()->overlap[0].dataBegin()->origin(i)
              == 0)+
            (_inside.gridlevel()->mg->begin()->overlap[0].dataBegin()->origin(i) +
                      _inside.gridlevel()->mg->begin()->overlap[0].dataBegin()->size(i)
                      ==
                      _inside.gridlevel()->mg->levelSize(0,i)));

        }
      }
      // global position of the cell on macro grid
      std::array<int, dim> pos = _inside.transformingsubiterator().coord();
      for(int i=0; i<dim; i++)
      {
        pos[i] = pos[i] / (1<<_inside.level());
        pos[i] = pos[i] - origin[i];
      }
      // compute unit-cube-face-sizes
      std::array<int, dim> fsize;
      {
        int vol = 1;
        for (int k=0; k<dim; k++)
          vol *= size[k];
        for (int k=0; k<dim; k++)
          fsize[k] = vol/size[k];
      }
      // compute index in the unit-cube-face
      int index = 0;
      {
        int localoffset = 1;
        for (int k=dim-1; k>=0; k--)
        {
          if (k == _dir) continue;
          index += (pos[k]) * localoffset;
          localoffset *= size[k];
        }
      }
      // add unit-cube-face-offsets
      {
        for (int k=0; k<_dir; k++)
          index += sides[k] * fsize[k];
        // add fsize if we are on the right face and there is a left-face-boundary on this processor
        index += _face * (sides[_dir]>1) * fsize[_dir];
      }

      return index;
    }

    //! return unit outer normal, this should be dependent on local coordinates for higher order boundary
    FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& /* local */) const
    {
      return centerUnitOuterNormal();
    }

    //! return unit outer normal, this should be dependent on local coordinates for higher order boundary
    FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& /* local */) const
    {
      return centerUnitOuterNormal();
    }

    //! return unit outer normal at center of intersection geometry
    FieldVector<ctype, dimworld> centerUnitOuterNormal () const
    {
      FieldVector<ctype, dimworld> normal(0);
      normal[_dir] = (_face==0) ? -1.0 : 1.0;
      return normal;
    }

    //! return unit outer normal, this should be dependent on
    //! local coordinates for higher order boundary
    //! the normal is scaled with the integration element of the intersection.
    FieldVector<ctype, dimworld> integrationOuterNormal ([[maybe_unused]] const FieldVector<ctype, dim-1>& local) const
    {
      return geometry().volume() * centerUnitOuterNormal();
    }

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of the element where iteration started.
     */
    LocalGeometry geometryInInside () const
    {
      // set of dimensions that span the intersection
      std::bitset<dim> s;
      s.set();
      s[_dir] = false;

      // lower-left and upper-right corners
      Dune::FieldVector<ctype, dim> ll(0.0);
      Dune::FieldVector<ctype, dim> ur(1.0);

      ll[_dir] = ur[_dir] = (_face==0) ? 0.0 : 1.0;

      return LocalGeometry(LocalGeometryImpl(ll,ur,s));
    }

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of neighbor
     */
    LocalGeometry geometryInOutside () const
    {
      // set of dimensions that span the intersection
      std::bitset<dim> s;
      s.set();
      s[_dir] = false;

      // lower-left and upper-right corners
      Dune::FieldVector<ctype, dim> ll(0.0);
      Dune::FieldVector<ctype, dim> ur(1.0);

      ll[_dir] = ur[_dir] = (_face==1) ? 0.0 : 1.0;

      return LocalGeometry(LocalGeometryImpl(ll,ur,s));
    }

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
     */
    Geometry geometry () const
    {

      std::bitset<dim> shift;
      shift.set();
      shift[_dir] = false;

      Dune::FieldVector<ctype,dimworld> ll, ur;
      for (int i=0; i<dimworld; i++)
      {
        int coord = _inside.transformingsubiterator().coord(i);

        if ((i == _dir) and (_face))
          coord++;

        ll[i] = _inside.transformingsubiterator().coordCont()->coordinate(i,coord);
        if (i != _dir)
          coord++;
        ur[i] = _inside.transformingsubiterator().coordCont()->coordinate(i,coord);

        // If on periodic overlap, transform coordinates by domain size
        if (_inside.gridlevel()->mg->isPeriodic(i)) {
          int coordPeriodic = _inside.transformingsubiterator().coord(i);
          if (coordPeriodic < 0) {
            auto size = _inside.gridlevel()->mg->domainSize()[i];
            ll[i] += size;
            ur[i] += size;
          } else if (coordPeriodic + 1 > _inside.gridlevel()->mg->levelSize(_inside.gridlevel()->level(),i)) {
            auto size = _inside.gridlevel()->mg->domainSize()[i];
            ll[i] -= size;
            ur[i] -= size;
          }
        }
      }

      GeometryImpl _is_global(ll,ur,shift);
      return Geometry( _is_global );
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return GeometryTypes::cube(dim-1);
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return _count;
    }

    //! local index of codim 1 entity in neighbor where intersection is contained in
    int indexInOutside () const
    {
      // flip the last bit
      return _count^1;
    }

    YaspIntersection()
      : _count(~std::uint8_t(0)) // Use as marker for invalid intersection
      , _dir(0)
      , _face(0)
    {}

    //! make intersection iterator from entity, initialize to first neighbor
    YaspIntersection (const YaspEntity<0,dim,GridImp>& myself, bool toend) :
      _inside(myself.gridlevel(),
              myself.transformingsubiterator()),
      _outside(myself.gridlevel(),
               myself.transformingsubiterator()),
      // initialize to first neighbor
      _count(0),
      _dir(0),
      _face(0)
    {
      if (toend)
      {
        // initialize end iterator
        _count = 2*dim;
        return;
      }
      _count = 0;

      // move transforming iterator
      _outside.transformingsubiterator().move(_dir,-1);
    }

    //! copy constructor -- use default

    //! copy operator - use default
    void assign (const YaspIntersection& it)
    {
      *this = it;
    }

    bool equals(const YaspIntersection& other) const
    {
      // compare counts first -- that's cheaper if the test fails
      return _count == other._count && _inside.equals(other._inside);
    }

  private:
    /* The two entities that make up the intersection */
    YaspEntity<0,GridImp::dimension,GridImp> _inside;  //!< entity from which the intersection was created from
    YaspEntity<0,GridImp::dimension,GridImp> _outside; //!< outside entity
    /* current position */
    std::uint8_t _count; //!< valid neighbor count in 0 .. 2*dim-1
    std::uint8_t _dir; //!< count/2
    std::uint8_t _face; //!< count%2
  };
}   // namespace Dune

#endif   // DUNE_GRID_YASPGRIDINTERSECTION_HH
