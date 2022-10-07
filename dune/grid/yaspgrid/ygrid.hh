// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRID_YGRID_HH
#define DUNE_GRID_YASPGRID_YGRID_HH

#include <array>
#include <vector>
#include <bitset>
#include <deque>

#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/streamoperators.hh>

/** \file
    \brief This provides a YGrid, the elemental component of the yaspgrid implementation
 */

namespace Dune {

 namespace Yasp {
  /** @returns an array containing the sizes of the grids associated with vectors in given array.
   *  Needed in this form due to the need of such functionality in class initializer lists.
   *  @param v the array of vectors to examine
   */
  template<int d, typename ct>
  std::array<int,d> sizeArray(const std::array<std::vector<ct>,d>& v)
  {
    std::array<int,d> tmp;
    for (int i=0; i<d; ++i)
      tmp[i] = v[i].size() - 1;
    return tmp;
  }
 } //namespace Yasp

  /**
     The YGrid considered here describes a finite set \f$d\f$-tupels of the form
     \f[ G = \{ (k_0,\ldots,k_{d-1}) | o_i \leq k_i < o_i+s_i \}  \f]

     together with an affine mapping.

     A YGrid is characterized by the following quantities:

     - The origin \f$ o=(o_0,\ldots,o_{d-1}) \in Z^d\f$,
     - the size \f$ s=(s_0,\ldots,s_{d-1}) \in Z^d\f$,
     - The shift \f$ r=(r_0,\ldots,r_{d-1}) \in R^d\f$.
     - a coordinate container, that gives the mapping of the index to global coordinates (see coordinates.hh)

     The shift can be used to interpret the points of a grid as midpoints of cells, faces, edges, etc.

     Here is a graphical illustration of a grid:

     \image html  grid.png "A YGrid."
     \image latex grid.eps "A YGrid." width=\textwidth

     A YGrid allows to iterate over all its cells with an Iterator class.

     A YGrid is always considered as being embedded in a larger grid.
     This embedding is characterized by an offset and an enclosing grid as
     shown in the following picture:

     \image html  subgrid.png "The SubYGrid is shown in red, blue is the enclosing grid."
     \image latex subgrid.eps "The SubYGrid is shown in red, blue is the enclosing grid." width=\textwidth

     The iterator provides also a mapping to the consecutive index in the enclosing grid.

     Note: as of november 2013 there are only YGrid and YGrid::Iterator. These represent
     the functionality of former SubYGrid and SubYGrid::TransformingSubIterator. All other
     classes in the hierarchy have not been used.
   */
  template<class Coordinates>
  class YGridComponent
  {
  public:
    //extract coordinate type and dimension from the coordinate container
    typedef typename Coordinates::ctype ct;
    static const int d = Coordinates::dimension;

    typedef std::array<int, d> iTupel;
    typedef FieldVector<ct,d> fTupel;

    //! make uninitialized ygrid
    YGridComponent () : _shift(0ULL)
    {
      std::fill(_origin.begin(), _origin.end(), 0);
      std::fill(_offset.begin(), _offset.end(), 0);
      std::fill(_size.begin(), _size.end(), 0);
    }

    /** @brief make ygrid without coordinate information
     *  @param origin origin of the grid in global coordinates
     *  @param size size of the grid
     *  Such grid has no coordinate information stored but can be
     *  used to determine an intersection with a grid with coordinate
     *  information. This avoids sending coordinates in the parallel case.
     */
    YGridComponent(iTupel origin, iTupel size)
      : _origin(origin), _size(size)
    {}

    /** @brief make a subgrid by taking coordinates from a larger grid
     *  @param origin origin of the grid to be constructed
     *  @param size size of the grid to be constructed
     *  @param enclosing the grid to take coordinates and shift vector from
     */
    YGridComponent (iTupel origin, iTupel size, const YGridComponent<Coordinates>& enclosing)
      :  _origin(origin), _shift(enclosing.shift()), _coords(enclosing.getCoords()), _size(size), _supersize(enclosing.supersize())
    {
      for (int i=0; i<d; i++)
        _offset[i] = origin[i] - enclosing.origin(i) + enclosing.offset(i);

      // compute superincrements
      int inc = 1;
      for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= _supersize[i];
        }
    }

    /** @brief Make YGridComponent by giving all parameters
     *  @param origin the origin of the grid in global coordinates
     *  @param shift the shift vector
     *  @param coords the coordinate vectors to be used
     *  @param size the size vector
     *  @param offset the offset in the enclosing grid
     *  @param supersize size of the enclosing grid
     */
    YGridComponent (iTupel origin, std::bitset<d> shift, Coordinates* coords, iTupel size, iTupel offset, iTupel supersize)
      : _origin(origin), _shift(shift), _coords(coords), _size(size), _offset(offset), _supersize(supersize)
    {
      // compute superincrements
      int inc = 1;
      for (int i=0; i<d; ++i)
        {
          _superincrement[i] = inc;
          inc *= _supersize[i];
        }
    }

    //! Return origin in direction i
    int origin (int i) const
    {
      return _origin[i];
    }

    //! return reference to origin
    const iTupel& origin () const
    {
      return _origin;
    }

    //! Return shift in direction i
    bool shift (int i) const
    {
      return _shift[i];
    }

    //! Return shift tupel
    const std::bitset<d>& shift () const
    {
      return _shift;
    }

    Coordinates* getCoords() const
    {
      return _coords;
    }

    //! Return offset to origin of enclosing grid
    int offset (int i) const
    {
      return _offset[i];
    }

    //! Return offset to origin of enclosing grid
    const iTupel & offset () const
    {
      return _offset;
    }

    //! return size of enclosing grid
    int supersize (int i) const
    {
      return _supersize[i];
    }

    //! return size of enclosing grid
    const iTupel & supersize () const
    {
      return _supersize;
    }

    //! return size in direction i
    int size (int i) const
    {
      return _size[i];
    }

    //! retrun size
    iTupel size () const
    {
      return _size;
    }

    //! Return total size of index set which is the product of all size per direction.
    int totalsize () const
    {
      int s=1;
      for (int i=0; i<d; ++i)
        s *= size(i);
      return s;
    }

    //! Return minimum index in direction i
    int min (int i) const
    {
      return _origin[i];
    }

    //! Return maximum index in direction i
    int max (int i) const
    {
      return _origin[i] + size(i) - 1;
    }

    //! Return true if YGrid is empty, i.e. has size 0 in all directions.
    bool empty () const
    {
      for (int i=0; i<d; ++i)
      {
        if (size(i) == 0)
          return true;
      }
      return false;
    }

    //! given a coordinate, return true if it is in the grid
    bool inside (const iTupel& coord) const
    {
      for (int i=0; i<d; i++)
      {
        if ((coord[i]<_origin[i]) || (coord[i]>=_origin[i]+_size[i]))
          return false;
      }
      return true;
    }

    //! given a tupel compute its index in the lexicographic numbering
    int index (const iTupel& coord) const
    {
      int index = (coord[d-1]-_origin[d-1]);

      for (int i=d-2; i>=0; i--)
        index = index*_size[i] + (coord[i]-_origin[i]);

      return index;
    }

    //! return grid moved by the vector v
    YGridComponent<Coordinates> move (iTupel v) const
    {
      for (int i=0; i<d; i++)
        v[i] += _origin[i];
      return YGridComponent<Coordinates>(v,_size,*this);
    }

    //! Return YGridComponent of supergrid of self which is the intersection of self and another YGridComponent
    YGridComponent<Coordinates> intersection (const YGridComponent<Coordinates>& r) const
    {
      for (int i=0; i<d; i++)
      {
        //empty coordinate vectors result in empty intersections
        if (empty() || r.empty())
          return YGridComponent<Coordinates>();
      }

      iTupel neworigin;
      iTupel newsize;
      for (int i=0; i<d; ++i)
      {
        neworigin[i] = std::max(origin(i),r.origin(i));
        newsize[i] = std::min(max(i),r.max(i)) - neworigin[i] + 1;
      }

      return YGridComponent<Coordinates>(neworigin,newsize,*this);
    }


    /** Iterator class allows one to run over all cells of a grid.
     *  The cells of the grid to iterate over are numbered consecutively starting
     *  with zero. Via the index() method the iterator provides a mapping of the
     *  cells of the grid to a one-dimensional array. The number of entries
     *  in this array must be the size of the grid.
     */
    class Iterator {
    public:
      // default constructor
      Iterator () = default;

      //! Make iterator pointing to first cell in a grid.
      Iterator (const YGridComponent<Coordinates>& r) : _grid(&r)
      {
        iTupel coord(r.origin());
        reinit(r,coord);
      }

      //! Make iterator pointing to given cell in a grid.
      Iterator (const YGridComponent<Coordinates>& r, const iTupel& coord)
      {
        reinit(r,coord);
      }

      //! reinitialize iterator to given position
      void reinit (const YGridComponent<Coordinates>& r, const iTupel& coord)
      {
        // initialize to given position in index set
        for (int i=0; i<d; ++i)
          _coord[i] = coord[i];

        // move superindex to first cell in subgrid
        _superindex = 0;
        for (int i=0; i<d; ++i)
          _superindex += (r.offset(i)+coord[i]-r.origin(i))*r.superincrement(i);

        _grid = &r;
      }

      //! Return true when two iterators over the same grid are equal (!).
      bool operator== (const Iterator& i) const
      {
        return _superindex == i._superindex;
      }

      //! Return true when two iterators over the same grid are not equal (!).
      bool operator!= (const Iterator& i) const
      {
        return _superindex != i._superindex;
      }

      //! Return consecutive index in enclosing grid
      int superindex () const
      {
        return _superindex;
      }

      //! Return coordinate of the cell in direction i.
      int coord (int i) const
      {
        return _coord[i];
      }

      //! Return coordinate of the cell as reference (do not modify).
      const iTupel& coord () const
      {
        return _coord;
      }

      //! move this iterator dist cells in direction i
      void move (int i, int dist)
      {
        _coord[i] += dist;
        _superindex += dist*_grid->superincrement(i);
      }

      //! move this iterator dist cells in direction i
      void move (const iTupel& dist)
      {
        for (int i = 0; i < d; ++i)
          {
            _coord[i] += dist[i];
            _superindex += dist[i]*_grid->superincrement(i);
          }
      }

      //! Increment iterator to next cell with position.
      Iterator& operator++ ()
      {
        for (int i=0; i<d; i++)         // check for wrap around
        {
          _superindex += _grid->superincrement(i);   // move on cell in direction i
          if (++_coord[i] <= _grid->max(i))
            return *this;
          else
          {
            _coord[i] = _grid->origin(i);         // move back to origin in direction i
            _superindex -= _grid->size(i) * _grid->superincrement(i);
          }
        }
        // if we wrapped around, back to to begin(), we must put the iterator to end()
        if (_coord == _grid->origin())
        {
          for (int i=0; i<d; i++)
            _superindex += (_grid->size(i)-1) * _grid->superincrement(i);
          _superindex += _grid->superincrement(0);
        }
        return *this;
      }

      //! Return ith component of lower left corner of the entity associated with the current coordinates and shift.
      ct lowerleft(int i) const
      {
        return _grid->getCoords()->coordinate(i,_coord[i]);
      }

      //! Return lower left corner of the entity associated with the current coordinates and shift.
      fTupel lowerleft() const
      {
        fTupel ll;
        for (int i=0; i<d; i++)
          ll[i] = lowerleft(i);
        return ll;
      }

      //! Return ith component of upper right corder of the entity associated with the current coordinates and shift.
      ct upperright(int i) const
      {
        int coord = _coord[i];
        if (shift(i))
          coord++;
        return _grid->getCoords()->coordinate(i,coord);
      }

      //! Return upper right corder of the entity associated with the current coordinates and shift.
      fTupel upperright() const
      {
        fTupel ur;
        for (int i=0; i<d; i++)
          ur[i] = upperright(i);
        return ur;
      }

      //! Return meshsize in direction i
      ct meshsize (int i) const
      {
        return _grid->getCoords()->meshsize(i,_coord[i]);
      }

      //! Return meshsize of current cell as reference.
      fTupel meshsize () const
      {
        fTupel h;
        for (int i=0; i<d; i++)
          h[i] = meshsize(i);
        return h;
      }

      bool shift (int i) const
      {
        return _grid->shift(i);
      }

      std::bitset<d> shift() const
      {
        return _grid->shift();
      }

      Coordinates* coordCont() const
      {
        return _grid->getCoords();
      }

    protected:
      iTupel _coord;              //!< current position in index set
      int _superindex = 0;        //!< consecutive index in enclosing grid
      const YGridComponent<Coordinates>* _grid = nullptr;
    };


    int superindex(iTupel coord) const
    {
      // move superindex to first cell in subgrid
      int si = 0;
      for (int i=0; i<d; ++i)
        si += (offset(i)+coord[i]-origin(i))*_superincrement[i];
      return si;
    }

    int superincrement(int i) const
    {
      return _superincrement[i];
    }

    //! return iterator to first element of index set
    Iterator begin () const
    {
      return Iterator(*this);
    }

    //! return iterator to given element of index set
    Iterator begin (const iTupel& co) const
    {
      return Iterator(*this,co);
    }

    //! return subiterator to last element of index set
    Iterator end () const
    {
      iTupel last;
      for (int i=0; i<d; i++)
        last[i] = max(i);
      last[0] += 1;
      return Iterator(*this,last);
    }

  private:
    iTupel _origin;
    std::bitset<d> _shift;
    Coordinates* _coords;
    iTupel _size;
    iTupel _offset;    //!< offset to origin of the enclosing grid
    iTupel _supersize; //!< size of the enclosing grid
    iTupel _superincrement; //!< moves consecutive index by one in this direction in supergrid

  };


  //! Output operator for ygrids
  template <class Coordinates>
  inline std::ostream& operator<< (std::ostream& s, YGridComponent<Coordinates> e)
  {
    s << "Printing YGridComponent structure:" << std::endl;
    s << "Origin: " << e.origin() << std::endl;
    s << "Shift: " << e.shift() << std::endl;
    s << "Size: " << e.size() << std::endl;
    s << "Offset: " << e.offset() << std::endl;
    s << "Supersize: " << e.supersize() << std::endl;
    return s;
  }

  //! Output operator for ygrids
  template <class Coordinates>
  inline std::ostream& operator<< (std::ostream& s, typename YGridComponent<Coordinates>::Iterator& e)
  {
    s << "Printing YGridComponent Iterator:" << std::endl << "Iterator at " << e.coord() << " (index ";
    s << e.index() << "), position " << e.position();
    return s;
  }

  /** \brief implements a collection of YGridComponents which form a codimension
   * Entities of given codimension c need to be represented by d choose c YgridComponents.
   * All entities in one such component share the same set of spanning unit vectors.
   * A YGrid is used to iterate over the entire set of components the codimension
   * consists of. It doesn't hold any data, but instead holds an iterator range into
   * an array of components (which is owned by YGridLevel).
   */
  template<class Coordinates>
  class YGrid
  {
    public:
    static const int dim = Coordinates::dimension;

    // define data array iterator
    typedef YGridComponent<Coordinates>* DAI;

    typedef typename std::array<int, dim> iTupel;

    //! set start iterator in the data array
    void setBegin(DAI begin)
    {
      _begin = begin;
    }

    //! get which component belongs to a given shift vector
    int shiftmapping(const std::bitset<dim>& shift) const
    {
      return _shiftmapping[shift.to_ulong()];
    }

    //! get start iterator in the data array
    DAI dataBegin() const
    {
      return _begin;
    }

    //! get end iterator in the data array
    DAI dataEnd() const
    {
      return _end;
    }

    //! decide whether a coordinate is in the grid (depending on the component)
    bool inside(const iTupel& coord, const std::bitset<dim>& shift = std::bitset<dim>()) const
    {
      return (_begin+_shiftmapping[shift.to_ulong()])->inside(coord);
    }

    /** \brief Iterator over a collection o YGrids
     * A YGrid::Iterator is the heart of an entity in YaspGrid.
     */
    class Iterator
    {
      public:

      //! default constructor
      Iterator () = default;

      //! construct an iterator from coordinates and component
      Iterator (const YGrid<Coordinates>& yg, const std::array<int,dim>& coords, int which = 0)
        : _which(which), _yg(&yg)
      {
        _it = typename YGridComponent<Coordinates>::Iterator(*(_yg->dataBegin()+which),coords);
      }

      //! create an iterator to start or end of the codimension
      Iterator (const YGrid<Coordinates>& yg, bool end=false) : _yg(&yg)
      {
        if (end)
        {
          _it = _yg->_itends.back();
          _which = _yg->_itends.size() - 1;
        }
        else
        {
          _it = _yg->_itbegins[0];
          _which = 0;
        }
      }

      //! reinitializes an iterator, as if it was just constructed.
      void reinit(const YGrid<Coordinates>& yg, const std::array<int,dim>& coords, int which = 0)
      {
        _yg = &yg;
        _which = which;
        _it = typename YGridComponent<Coordinates>::Iterator(*(_yg->dataBegin()+which),coords);
      }

      //! return coordinate at the current position (direction i)
      int coord (int i) const
      {
        return _it.coord(i);
      }

      //! return coordinate array at the current position
      const std::array<int, dim>& coord () const
      {
        return _it.coord();
      }

      typename Coordinates::ctype lowerleft(int i) const
      {
        return _it.lowerleft(i);
      }

      Dune::FieldVector<typename Coordinates::ctype,dim> lowerleft() const
      {
        return _it.lowerleft();
      }

      typename Coordinates::ctype upperright(int i) const
      {
        return _it.upperright(i);
      }

      Dune::FieldVector<typename Coordinates::ctype,dim> upperright() const
      {
        return _it.upperright();
      }

      //! return the current meshsize in direction i
      typename Coordinates::ctype meshsize (int i) const
      {
        return _it.meshsize(i);
      }

      //! return the current meshsize vector
      Dune::FieldVector<typename Coordinates::ctype,dim> meshsize() const
      {
        return _it.meshsize();
      }

      //! return the shift in direction i
      bool shift (int i) const
      {
        return _it.shift(i);
      }

      //! return the shift vector
      std::bitset<dim> shift () const
      {
        return _it.shift();
      }

      //! return the superindex
      int superindex() const
      {
        // the offset of the current component has to be taken into account
          return _yg->_indexOffset[_which] + _it.superindex();
      }

      //! increment to the next entity jumping to next component if necessary
      Iterator& operator++ ()
      {
        if ((++_it == _yg->_itends[_which]) && (_which < _yg->_itends.size()-1))
          _it = _yg->_itbegins[++_which];
        return *this;
      }

      //! compare two iterators: component has to match
      bool operator==(const Iterator& i) const
      {
        if (_which != i._which)
          return false;
        return _it == i._it;
      }

      //! compare two iterators: component has to match
      bool operator!=(const Iterator& i) const
      {
        if (_it != i._it)
          return true;
        return _which != i._which;
      }

      //! return the current component number
      int which() const
      {
        return _which;
      }

      //! move the grid, this is only done and needed for codim 0
      void move(int i, int dist)
      {
        _it.move(i,dist);
      }

      void move(const iTupel& dist)
      {
        _it.move(dist);
      }

      Coordinates* coordCont() const
      {
        return _it.coordCont();
      }


      private:
      unsigned int _which = 0;
      const YGrid<Coordinates>* _yg = nullptr;
      typename YGridComponent<Coordinates>::Iterator _it;
    };

    //! return begin iterator for the codimension and partition the ygrid represents
    Iterator begin() const
    {
      return Iterator(*this);
    }

    //! return iterator pointint to a specified position
    Iterator begin(const std::array<int, dim>& coord, int which = 0) const
    {
      return Iterator(*this, coord, which);
    }

    //! return end iterator for the codimension and partition the ygrid represents
    Iterator end() const
    {
      return Iterator(*this,true);
    }

    int superindex(const iTupel& coord, int which) const
    {
      return _indexOffset[which] + (dataBegin()+which)->superindex(coord);
    }


    // finalize the ygrid construction by storing component iterators
    void finalize(const DAI& end, int artificialOffset = 0)
    {
      // set the end iterator in the ygrid component array
      _end = end;

      _indexOffset.push_back(artificialOffset);
      int k = 0;
      for (DAI i=_begin; i != _end; ++i, ++k)
      {
        //store begin and end iterators
        _itbegins.push_back(i->begin());
        _itends.push_back(i->end());

        // store index offset
        _indexOffset.push_back(_indexOffset.back() + i->totalsize());

        // store shift to component mapping
        _shiftmapping[i->shift().to_ulong()] = k;
      }
      _indexOffset.resize(_itends.size());
    }

    private:

    friend class YGrid<Coordinates>::Iterator;
    DAI _begin;
    DAI _end;
    std::array<int,Dune::power(2,dim)> _shiftmapping;
    std::vector<typename YGridComponent<Coordinates>::Iterator> _itbegins;
    std::vector<typename YGridComponent<Coordinates>::Iterator> _itends;
    std::vector<int> _indexOffset;
  };

  //! Output operator for ygrids
  template <class Coordinates>
  inline std::ostream& operator<< (std::ostream& s, const YGrid<Coordinates>& e)
  {
    s << "Printing YGrid structure:" << std::endl;
    for (auto it = e.dataBegin(); it != e.dataEnd(); ++it)
      s << *it << std::endl;
    return s;
  }

  /** \brief implements a collection of multiple std::deque<Intersection>
   * Intersections with neighboring processors are stored as std::deque<Intersection>.
   * Eachsuch intersection only holds one YGridComponent. To do all communication
   * associated with one codimension, multiple such deques have to be concatenated.
   * YGridList manges this concatenation. As for YGrids, YGridList doesn't hold any
   * data, but an iterator range into a data array owned by YGridLevel.
   */
  template<class Coordinates>
  class YGridList
  {
    public:
    static const int dim = Coordinates::dimension;

    /** \brief type describing an intersection with a neighboring processor */
    struct Intersection
    {
      /** \brief The intersection as a subgrid of the local grid */
      YGridComponent<Coordinates> grid;
      /** \brief Rank of the process where the other grid is stored */
      int rank;
      /** \brief Manhattan distance to the other grid */
      int distance;
      /** \brief a YGrid stub, that acts wraps above YGrid Component and handels the index offset */
      YGrid<Coordinates> yg;
    };

    // define data array iterator type
    typedef typename std::array<std::deque<Intersection>, Dune::power(2,dim)>::iterator DAI;

    // iterator that allows to iterate over a concatenation of deques. namely those
    // that belong to the same codimension.
    class Iterator
    {
      public:

      //! return iterator to begin and end of the container
        Iterator(const YGridList<Coordinates>& ygl, bool end=false) : _end(ygl.dataEnd()), _which(ygl.dataBegin())
      {
        _it = _which->begin();

        // advance the iterator to the first element that exists.
        // some deques might be empty and should be skipped
        while ((_which != _end) && (_it == _which->end()))
        {
          ++_which;
          if (_which != _end)
            _it = _which->begin();
        }
        // the iterator is at the end if and only if _which==_end
        if (end)
        {
          _which = _end;
        }
      }

      //! increment iterator
      Iterator& operator++ ()
      {
        ++_it;
        // advance the iterator to the next element that exists.
        // some deques might be empty and should be skipped
        while ((_which != _end) && (_it == _which->end()))
        {
          ++_which;
          if (_which != _end)
            _it = _which->begin();
        }
        return *this;
      }

      //! dereference iterator
      typename std::deque<Intersection>::iterator  operator->() const
      {
        return _it;
      }

      //! dereference iterator
      typename std::deque<Intersection>::iterator  operator*() const
      {
        return _it;
      }

      //! compare two iterators
      bool operator== (const Iterator& i) const
      {
        if (_which != i._which)
          return false;
        if (_which == _end)
          return true;
        return _it == i._it;
      }

      //! compare two iterators
      bool operator!= (const Iterator& i) const
      {
        if (_which != i._which)
          return true;
        if (_which == _end)
          return false;
        return _it != i._it;
      }

      private:
      typename std::deque<Intersection>::iterator _it;
      DAI _end;
      DAI _which;
    };

    //! return iterator pointing to the begin of the container
    Iterator begin() const
    {
      return Iterator(*this);
    }

    //! return iterator pointing to the end of the container
    Iterator end() const
    {
      return Iterator(*this,true);
    }

    //! set start iterator in the data array
    void setBegin(typename std::array<std::deque<Intersection>, Dune::power(2,dim)>::iterator begin)
    {
      _begin = begin;
    }

    //! get start iterator in the data array
    DAI dataBegin() const
    {
      return _begin;
    }

    //! get end iterator in the data array
    DAI dataEnd() const
    {
      return _end;
    }

    //! return the size of the container, this is the sum of the sizes of all deques
    int size() const
    {
      int count = 0;
      for (DAI it = _begin; it != _end; ++it)
        count += it->size();
      return count;
    }

    //! finalize the YGridLIst
    void finalize(DAI end, const YGrid<Coordinates>& ygrid)
    {
      // Instead of directly iterating over the intersection deques, this code
      // iterates over the components of an associated ygrid and works its way
      // through the list of intersection deques in parallel.
      // The reason for this convoluted iteration technique is that there are not
      // necessarily intersections for all possible shifts, but we have to make
      // sure that we stop at each shift to update the per-component index shift.
      // This is ensured by iterating over the ygrid, which is guaranteed to have
      // a component for each shift vector.

      // set end iterator in the data array
      _end = end;

      //! set offsets allow the YGridComponents in the Intersctions to behave as YGrids
      int offset = 0;

      DAI i = _begin;

      // make sure that we have a valid deque (i.e. a non-empty one)
      while (i != _end && i->begin() == i->end())
        ++i;

      for (auto yit = ygrid.dataBegin(); yit != ygrid.dataEnd(); ++yit)
      {
        if (i == _end)
          break;
        auto it = i->begin();
        if (it->grid.shift() == yit->shift())
        {
          // iterate over the intersections in the deque and set the offset
          for (; it != i->end(); ++it)
          {
            it->yg.setBegin(&(it->grid));
            it->yg.finalize(&(it->grid)+1, offset);
          }

          // advance to next non-empty deque
          ++i;
          while (i != _end && i->begin() == i->end())
            ++i;
        }

        // update the offset from the ygrid component
        int add = 1;
        for (int j=0; j<dim; j++)
          add *= yit->supersize(j);
        offset += add;
      }
      assert (i == end);
    }

    private:
    DAI _begin;
    DAI _end;
  };

} // namespace Dune

#endif
