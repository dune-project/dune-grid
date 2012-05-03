// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_YASPGRID_HH
#define DUNE_YASPGRID_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>

// either include stdint.h or provide fallback for uint8_t
#if HAVE_STDINT_H
#include <stdint.h>
#else
typedef unsigned char uint8_t;
#endif

#include <dune/grid/common/grid.hh>     // the grid base classes
#include <dune/grid/yaspgrid/grids.hh>  // the yaspgrid base classes
#include <dune/grid/common/capabilities.hh> // the capabilities
#include <dune/common/misc.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpihelper.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/datahandleif.hh>


#if HAVE_MPI
#include <dune/common/mpicollectivecommunication.hh>
#endif

/*! \file yaspgrid.hh
   YaspGrid stands for yet another structured parallel grid.
   It will implement the dune grid interface for structured grids with codim 0
   and dim, with arbitrary overlap, parallel features with two overlap
   models, periodic boundaries and fast a implementation allowing on-the-fly computations.
 */

namespace Dune {

  //************************************************************************
  /*! define name for floating point type used for coordinates in yaspgrid.
     You can change the type for coordinates by changing this single typedef.
   */
  typedef double yaspgrid_ctype;
  static const yaspgrid_ctype yasptolerance=1E-13; // tolerance in coordinate computations

  /* some sizes for building global ids
   */
  const int yaspgrid_dim_bits = 24; // bits for encoding each dimension
  const int yaspgrid_level_bits = 6; // bits for encoding level number
  const int yaspgrid_codim_bits = 4; // bits for encoding codimension


  //************************************************************************
  // forward declaration of templates

  template<int dim>                             class YaspGrid;
  template<int mydim, int cdim, class GridImp>  class YaspGeometry;
  template<int codim, int dim, class GridImp>   class YaspEntity;
  template<int codim, class GridImp>            class YaspEntityPointer;
  template<int codim, class GridImp>            class YaspEntitySeed;
  template<int codim, PartitionIteratorType pitype, class GridImp> class YaspLevelIterator;
  template<class GridImp>            class YaspIntersectionIterator;
  template<class GridImp>            class YaspIntersection;
  template<class GridImp>            class YaspHierarchicIterator;
  template<class GridImp>            class YaspLevelIndexSet;
  template<class GridImp>            class YaspLeafIndexSet;
  template<class GridImp>            class YaspGlobalIdSet;

  namespace FacadeOptions
  {

    template<int dim, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim, YaspGrid<dim>, YaspGeometry>
    {
      static const bool v = false;
    };

    template<int dim, int mydim, int cdim>
    struct StoreGeometryReference<mydim, cdim, const YaspGrid<dim>, YaspGeometry>
    {
      static const bool v = false;
    };

  }

  //========================================================================
  // The transformation describing the refinement rule

  template<int dim, class GridImp>
  struct YaspFatherRelativeLocalElement {
    static const array<YaspGeometry<dim,dim,GridImp>, (1<<dim) > _geo;
    static array<YaspGeometry<dim,dim,GridImp>, (1<<dim) > initSons()
    {
      array<YaspGeometry<dim,dim,GridImp>, (1<<dim) > geo;
      FieldVector<yaspgrid_ctype,dim> midpoint(0.25);
      FieldVector<yaspgrid_ctype,dim> extension(0.5);
      for (int i=0; i<(1<<dim); i++)
      {
        midpoint = 0.25;
        for (int k=0; k<dim; k++)
        {
          if (i&(1<<k))
            midpoint[k] = 0.75;
        }
        geo[i] = YaspGeometry<dim,dim,GridImp>(midpoint, extension);
      }
      return geo;
    }
  };

  template<int dim, class GridImp>
  const array<YaspGeometry<dim,dim,GridImp>, (1<<dim)>
  YaspFatherRelativeLocalElement<dim, GridImp>::_geo =
    YaspFatherRelativeLocalElement<dim, GridImp>::initSons();

  //========================================================================
  /*!
     YaspGeometry realizes the concept of the geometric part of a mesh entity.

     We have specializations for dim == dimworld (elements) and dim == 0
     (vertices).  The general version implements dim == dimworld-1 (faces)
     and otherwise throws a GridError.
   */
  //========================================================================

  //! The general version implements dim==dimworld-1. If this is not the case an error is thrown
  template<int mydim,int cdim, class GridImp>
  class YaspGeometry : public GeometryDefaultImplementation<mydim,cdim,GridImp,YaspGeometry>
  {
  public:
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! return the element type identifier
    GeometryType type () const
    {
      return GeometryType(GeometryType::cube,mydim);
    }

    //! here we have always an affine geometry
    bool affine() const { return true; }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const
    {
      return 1<<mydim;
    }

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector< ctype, cdim > corner ( const int i ) const
    {
      assert( i >= 0 && i < (int) coord_.N() );
      FieldVector<ctype, cdim>& c = coord_[i];
      int bit=0;
      for (int k=0; k<cdim; k++) // run over all directions in world
      {
        if (k==missing)
        {
          c[k] = midpoint[k];
          continue;
        }
        //k is not the missing direction
        if (i&(1<<bit)) // check whether bit is set or not
          c[k] = midpoint[k]+0.5*extension[k]; // bit is 1 in i
        else
          c[k] = midpoint[k]-0.5*extension[k]; // bit is 0 in i
        bit++; // we have processed a direction
      }

      return c;
    }

    //! access to the center/centroid
    FieldVector< ctype, cdim > center ( ) const
    {
      return midpoint;
    }

    //! maps a local coordinate within reference element to global coordinate in element
    FieldVector<ctype, cdim> global (const FieldVector<ctype, mydim>& local) const
    {
      FieldVector<ctype, cdim> g;
      int bit=0;
      for (int k=0; k<cdim; k++)
        if (k==missing)
          g[k] = midpoint[k];
        else
        {
          g[k] = midpoint[k] + (local[bit]-0.5)*extension[k];
          bit++;
        }
      return g;
    }

    //! maps a global coordinate within the element to a local coordinate in its reference element
    FieldVector<ctype, mydim> local (const FieldVector<ctype, cdim>& global) const
    {
      FieldVector<ctype, mydim> l; // result
      int bit=0;
      for (int k=0; k<cdim; k++)
        if (k!=missing)
        {
          l[bit] = (global[k]-midpoint[k])/extension[k] + 0.5;
          bit++;
        }
      return l;
    }

    //! return volume of geometry
    ctype volume () const
    {
      ctype volume=1.0;
      for (int k=0; k<cdim; k++)
        if (k!=missing) volume *= extension[k];
      return volume;
    }

    /*! determinant of the jacobian of the mapping
     */
    ctype integrationElement (const FieldVector<ctype, mydim>& local) const
    {
      return volume();
    }

    //! Compute the transposed of the jacobi matrix
    FieldMatrix<ctype,mydim,cdim>& jacobianTransposed (const FieldVector<ctype, mydim>& local) const
    {
      JT = 0.0;
      int k=0;
      for (int i=0; i<cdim; ++i)
      {
        if (i != missing)
        {
          JT[k][i] = extension[i]; // set diagonal element
          k++;
        }
      }
      return JT;
    }
    //! Compute the transposed of the inverse jacobi matrix
    FieldMatrix<ctype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const
    {
      Jinv = 0.0;
      int k=0;
      for (int i=0; i<cdim; ++i)
      {
        if (i != missing)
        {
          Jinv[i][k] = 1.0/extension[i]; // set diagonal element
          k++;
        }
      }
      return Jinv;
    }

    //! default constructor
    YaspGeometry () {}

    //! constructor from midpoint and extension and missing direction number
    YaspGeometry (const FieldVector<ctype, cdim>& p, const FieldVector<ctype, cdim>& h, uint8_t& m)
      : midpoint(p), extension(h), missing(m)
    {
      if (cdim!=mydim+1)
        DUNE_THROW(GridError, "general YaspGeometry assumes cdim=mydim+1");
    }

    //! copy constructor (skipping temporary variables)
    YaspGeometry (const YaspGeometry& other)
      : midpoint(other.midpoint),
        extension(other.extension),
        missing(other.missing)
    {}

    //! print function
    void print (std::ostream& s) const
    {
      s << "YaspGeometry<"<<mydim<<","<<cdim<< "> ";
      s << "midpoint";
      for (int i=0; i<cdim; i++)
        s << " " << midpoint[i];
      s << " extension";
      for (int i=0; i<cdim; i++)
        s << " " << extension[i];
      s << " missing is " << missing;
    }

    // const YaspGeometry<mydim,cdim,GridImp>&
    // operator = (const YaspGeometry<mydim,cdim,GridImp>& g);

  private:
    // the element is fully defined by its midpoint the extension
    // in each direction and the missing direction.
    // Note cdim == mydim+1

    FieldVector<ctype, cdim> midpoint;  // the midpoint
    FieldVector<ctype, cdim> extension; // the extension
    uint8_t missing;                    // the missing, i.e. constant direction

    // In addition we need memory in order to return references.
    // Possibly we should change this in the interface ...
    mutable FieldMatrix<ctype, mydim, cdim> JT;   // the transposed of the jacobian
    mutable FieldMatrix<ctype, cdim, mydim> Jinv; // the transposed of the jacobian inverse
    mutable FieldMatrix<ctype, Power_m_p<2,mydim>::power, cdim> coord_; // the coordinates

  };



  //! specialize for dim=dimworld, i.e. a volume element
  template<int mydim, class GridImp>
  class YaspGeometry<mydim,mydim,GridImp> : public GeometryDefaultImplementation<mydim,mydim,GridImp,YaspGeometry>
  {
  public:
    typedef typename GridImp::ctype ctype;

    //! return the element type identifier
    GeometryType type () const
    {
      return GeometryType(GeometryType::cube,mydim);
    }

    //! here we have always an affine geometry
    bool affine() const { return true; }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const
    {
      return 1<<mydim;
    }

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<ctype, mydim>& operator[] (int i) const
    {
      return corner(i);
    }

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector< ctype, mydim > corner ( const int i ) const
    {
      assert( i >= 0 && i < (int) coord_.N() );
      FieldVector<ctype, mydim>& c = coord_[i];
      for (int k=0; k<mydim; k++)
        if (i&(1<<k))
          c[k] = midpoint[k]+0.5*extension[k]; // kth bit is 1 in i
        else
          c[k] = midpoint[k]-0.5*extension[k]; // kth bit is 0 in i
      return c;
    }

    //! access to the center/centroid
    FieldVector< ctype, mydim > center ( ) const
    {
      return midpoint;
    }

    //! maps a local coordinate within reference element to global coordinate in element
    FieldVector<ctype, mydim> global (const FieldVector<ctype, mydim>& local) const
    {
      FieldVector<ctype,mydim> g;
      for (int k=0; k<mydim; k++)
        g[k] = midpoint[k] + (local[k]-0.5)*extension[k];
      return g;
    }

    //! maps a global coordinate within the element to a local coordinate in its reference element
    FieldVector<ctype, mydim> local (const FieldVector<ctype,mydim>& global) const
    {
      FieldVector<ctype, mydim> l; // result
      for (int k=0; k<mydim; k++)
        l[k] = (global[k]-midpoint[k])/extension[k] + 0.5;
      return l;
    }

    /*! determinant of the jacobian of the mapping
     */
    ctype integrationElement (const FieldVector<ctype, mydim>& local) const
    {
      return volume();
    }

    //! return volume of geometry
    ctype volume () const
    {
      ctype vol=1.0;
      for (int k=0; k<mydim; k++) vol *= extension[k];
      return vol;
    }

    //! Compute the transposed of the jacobi matrix
    FieldMatrix<ctype,mydim,mydim>& jacobianTransposed (const FieldVector<ctype, mydim>& local) const
    {
      for (int i=0; i<mydim; ++i)
      {
        JT[i] = 0.0;              // set column to zero
        JT[i][i] = extension[i]; // set diagonal element
      }
      return JT;
    }
    //! Compute the transposed of the inverse jacobi matrix
    FieldMatrix<ctype,mydim,mydim>& jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const
    {
      for (int i=0; i<mydim; ++i)
      {
        Jinv[i] = 0.0;              // set column to zero
        Jinv[i][i] = 1.0/extension[i]; // set diagonal element
      }
      return Jinv;
    }

    //! default constructor
    YaspGeometry () {}

    //! constructor from midpoint and extension
    YaspGeometry (const FieldVector<ctype, mydim>& p, const FieldVector<ctype, mydim>& h)
      : midpoint(p), extension(h)
    {}

    //! copy constructor (skipping temporary variables)
    YaspGeometry (const YaspGeometry& other)
      : midpoint(other.midpoint),
        extension(other.extension)
    {}

    //! print function
    void print (std::ostream& s) const
    {
      s << "YaspGeometry<"<<mydim<<","<<mydim<< "> ";
      s << "midpoint";
      for (int i=0; i<mydim; i++)
        s << " " << midpoint[i];
      s << " extension";
      for (int i=0; i<mydim; i++)
        s << " " << extension[i];
    }

    // const YaspGeometry<mydim,mydim,GridImp>&
    // operator = (const YaspGeometry<mydim,mydim,GridImp>& g);

  private:
    // the element is fully defined by midpoint and the extension
    // in each direction. References are used because this information
    // is known outside the element in many cases.
    // Note mydim==cdim

    FieldVector<ctype, mydim> midpoint; // the midpoint
    FieldVector<ctype, mydim> extension; // the extension

    // In addition we need memory in order to return references.
    // Possibly we should change this in the interface ...
    mutable FieldMatrix<ctype, mydim, mydim> Jinv,JT; // the transpose of the jacobian and its inverse inverse
    mutable FieldMatrix<ctype, Power_m_p<2,mydim>::power, mydim> coord_; // the coordinates
  };

  //! specialization for dim=0, this is a vertex
  template<int cdim, class GridImp>
  class YaspGeometry<0,cdim,GridImp> : public GeometryDefaultImplementation<0,cdim,GridImp,YaspGeometry>
  {
  public:
    typedef typename GridImp::ctype ctype;

    //! return the element type identifier
    GeometryType type () const
    {
      return GeometryType(GeometryType::cube,0);
    }

    //! here we have always an affine geometry
    bool affine() const { return true; }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const
    {
      return 1;
    }

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<ctype, cdim>& operator[] (int i) const
    {
      return position;
    }

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector< ctype, cdim > corner ( const int i ) const
    {
      return position;
    }

    //! access to the center/centroid
    FieldVector< ctype, cdim > center ( ) const
    {
      return position;
    }

    /*! determinant of the jacobian of the mapping
     */
    ctype integrationElement (const FieldVector<ctype, 0>& local) const
    {
      return 1.0;
    }

    //! Compute the transposed of the jacobi matrix
    FieldMatrix<ctype,0,cdim>& jacobianTransposed (const FieldVector<ctype, 0>& local) const
    {
      static FieldMatrix<ctype,0,cdim> JT(0.0);
      return JT;
    }
    //! Compute the transposed of the inverse jacobi matrix
    FieldMatrix<ctype,cdim,0>& jacobianInverseTransposed (const FieldVector<ctype, 0>& local) const
    {
      static FieldMatrix<ctype,cdim,0> Jinv(0.0);
      return Jinv;
    }

    //! default constructor
    YaspGeometry ()
    {}

    //! constructor
    explicit YaspGeometry ( const FieldVector< ctype, cdim > &p )
      : position( p )
    {}

    YaspGeometry ( const FieldVector< ctype, cdim > &p, const FieldVector< ctype, cdim > &, uint8_t &)
      : position( p )
    {}

    //! print function
    void print (std::ostream& s) const
    {
      s << "YaspGeometry<"<<0<<","<<cdim<< "> ";
      s << "position " << position;
    }

    // const YaspGeometry<0,cdim,GridImp>&
    // operator = (const YaspGeometry<0,cdim,GridImp>& g);

  private:
    FieldVector<ctype, cdim> position; //!< position of the vertex
  };

  // operator<< for all YaspGeometrys
  template <int mydim, int cdim, class GridImp>
  inline
  std::ostream& operator<< (std::ostream& s, YaspGeometry<mydim,cdim,GridImp>& e)
  {
    e.print(s);
    return s;
  }

  //========================================================================
  /*!
     YaspEntity realizes the concept a mesh entity.

     We have specializations for codim==0 (elements) and
     codim=dim (vertices).
     The general version throws a GridError.
   */
  //========================================================================

  template<int codim, int dim, class GridImp>
  class YaspSpecialEntity :
    public GridImp::template Codim<codim>::Entity
  {
  public:
    typedef typename GridImp::ctype ctype;

    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;

    YaspSpecialEntity(const GridImp* yg, const YGLI& g, const TSI& it) :
      GridImp::template Codim<codim>::Entity (YaspEntity<codim, dim, GridImp>(yg,g,it))
    {};
    YaspSpecialEntity(const YaspEntity<codim, dim, GridImp>& e) :
      GridImp::template Codim<codim>::Entity (e)
    {};
    const TSI& transformingsubiterator () const
    {
      return this->realEntity.transformingsubiterator();
    }
    const YGLI& gridlevel () const
    {
      return this->realEntity.gridlevel();
    }
    const GridImp * yaspgrid() const
    {
      return this->realEntity.yaspgrid();
    }
  };

  template<int codim, int dim, class GridImp>
  class YaspEntity
    :  public EntityDefaultImplementation <codim,dim,GridImp,YaspEntity>
  {
  public:
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;
    //! level of this element
    int level () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! index is unique and consecutive per level and codim used for access to degrees of freedom
    int index () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! return partition type attribute
    PartitionType partitionType () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    const GridImp * yaspgrid() const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    YaspEntity (const GridImp* yg, const YGLI& g, const TSI& it)
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    // IndexSets needs access to the private index methods
    friend class Dune::YaspLevelIndexSet<GridImp>;
    friend class Dune::YaspLeafIndexSet<GridImp>;
    friend class Dune::YaspGlobalIdSet<GridImp>;
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      DUNE_THROW(GridError, "YaspEntity not implemented");
    }
  };


  // specialization for codim=0
  template<int dim, class GridImp>
  class YaspEntity<0,dim,GridImp>
    : public EntityDefaultImplementation <0,dim,GridImp,YaspEntity>
  {
    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl GeometryImpl;

  public:
    typedef typename GridImp::ctype ctype;

    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;

    typedef typename GridImp::template Codim< 0 >::Geometry Geometry;
    typedef typename GridImp::template Codim< 0 >::LocalGeometry LocalGeometry;

    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;
    typedef typename GridImp::LevelIntersectionIterator IntersectionIterator;
    typedef typename GridImp::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename GridImp::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename GridImp::HierarchicIterator HierarchicIterator;

    //! define the type used for persisitent indices
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! define type used for coordinates in grid module
    typedef typename YGrid<dim,ctype>::iTupel iTupel;

    // constructor
    YaspEntity (const GridImp * yg, const YGLI& g, const TSI& it)
      : _yg(yg), _it(it), _g(g)
    {}

    //! level of this element
    int level () const { return _g.level(); }

    //! index is unique and consecutive per level
    int index () const { return _it.superindex(); } // superindex works also for iteration over subgrids

    //! globalIndex is unique and consecutive per global level
    int globalIndex () const {
      return _g.cell_global().index(_it.coord());
    }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as less memory as possible
     */
    EntitySeed seed () const {
      return EntitySeed(_g.level(), _it.coord());
    }

    //! return partition type attribute
    PartitionType partitionType () const
    {
      if (_g.cell_interior().inside(_it.coord())) return InteriorEntity;
      if (_g.cell_overlap().inside(_it.coord())) return OverlapEntity;
      DUNE_THROW(GridError, "Impossible GhostEntity " << _it.coord() << "\t"
                                                      << _g.cell_interior().origin() << "/" << _g.cell_interior().size());
      return GhostEntity;
    }

    //! geometry of this entity
    Geometry geometry () const {
      // the element geometry
      GeometryImpl _geometry(_it.position(),_it.meshsize());
      return Geometry( _geometry );
    }

    /*! Return number of subentities with codimension cc.
     */
    template<int cc> int count () const
    {
      if (cc==dim) return 1<<dim;
      if (cc==1) return 2*dim;
      if (cc==dim-1) return dim*(1<<(dim-1));
      if (cc==0) return 1;
      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    /*! Intra-element access to subentities of codimension cc > codim.
     */
    template<int cc>
    typename Codim<cc>::EntityPointer subEntity (int i) const
    {
      dune_static_assert( cc == dim || cc == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      // coordinates of the cell == coordinates of lower left corner
      if (cc==dim)
      {
        iTupel coord = _it.coord();

        // get corner from there
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        return YaspEntityPointer<cc,GridImp>(_yg,_g,_g.vertex_overlapfront().tsubbegin(coord));
      }
      if (cc==0)
      {
        return YaspEntityPointer<cc,GridImp>(_yg,_g,_it);
      }
      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    //! Inter-level access to father element on coarser grid. Assumes that meshes are nested.
    EntityPointer father () const
    {
      // check if coarse level exists
      if (_g.level()<=0)
        DUNE_THROW(GridError, "tried to call father on level 0");

      // yes, get iterator to it
      YGLI cg = _g.coarser();

      // coordinates of the cell
      iTupel coord = _it.coord();

      // get coordinates on next coarser level
      for (int k=0; k<dim; k++) coord[k] = coord[k]/2;

      return YaspEntityPointer<0,GridImp>(_yg,cg,cg.cell_overlap().tsubbegin(coord));
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return (_g.level()>0);
    }

    /*! Location of this element relative to the reference element element of the father.
          This is sufficient to interpolate all dofs in conforming case.
          Nonconforming case may require access to neighbors of father and
          computations with local coordinates.
          On the fly case is somewhat inefficient since dofs  are visited several times.
          If we store interpolation matrices, this is tolerable. We assume that on-the-fly
          implementation of numerical algorithms is only done for simple discretizations.
          Assumes that meshes are nested.
     */
    LocalGeometry geometryInFather () const
    {
      // determine which son we are
      int son = 0;
      for (int k=0; k<dim; k++)
        if (_it.coord(k)%2)
          son += (1<<k);

      // configure one of the 2^dim transformations
      return LocalGeometry( YaspFatherRelativeLocalElement<dim,GridImp>::_geo[son] );
    }

    const TSI& transformingsubiterator () const
    {
      return _it;
    }

    const YGLI& gridlevel () const
    {
      return _g;
    }

    const GridImp* yaspgrid () const
    {
      return _yg;
    }

    bool isLeaf() const
    {
      return (_g.level() == _g.mg()->maxlevel());
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return _yg->adaptRefCount > 0 && _g.mg()->maxlevel() < _g.level() + _yg->adaptRefCount; }

    /**\brief Returns true, if entity might disappear during the next call to adapt()
     */
    bool mightVanish () const { return false; }
    // { return _yg->adaptRefCount < 0 && _g.mg()->maxlevel() < _g.level() - _yg->adaptRefCount; }

    //! returns intersection iterator for first intersection
    IntersectionIterator ibegin () const
    {
      return YaspIntersectionIterator<GridImp>(*this,false);
    }

    //! returns intersection iterator for first intersection
    LeafIntersectionIterator ileafbegin () const
    {
      // only if entity is leaf this iterator delivers intersections
      return YaspIntersectionIterator<GridImp>(*this, ! isLeaf() );
    }

    //! returns intersection iterator for first intersection
    LevelIntersectionIterator ilevelbegin () const
    {
      return ibegin();
    }

    //! Reference to one past the last neighbor
    IntersectionIterator iend () const
    {
      return YaspIntersectionIterator<GridImp>(*this,true);
    }

    //! Reference to one past the last neighbor
    LeafIntersectionIterator ileafend () const
    {
      return iend();
    }

    //! Reference to one past the last neighbor
    LevelIntersectionIterator ilevelend () const
    {
      return iend();
    }

    /*! Inter-level access to son elements on higher levels<=maxlevel.
          This is provided for sparsely stored nested unstructured meshes.
          Returns iterator to first son.
     */
    HierarchicIterator hbegin (int maxlevel) const
    {
      return YaspHierarchicIterator<GridImp>(_yg,_g,_it,maxlevel);
    }

    //! Returns iterator to one past the last son
    HierarchicIterator hend (int maxlevel) const
    {
      return YaspHierarchicIterator<GridImp>(_yg,_g,_it,_g.level());
    }

  private:
    // IndexSets needs access to the private index methods
    friend class Dune::YaspLevelIndexSet<GridImp>;
    friend class Dune::YaspLeafIndexSet<GridImp>;
    friend class Dune::YaspGlobalIdSet<GridImp>;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      // get size of global grid
      const iTupel& size =  _g.cell_global().size();

      // get coordinate correction for periodic boundaries
      int coord[dim];
      for (int i=0; i<dim; i++)
      {
        coord[i] = _it.coord(i);
        if (coord[i]<0) coord[i] += size[i];
        if (coord[i]>=size[i]) coord[i] -= size[i];
      }

      // encode codim
      PersistentIndexType id(0);

      // encode level
      id = id << yaspgrid_level_bits;
      id = id+PersistentIndexType(_g.level());


      // encode coordinates
      for (int i=dim-1; i>=0; i--)
      {
        id = id << yaspgrid_dim_bits;
        id = id+PersistentIndexType(coord[i]);
      }

      return id;
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const
    {
      return _it.superindex();
    }

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      return _it.superindex();
    }

    //! subentity persistent index
    PersistentIndexType subPersistentIndex (int i, int cc) const
    {
      if (cc==0)
        return persistentIndex();

      // get position of cell, note that global origin is zero
      // adjust for periodic boundaries
      int coord[dim];
      for (int k=0; k<dim; k++)
      {
        coord[k] = _it.coord(k);
        if (coord[k]<0) coord[k] += _g.cell_global().size(k);
        if (coord[k]>=_g.cell_global().size(k)) coord[k] -= _g.cell_global().size(k);
      }

      if (cc==dim)
      {
        // transform to vertex coordinates
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        // determine min number of trailing zeroes
        int trailing = 1000;
        for (int i=0; i<dim; i++)
        {
          // count trailing zeros
          int zeros = 0;
          for (int j=0; j<_g.level(); j++)
            if (coord[i]&(1<<j))
              break;
            else
              zeros++;
          trailing = std::min(trailing,zeros);
        }

        // determine the level of this vertex
        int level = _g.level()-trailing;

        // encode codim
        PersistentIndexType id(dim);

        // encode level
        id = id << yaspgrid_level_bits;
        id = id+PersistentIndexType(level);

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << yaspgrid_dim_bits;
          id = id+PersistentIndexType(coord[i]>>trailing);
        }

        return id;
      }

      if (cc==1) // faces, i.e. for dim=2 codim=1 is treated as a face
      {
        // Idea: Use the doubled grid to assign coordinates to faces

        // ivar is the direction that varies
        int ivar=i/2;

        // compute position from cell position
        for (int k=0; k<dim; k++)
          coord[k] = coord[k]*2 + 1; // the doubled grid
        if (i%2)
          coord[ivar] += 1;
        else
          coord[ivar] -= 1;

        // encode codim
        PersistentIndexType id(1);

        // encode level
        id = id << yaspgrid_level_bits;
        id = id+PersistentIndexType(_g.level());

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << yaspgrid_dim_bits;
          id = id+PersistentIndexType(coord[i]);
        }

        return id;
      }

      // map to old numbering
      static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
      i = edge[i];

      if (cc==dim-1) // edges, exist only for dim>2
      {
        // Idea: direction i is fixed, all others are vary, i.e. 2^(dim-1) possibilities per direction

        // number of entities per direction
        int m=1<<(dim-1);

        // ifix is the direction that is fixed
        int ifix=(dim-1)-(i/m);

        // compute position from cell position
        int bit=1;
        for (int k=0; k<dim; k++)
        {
          coord[k] = coord[k]*2+1;   // cell position in doubled grid
          if (k==ifix) continue;
          if ((i%m)&bit) coord[k] += 1;else coord[k] -= 1;
          bit *= 2;
        }

        // encode codim
        PersistentIndexType id(dim-1);

        // encode level
        id = id << yaspgrid_level_bits;
        id = id+PersistentIndexType(_g.level());

        // encode coordinates
        for (int i=dim-1; i>=0; i--)
        {
          id = id << yaspgrid_dim_bits;
          id = id+PersistentIndexType(coord[i]);
        }

        return id;
      }

      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    //! subentity compressed index
    int subCompressedIndex (int i, int cc) const
    {
      if (cc==0)
        return compressedIndex();

      // get cell position relative to origin of local cell grid
      iTupel coord;
      for (int k=0; k<dim; ++k)
        coord[k] = _it.coord(k)-_g.cell_overlap().origin(k);

      if (cc==dim) // vertices
      {
        // transform cell coordinate to corner coordinate
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          index = (index*(_g.cell_overlap().size(k)+1))+coord[k];
        return index;
      }

      if (cc==1) // faces, i.e. for dim=2 codim=1 is treated as a face
      {
        // Idea: direction ivar varies, all others are fixed, i.e. 2 possibilities per direction

        // ivar is the direction that varies
        int ivar=i/2;

        // compute position from cell position
        if (i%2) coord[ivar] += 1;

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k==ivar)
            index = (index*(_g.cell_overlap().size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g.cell_overlap().size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=0; j<ivar; j++)
        {
          int n=_g.cell_overlap().size(j)+1;
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g.cell_overlap().size(l);
          index += n;
        }

        return index;
      }

      // map to old numbering
      static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
      i = edge[i];

      if (cc==dim-1) // edges, exist only for dim>2
      {
        // Idea: direction i is fixed, all others are vary, i.e. 2^(dim-1) possibilities per direction

        // number of entities per direction
        int m=1<<(dim-1);

        // ifix is the direction that is fixed
        int ifix=(dim-1)-(i/m);

        // compute position from cell position
        int bit=1;
        for (int k=0; k<dim; k++)
        {
          if (k==ifix) continue;
          if ((i%m)&bit) coord[k] += 1;
          bit *= 2;
        }

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k!=ifix)
            index = (index*(_g.cell_overlap().size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g.cell_overlap().size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=dim-1; j>ifix; j--)
        {
          int n=_g.cell_overlap().size(j);
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g.cell_overlap().size(l)+1;
          index += n;
        }

        return index;
      }

      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    //! subentity compressed index
    int subCompressedLeafIndex (int i, int cc) const
    {
      if (cc==0)
        return compressedIndex();

      // get cell position relative to origin of local cell grid
      iTupel coord;
      for (int k=0; k<dim; ++k)
        coord[k] = _it.coord(k)-_g.cell_overlap().origin(k);

      if (cc==dim) // vertices
      {
        // transform cell coordinate to corner coordinate
        for (int k=0; k<dim; k++)
          if (i&(1<<k)) (coord[k])++;

        // move coordinates up to maxlevel
        for (int k=0; k<dim; k++)
          coord[k] = coord[k]<<(_g.mg()->maxlevel()-_g.level());

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          index = (index*(_g.mg()->rbegin().cell_overlap().size(k)+1))+coord[k];
        return index;
      }

      if (cc==1) // faces, i.e. for dim=2 codim=1 is treated as a face
      {
        // Idea: direction ivar varies, all others are fixed, i.e. 2 possibilities per direction

        // ivar is the direction that varies
        int ivar=i/2;

        // compute position from cell position
        if (i%2) coord[ivar] += 1;

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k==ivar)
            index = (index*(_g.cell_overlap().size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g.cell_overlap().size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=0; j<ivar; j++)
        {
          int n=_g.cell_overlap().size(j)+1;
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g.cell_overlap().size(l);
          index += n;
        }

        return index;
      }

      // map to old numbering
      static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
      i = edge[i];

      if (cc==dim-1) // edges, exist only for dim>2
      {
        // Idea: direction i is fixed, all others are vary, i.e. 2^(dim-1) possibilities per direction

        // number of entities per direction
        int m=1<<(dim-1);

        // ifix is the direction that is fixed
        int ifix=(dim-1)-(i/m);

        // compute position from cell position
        int bit=1;
        for (int k=0; k<dim; k++)
        {
          if (k==ifix) continue;
          if ((i%m)&bit) coord[k] += 1;
          bit *= 2;
        }

        // do lexicographic numbering
        int index = coord[dim-1];
        for (int k=dim-2; k>=0; --k)
          if (k!=ifix)
            index = (index*(_g.cell_overlap().size(k)+1))+coord[k]; // one more
          else
            index = (index*(_g.cell_overlap().size(k)))+coord[k];

        // add size of all subsets for smaller directions
        for (int j=dim-1; j>ifix; j--)
        {
          int n=_g.cell_overlap().size(j);
          for (int l=0; l<dim; l++)
            if (l!=j) n *= _g.cell_overlap().size(l)+1;
          index += n;
        }

        return index;
      }

      DUNE_THROW(GridError, "codim " << cc << " (dim=" << dim << ") not (yet) implemented");
    }

    const GridImp * _yg;    // access to YaspGrid
    const TSI& _it;         // position in the grid level
    const YGLI& _g;         // access to grid level
  };


  // specialization for codim=dim (vertex)
  template<int dim, class GridImp>
  class YaspEntity<dim,dim,GridImp>
    : public EntityDefaultImplementation <dim,dim,GridImp,YaspEntity>
  {
    enum { dimworld = GridImp::dimensionworld };

    typedef typename GridImp::Traits::template Codim<dim>::GeometryImpl GeometryImpl;

  public:
    typedef typename GridImp::ctype ctype;

    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;

    typedef typename GridImp::template Codim<dim>::Geometry Geometry;

    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    typedef typename GridImp::template Codim<dim>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<dim>::EntitySeed EntitySeed;

    //! define the type used for persisitent indices
    typedef typename GridImp::PersistentIndexType PersistentIndexType;

    //! define type used for coordinates in grid module
    typedef typename YGrid<dim,ctype>::iTupel iTupel;

    // constructor
    YaspEntity (const GridImp* yg, const YGLI& g, const TSI& it)
      : _yg(yg), _it(it), _g(g)
    {}

    //! level of this element
    int level () const {return _g.level();}

    //! index is unique and consecutive per level
    int index () const {return _it.superindex();}

    //! globally unique, persistent index
    int globalIndex () const { return _g.cell_global().index(_it.coord()); }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as less memory as possible
     */
    EntitySeed seed () const {
      return EntitySeed(_g.level(), _it.coord());
    }

    //! geometry of this entity
    Geometry geometry () const {
      GeometryImpl _geometry(_it.position());
      return Geometry( _geometry );
    }

    //! return partition type attribute
    PartitionType partitionType () const
    {
      if (_g.vertex_interior().inside(_it.coord())) return InteriorEntity;
      if (_g.vertex_interiorborder().inside(_it.coord())) return BorderEntity;
      if (_g.vertex_overlap().inside(_it.coord())) return OverlapEntity;
      if (_g.vertex_overlapfront().inside(_it.coord())) return FrontEntity;
      return GhostEntity;
    }

  private:
    // IndexSets needs access to the private index methods
    friend class Dune::YaspLevelIndexSet<GridImp>;
    friend class Dune::YaspLeafIndexSet<GridImp>;
    friend class Dune::YaspGlobalIdSet<GridImp>;

    //! globally unique, persistent index
    PersistentIndexType persistentIndex () const
    {
      // get coordinate and size of global grid
      const iTupel& size =  _g.vertex_global().size();
      int coord[dim];

      // correction for periodic boundaries
      for (int i=0; i<dim; i++)
      {
        coord[i] = _it.coord(i);
        if (coord[i]<0) coord[i] += size[i];
        if (coord[i]>=size[i]) coord[i] -= size[i];
      }

      // determine min number of trailing zeroes
      int trailing = 1000;
      for (int i=0; i<dim; i++)
      {
        // count trailing zeros
        int zeros = 0;
        for (int j=0; j<_g.level(); j++)
          if (coord[i]&(1<<j))
            break;
          else
            zeros++;
        trailing = std::min(trailing,zeros);
      }

      // determine the level of this vertex
      int level = _g.level()-trailing;

      // encode codim
      PersistentIndexType id(dim);

      // encode level
      id = id << yaspgrid_level_bits;
      id = id+PersistentIndexType(level);

      // encode coordinates
      for (int i=dim-1; i>=0; i--)
      {
        id = id << yaspgrid_dim_bits;
        id = id+PersistentIndexType(coord[i]>>trailing);
      }

      return id;
    }

    //! consecutive, codim-wise, level-wise index
    int compressedIndex () const { return _it.superindex();}

    //! consecutive, codim-wise, level-wise index
    int compressedLeafIndex () const
    {
      if (_g.level()==_g.mg()->maxlevel())
        return _it.superindex();

      // not on leaf level, interpolate to finest grid
      int coord[dim];
      for (int i=0; i<dim; i++) coord[i] = _it.coord(i)-(_g).vertex_overlap().origin(i);

      // move coordinates up to maxlevel (multiply by 2 for each level
      for (int k=0; k<dim; k++)
        coord[k] = coord[k]*(1<<(_g.mg()->maxlevel()-_g.level()));

      // do lexicographic numbering
      int index = coord[dim-1];
      for (int k=dim-2; k>=0; --k)
        index = (index*(_g.mg()->rbegin().cell_overlap().size(k)+1))+coord[k];
      return index;
    }

  public:
    const TSI& transformingsubiterator() const { return _it; }
    const YGLI& gridlevel() const { return _g; }
    const GridImp * yaspgrid() const { return _yg; }
  protected:
    const GridImp * _yg;          // access to YaspGrid
    const TSI& _it;               // position in the grid level
    const YGLI& _g;               // access to grid level
    // temporary object
    mutable FieldVector<ctype, dim> loc; // always computed before being returned
  };

  //========================================================================
  /*!
     YaspIntersection provides data about intersection with
     neighboring codim 0 entities.
   */
  //========================================================================

  template<class GridImp>
  class YaspIntersection
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
    YaspIntersection();
    YaspIntersection& operator = (const YaspIntersection&);

    typedef typename GridImp::Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

  public:
    // types used from grids
    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef YaspSpecialEntity<0,dim,GridImp> SpecialEntity;
    typedef Dune::Intersection<const GridImp, Dune::YaspIntersectionIterator> Intersection;

    // void update() const {
    //     const_cast<YaspIntersection*>(this)->update();
    // }
    void update() const {
      if (_count == 2*_dir + _face || _count >= 2*dim)
        return;

      // cleanup old stuff
      _outside.transformingsubiterator().move(_dir,1-2*_face);   // move home
      _pos_world[_dir] = _inside.transformingsubiterator().position(_dir);

      // update face info
      _dir = _count / 2;
      _face = _count % 2;

      // move transforming iterator
      _outside.transformingsubiterator().move(_dir,-1+2*_face);

      // make up faces
      _pos_world[_dir] += (-0.5+_face)*_inside.transformingsubiterator().meshsize(_dir);
    }

    //! increment
    void increment() {
      _count += (_count < 2*dim);
    }

    //! equality
    bool equals (const YaspIntersection& other) const
    {
      return _inside.equals(other._inside) && _count == other._count;
    }

    /*! return true if neighbor ist outside the domain. Still the neighbor might
       exist in case of periodic boundary conditions, i.e. true is returned
       if the neighbor is outside the periodic unit cell
     */
    bool boundary () const
    {
#if 1
      return (_inside.transformingsubiterator().coord(_count/2) + 2*(_count%2) - 1 < _inside.gridlevel().cell_global().min(_count/2)
              ||
              _inside.transformingsubiterator().coord(_count/2) + 2*(_count%2) - 1 > _inside.gridlevel().cell_global().max(_count/2));
#else
      update();
      // The transforming iterator can be safely moved beyond the boundary.
      // So we only have to compare against the cell_global grid
      return (_outside.transformingsubiterator().coord(_dir) < _inside.gridlevel().cell_global().min(_dir)
              ||
              _outside.transformingsubiterator().coord(_dir) > _inside.gridlevel().cell_global().max(_dir));
#endif
    }

    //! return true if neighbor across intersection exists in this processor
    bool neighbor () const
    {
#if 1
      return (_inside.transformingsubiterator().coord(_count/2) + 2*(_count%2) - 1 >= _inside.gridlevel().cell_overlap().min(_count/2)
              &&
              _inside.transformingsubiterator().coord(_count/2) + 2*(_count%2) - 1 <= _inside.gridlevel().cell_overlap().max(_count/2));
#else
      update();
      return (_outside.transformingsubiterator().coord(_dir) >= _inside.gridlevel().cell_overlap().min(_dir)
              &&
              _outside.transformingsubiterator().coord(_dir) <= _inside.gridlevel().cell_overlap().max(_dir));
#endif
    }

    //! Yasp is always conform
    bool conforming () const
    {
      return true;
    }

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const
    {
      return _inside;
    }

    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const
    {
      update();
      return _outside;
    }

    //! identifier for boundary segment from macro grid
    //! (attach your boundary condition as needed)
    int boundaryId() const
    {
      if(boundary()) return indexInInside()+1;
      return 0;
    }

    //! identifier for boundary segment from macro grid
    //! (attach your boundary condition as needed)
    int boundarySegmentIndex() const
    {
      if(! boundary())
        DUNE_THROW(GridError, "called boundarySegmentIndex while boundary() == false");
      update();
      // size of local macro grid
      const FieldVector<int, dim> & size = _inside.gridlevel().mg()->begin().cell_overlap().size();
      const FieldVector<int, dim> & origin = _inside.gridlevel().mg()->begin().cell_overlap().origin();
      FieldVector<int, dim> sides;
      {
        for (int i=0; i<dim; i++)
        {
          sides[i] =
            ((_inside.gridlevel().mg()->begin().cell_overlap().origin(i)
              == _inside.gridlevel().mg()->begin().cell_global().origin(i))+
             (_inside.gridlevel().mg()->begin().cell_overlap().origin(i) +
                      _inside.gridlevel().mg()->begin().cell_overlap().size(i)
                      == _inside.gridlevel().mg()->begin().cell_global().origin(i) +
                      _inside.gridlevel().mg()->begin().cell_global().size(i)));
        }
      }
      // gobal position of the cell on macro grid
      FieldVector<int, dim> pos = _inside.transformingsubiterator().coord();
      pos /= (1<<_inside.level());
      pos -= origin;
      // compute unit-cube-face-sizes
      FieldVector<int, dim> fsize;
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

      // int rank = 0;
      // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      // std::cout << rank << "... size: " << size << " sides: " << sides
      //           << " fsize: " << fsize
      //           << " pos: " << pos << " face: " << int(_dir) << "/" << int(_face)
      //           << " index: " << index << std::endl;

      return index;
    }

    //! return unit outer normal, this should be dependent on local coordinates for higher order boundary
    FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const
    {
      return _faceInfo[_count].normal;
    }

    //! return unit outer normal, this should be dependent on local coordinates for higher order boundary
    FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const
    {
      return _faceInfo[_count].normal;
    }

    //! return unit outer normal at center of intersection geometry
    FieldVector<ctype, dimworld> centerUnitOuterNormal () const
    {
      return _faceInfo[_count].normal;
    }

    //! return unit outer normal, this should be dependent on
    //! local coordinates for higher order boundary
    //! the normal is scaled with the integration element of the intersection.
    FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const
    {
      FieldVector<ctype, dimworld> n = _faceInfo[_count].normal;
      n *= geometry().volume();
      return n;
    }

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of the element where iteration started.
     */
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( _faceInfo[_count].geom_inside );
    }

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of neighbor
     */
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( _faceInfo[_count].geom_outside );
    }

    /*! intersection of codimension 1 of this neighbor with element where iteration started.
     */
    Geometry geometry () const
    {
      update();
      GeometryImpl
      _is_global(_pos_world,_inside.transformingsubiterator().meshsize(),_dir);
      return Geometry( _is_global );
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      static const GeometryType cube(GeometryType::cube, dim-1);
      return cube;
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

    //! make intersection iterator from entity, initialize to first neighbor
    YaspIntersection (const YaspEntity<0,dim,GridImp>& myself, bool toend) :
      _inside(myself.yaspgrid(), myself.gridlevel(),
              myself.transformingsubiterator()),
      _outside(myself.yaspgrid(), myself.gridlevel(),
               myself.transformingsubiterator()),
      // initialize to first neighbor
      _count(0),
      _dir(0),
      _face(0),
      _pos_world(myself.transformingsubiterator().position())
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

      // make up faces
      _pos_world[0] -= 0.5*_inside.transformingsubiterator().meshsize(0);
    }

    //! copy constructor
    YaspIntersection (const YaspIntersection& it) :
      _inside(it._inside),
      _outside(it._outside),
      _count(it._count),
      _dir(it._dir),
      _face(it._face),
      _pos_world(it._pos_world)
    {}

    //! copy operator
    void assign (const YaspIntersection& it)
    {
      _inside = it._inside;
      _outside = it._outside;
      _count = it._count;
      _dir = it._dir;
      _face = it._face;
      _pos_world = it._pos_world;
    }

  private:
    /* EntityPointers (get automatically updated) */
    mutable YaspEntityPointer<0,GridImp> _inside;  //!< entitypointer to myself
    mutable YaspEntityPointer<0,GridImp> _outside; //!< outside entitypointer
    /* current position */
    uint8_t _count;                                //!< valid neighbor count in 0 .. 2*dim-1
    mutable uint8_t _dir;                          //!< count/2
    mutable uint8_t _face;                         //!< count%2
    /* current position */
    mutable FieldVector<ctype, dimworld> _pos_world;       //!< center of face in world coordinates

    /* static data */
    struct faceInfo
    {
      FieldVector<ctype, dimworld> normal;
      LocalGeometryImpl geom_inside;           //!< intersection in own local coordinates
      LocalGeometryImpl geom_outside;          //!< intersection in neighbors local coordinates
    };

    /* static face info */
    static const array<faceInfo, 2*GridImp::dimension> _faceInfo;

    static array<faceInfo, 2*dim> initFaceInfo()
    {
      const FieldVector<typename GridImp::ctype, GridImp::dimension> ext_local(1.0);
      array<faceInfo, 2*dim> I;
      for (uint8_t i=0; i<dim; i++)
      {
        // center of face
        FieldVector<ctype, dim> a(0.5); a[i] = 0.0;
        FieldVector<ctype, dim> b(0.5); b[i] = 1.0;
        // normal vectors
        I[2*i].normal = 0.0;
        I[2*i+1].normal = 0.0;
        I[2*i].normal[i] = -1.0;
        I[2*i+1].normal[i] = +1.0;
        // geometries
        I[2*i].geom_inside =
          LocalGeometryImpl(a, ext_local, i);
        I[2*i].geom_outside =
          LocalGeometryImpl(b, ext_local, i);
        I[2*i+1].geom_inside =
          LocalGeometryImpl(b, ext_local, i);
        I[2*i+1].geom_outside =
          LocalGeometryImpl(a, ext_local, i);
      }
      return I;
    }
  };

  template<class GridImp>
  const array<typename YaspIntersection<GridImp>::faceInfo, 2*GridImp::dimension>
  YaspIntersection<GridImp>::_faceInfo =
    YaspIntersection<GridImp>::initFaceInfo();

  //========================================================================
  /*!
     YaspIntersectionIterator enables iteration over intersection with
     neighboring codim 0 entities.
   */
  //========================================================================

  template<class GridImp>
  class YaspIntersectionIterator : public MakeableInterfaceObject< Dune::Intersection<const GridImp, Dune::YaspIntersection > >
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
    YaspIntersectionIterator();
  public:
    // types used from grids
    typedef Dune::Intersection<const GridImp, Dune::YaspIntersection> Intersection;
    typedef MakeableInterfaceObject<Intersection> MakeableIntersection;

    //! increment
    void increment()
    {
      GridImp::getRealImplementation(*this).increment();
    }

    //! equality
    bool equals (const YaspIntersectionIterator& other) const
    {
      return GridImp::getRealImplementation(*this).equals(
               GridImp::getRealImplementation(other));
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      return *this;
    }

    //! make intersection iterator from entity
    YaspIntersectionIterator (const YaspEntity<0,dim,GridImp>& myself, bool toend) :
      MakeableIntersection(YaspIntersection<GridImp>(myself, toend))
    {}

    //! copy constructor
    YaspIntersectionIterator (const YaspIntersectionIterator& it) :
      MakeableIntersection(it)
    {}

    //! assignment
    YaspIntersectionIterator & operator = (const YaspIntersectionIterator& it)
    {
      GridImp::getRealImplementation(*this).assign(
        GridImp::getRealImplementation(it));
      return *this;
    }
  };

  //========================================================================
  /*!
     YaspHierarchicIterator enables iteration over son entities of codim 0
   */
  //========================================================================

  template<class GridImp>
  class YaspHierarchicIterator :
    public YaspEntityPointer<0,GridImp>
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
  public:
    // types used from grids
    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef YaspSpecialEntity<0,dim,GridImp> SpecialEntity;

    //! define type used for coordinates in grid module
    typedef typename YGrid<dim,ctype>::iTupel iTupel;

    //! constructor
    YaspHierarchicIterator (const GridImp* yg, const YGLI& g, const TSI& it, int maxlevel) :
      YaspEntityPointer<0,GridImp>(yg,g,it)
    {
      // now iterator points to current cell
      StackElem se(this->_g);
      se.coord = this->_it.coord();
      stack.push(se);

      // determine maximum level
      _maxlevel = std::min(maxlevel,this->_g.mg()->maxlevel());

      // if maxlevel not reached then push yourself and sons
      if (this->_g.level()<_maxlevel)
      {
        push_sons();
      }

      // and make iterator point to first son if stack is not empty
      if (!stack.empty())
        pop_tos();
    }

    //! constructor
    YaspHierarchicIterator (const YaspHierarchicIterator& it) :
      YaspEntityPointer<0,GridImp>(it),
      _maxlevel(it._maxlevel), stack(it.stack)
    {}

    //! increment
    void increment ()
    {
      // sanity check: do nothing when stack is empty
      if (stack.empty()) return;

      // if maxlevel not reached then push sons
      if (this->_g.level()<_maxlevel)
        push_sons();

      // in any case pop one element
      pop_tos();
    }

    void print (std::ostream& s) const
    {
      s << "HIER: " << "level=" << this->_g.level()
        << " position=" << this->_it.coord()
        << " superindex=" << this->_it.superindex()
        << " maxlevel=" << this->_maxlevel
        << " stacksize=" << stack.size()
        << std::endl;
    }

  private:
    int _maxlevel;       //!< maximum level of elements to be processed

    struct StackElem {
      YGLI g;         // grid level of the element
      iTupel coord;   // and the coordinates
      StackElem(YGLI gg) : g(gg) {}
    };
    std::stack<StackElem> stack;    //!< stack holding elements to be processed

    // push sons of current element on the stack
    void push_sons ()
    {
      // yes, process all 1<<dim sons
      StackElem se(this->_g.finer());
      for (int i=0; i<(1<<dim); i++)
      {
        for (int k=0; k<dim; k++)
          if (i&(1<<k))
            se.coord[k] = this->_it.coord(k)*2+1;
          else
            se.coord[k] = this->_it.coord(k)*2;
        stack.push(se);
      }
    }

    // make TOS the current element
    void pop_tos ()
    {
      StackElem se = stack.top();
      stack.pop();
      this->_g = se.g;
      this->_it.reinit(this->_g.cell_overlap(),se.coord);
    }
  };

  //========================================================================
  /*!
     YaspEntitySeed describes the minimal information necessary to create a fully functional YaspEntity
   */
  //========================================================================
  template<int codim, class GridImp>
  class YaspEntitySeed
  {
    //! know your own dimension
    enum { dim=GridImp::dimension };

  public:
    //! codimension of entity pointer
    enum { codimension = codim };

    //! constructor
    YaspEntitySeed (int level, FieldVector<int, dim> coord)
      : _l(level), _c(coord)
    {}

    //! copy constructor
    YaspEntitySeed (const YaspEntitySeed& rhs)
      : _l(rhs._l), _c(rhs._c)
    {}

    int level () const { return _l; }
    const FieldVector<int, dim> & coord() const { return _c; }

  protected:
    int _l;                  // grid level
    FieldVector<int, dim> _c; // coord in the global grid
  };

  //========================================================================
  /*!
     YaspEntityPointer serves as a Reference or Pointer to a YaspGrid::Entity.
     It can also be initialized from Yasp::LevelIterator, Yasp::LeafIterator,
     Yasp::HierarchicIterator and Yasp::IntersectionIterator.

     We have specializations for codim==0 (elements) and
     codim=dim (vertices).
     The general version throws a GridError.
   */
  //========================================================================
  template<int codim, class GridImp>
  class YaspEntityPointer
  {
    //! know your own dimension
    enum { dim=GridImp::dimension };
    //! know your own dimension of world
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef YaspSpecialEntity<codim,dim,GridImp> SpecialEntity;
    typedef YaspEntityPointer<codim,GridImp> EntityPointerImp;
  protected:
    typedef YaspEntity<codim, dim, GridImp> YaspEntityImp;

  public:
    //! codimension of entity pointer
    enum { codimension = codim };

    //! constructor
    YaspEntityPointer (const GridImp * yg, const YGLI & g, const TSI & it)
      : _g(g), _it(it), _entity(yg, _g,_it)
    {
      if (codim>0 && codim<dim)
      {
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! copy constructor
    YaspEntityPointer (const YaspEntityImp& entity)
      : _g(entity.gridlevel()),
        _it(entity.transformingsubiterator()),
        _entity(entity.yaspgrid(),_g,_it)
    {
      if (codim>0 && codim<dim)
      {
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! copy constructor
    YaspEntityPointer (const YaspEntityPointer& rhs)
      : _g(rhs._g), _it(rhs._it), _entity(rhs._entity.yaspgrid(),_g,_it)
    {
      if (codim>0 && codim<dim)
      {
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! equality
    bool equals (const YaspEntityPointer& rhs) const
    {
      return (_it==rhs._it && _g == rhs._g);
    }

    //! dereferencing
    Entity& dereference() const
    {
      return _entity;
    }

    //! ask for level of entity
    int level () const {return _g.level();}

    const YaspEntityPointer&
    operator = (const YaspEntityPointer& rhs)
    {
      _g = rhs._g;
      _it = rhs._it;
      /* _entity = i._entity
       * is done implicitely, as the entity is completely
       * defined via the interator it belongs to
       */
      return *this;
    }

    const TSI& transformingsubiterator () const
    {
      return _it;
    }

    const YGLI& gridlevel () const
    {
      return _g;
    }

    TSI& transformingsubiterator ()
    {
      return _it;
    }

    YGLI& gridlevel ()
    {
      return _g;
    }

  protected:
    YGLI _g;             // access to grid level
    TSI _it;             // position in the grid level
    mutable SpecialEntity _entity; //!< virtual entity
  };

  //========================================================================
  /*!
     YaspLevelIterator enables iteration over entities of one grid level

     We have specializations for codim==0 (elements) and
     codim=dim (vertices).
     The general version throws a GridError.
   */
  //========================================================================

  template<int codim, PartitionIteratorType pitype, class GridImp>
  class YaspLevelIterator :
    public YaspEntityPointer<codim,GridImp>
  {
    //! know your own dimension
    enum { dim=GridImp::dimension };
    //! know your own dimension of world
    enum { dimworld=GridImp::dimensionworld };
    typedef typename GridImp::ctype ctype;
  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;
    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef YaspSpecialEntity<codim,dim,GridImp> SpecialEntity;

    //! constructor
    YaspLevelIterator (const GridImp * yg, const YGLI & g, const TSI & it) :
      YaspEntityPointer<codim,GridImp>(yg,g,it) {}

    //! copy constructor
    YaspLevelIterator (const YaspLevelIterator& i) :
      YaspEntityPointer<codim,GridImp>(i) {}

    //! increment
    void increment()
    {
      ++(this->_it);
    }
  };


  //========================================================================
  /*!
     \brief level-wise, non-persistent, consecutive index

   */
  //========================================================================

  template<class GridImp>
  class YaspLevelIndexSet
    : public IndexSet< GridImp, YaspLevelIndexSet< GridImp >, unsigned int >
  {
    typedef YaspLevelIndexSet< GridImp > This;
    typedef IndexSet< GridImp, This, unsigned int > Base;

  public:
    typedef typename Base::IndexType IndexType;

    using Base::subIndex;

    //! constructor stores reference to a grid and level
    YaspLevelIndexSet ( const GridImp &g, int l )
      : grid( g ),
        level( l )
    {
      // contains a single element type;
      for (int codim=0; codim<=GridImp::dimension; codim++)
        mytypes[codim].push_back(GeometryType(GeometryType::cube,GridImp::dimension-codim));
    }

    //! get index of an entity
    template<int cc>
    IndexType index (const typename GridImp::Traits::template Codim<cc>::Entity& e) const
    {
      assert( cc == 0 || cc == GridImp::dimension );
      return grid.getRealImplementation(e).compressedIndex();
    }

    //! get index of subentity of an entity
    template< int cc >
    IndexType subIndex ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      assert( cc == 0 || cc == GridImp::dimension );
      if( cc == GridImp::dimension )
        return grid.getRealImplementation(e).compressedIndex();
      else
        return grid.getRealImplementation(e).subCompressedIndex(i,codim);
    }

    //! get number of entities of given type and level (the level is known to the object)
    int size (GeometryType type) const
    {
      return grid.size( level, type );
    }

    //! return size of set for a given codim
    int size (int codim) const
    {
      return grid.size( level, codim );
    }

    //! return true if the given entity is contained in \f$E\f$.
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return e.level() == level;
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return mytypes[codim];
    }

  private:
    const GridImp& grid;
    int level;
    std::vector<GeometryType> mytypes[GridImp::dimension+1];
  };


  // Leaf Index Set

  template<class GridImp>
  class YaspLeafIndexSet
    : public IndexSet< GridImp, YaspLeafIndexSet< GridImp >, unsigned int >
  {
    typedef YaspLeafIndexSet< GridImp > This;
    typedef IndexSet< GridImp, This, unsigned int > Base;

  public:
    typedef typename Base::IndexType IndexType;

    using Base::subIndex;

    //! constructor stores reference to a grid
    explicit YaspLeafIndexSet ( const GridImp &g )
      : grid( g )
    {
      // contains a single element type;
      for (int codim=0; codim<=GridImp::dimension; codim++)
        mytypes[codim].push_back(GeometryType(GeometryType::cube,GridImp::dimension-codim));
    }

    //! get index of an entity
    template<int cc>
    IndexType index (const typename GridImp::Traits::template Codim<cc>::Entity& e) const
    {
      assert( cc == 0 || cc == GridImp::dimension );
      return grid.getRealImplementation(e).compressedIndex();
    }

    //! get index of subentity of an entity
    template< int cc >
    IndexType subIndex ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      assert( cc == 0 || cc == GridImp::dimension );
      if( cc == GridImp::dimension )
        return grid.getRealImplementation(e).compressedIndex();
      else
        return grid.getRealImplementation(e).subCompressedIndex(i,codim);
    }

    //! get number of entities of given type
    int size (GeometryType type) const
    {
      return grid.size( grid.maxLevel(), type );
    }

    //! return size of set for a given codim
    int size (int codim) const
    {
      return grid.size( grid.maxLevel(), codim );
    }

    //! return true if the given entity is contained in \f$E\f$.
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      //return e.isLeaf();
      return (e.level() == grid.maxLevel());
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return mytypes[codim];
    }

  private:
    const GridImp& grid;
    enum { ncodim = remove_const<GridImp>::type::dimension+1 };
    std::vector<GeometryType> mytypes[ncodim];
  };




  //========================================================================
  /*!
     \brief persistent, globally unique Ids

   */
  //========================================================================

  template<class GridImp>
  class YaspGlobalIdSet : public IdSet<GridImp,YaspGlobalIdSet<GridImp>,
                              typename remove_const<GridImp>::type::PersistentIndexType >
                          /*
                             We used the remove_const to extract the Type from the mutable class,
                             because the const class is not instantiated yet.
                           */
  {
    typedef YaspGlobalIdSet< GridImp > This;

  public:
    //! define the type used for persisitent indices
    typedef typename remove_const<GridImp>::type::PersistentIndexType IdType;

    using IdSet<GridImp, This, IdType>::subId;

    //! constructor stores reference to a grid
    explicit YaspGlobalIdSet ( const GridImp &g )
      : grid( g )
    {}

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return grid.getRealImplementation(e).persistentIndex();
    }

    //! get id of subentity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    IdType subId (const typename remove_const<GridImp>::type::Traits::template Codim< 0 >::Entity &e,
                  int i, unsigned int codim ) const
    {
      return grid.getRealImplementation(e).subPersistentIndex(i,codim);
    }

  private:
    const GridImp& grid;
  };


  template<int dim, int dimworld>
  struct YaspGridFamily
  {
#if HAVE_MPI
    typedef CollectiveCommunication<MPI_Comm> CCType;
#else
    typedef CollectiveCommunication<Dune::YaspGrid<dim> > CCType;
#endif

    typedef GridTraits<dim,dimworld,Dune::YaspGrid<dim>,
        YaspGeometry,YaspEntity,
        YaspEntityPointer,YaspLevelIterator,
        YaspIntersection,              // leaf  intersection
        YaspIntersection,              // level intersection
        YaspIntersectionIterator,              // leaf  intersection iter
        YaspIntersectionIterator,              // level intersection iter
        YaspHierarchicIterator,
        YaspLevelIterator,
        YaspLevelIndexSet< const YaspGrid< dim > >,
        YaspLeafIndexSet< const YaspGrid< dim > >,
        YaspGlobalIdSet<const YaspGrid<dim> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits>,
        YaspGlobalIdSet<const YaspGrid<dim> >,
        bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits>,
        CCType,
        DefaultLevelGridViewTraits, DefaultLeafGridViewTraits,
        YaspEntitySeed>
    Traits;
  };

  template<int dim, int codim>
  struct YaspCommunicateMeta {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,codim))
      {
        DUNE_THROW(GridError, "interface communication not implemented");
      }
      YaspCommunicateMeta<dim,codim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,dim> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,dim))
        g.template communicateCodim<DataHandle,dim>(data,iftype,dir,level);
      YaspCommunicateMeta<dim,dim-1>::comm(g,data,iftype,dir,level);
    }
  };

  template<int dim>
  struct YaspCommunicateMeta<dim,0> {
    template<class G, class DataHandle>
    static void comm (const G& g, DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level)
    {
      if (data.contains(dim,0))
        g.template communicateCodim<DataHandle,0>(data,iftype,dir,level);
    }
  };


  //************************************************************************
  /*!
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief Provides a distributed structured cube mesh.
     \ingroup GridImplementations

     YaspGrid stands for yet another structured parallel grid.
     It implements the dune grid interface for structured grids with codim 0
     and dim, with arbitrary overlap (including zero),
     periodic boundaries and fast implementation allowing on-the-fly computations.

     \tparam dim The dimension of the grid and its surrounding world

     \par History:
     \li started on July 31, 2004 by PB based on abstractions developed in summer 2003
   */
  template<int dim>
  class YaspGrid :
    public GridDefaultImplementation<dim,dim,yaspgrid_ctype,YaspGridFamily<dim,dim> >,
    public MultiYGrid<dim,yaspgrid_ctype>
  {
    typedef const YaspGrid<dim> GridImp;

    void init()
    {
      setsizes();
      indexsets.push_back( new YaspLevelIndexSet<const YaspGrid<dim> >(*this,0) );
      theleafindexset.push_back( new YaspLeafIndexSet<const YaspGrid<dim> >(*this) );
      theglobalidset.push_back( new YaspGlobalIdSet<const YaspGrid<dim> >(*this) );
      boundarysegmentssize();
    }

    void boundarysegmentssize()
    {
      // sizes of local macro grid
      const FieldVector<int, dim> & size = MultiYGrid<dim,ctype>::begin().cell_overlap().size();
      FieldVector<int, dim> sides;
      {
        for (int i=0; i<dim; i++)
        {
          sides[i] =
            ((MultiYGrid<dim,ctype>::begin().cell_overlap().origin(i)
              == MultiYGrid<dim,ctype>::begin().cell_global().origin(i))+
             (MultiYGrid<dim,ctype>::begin().cell_overlap().origin(i) +
                    MultiYGrid<dim,ctype>::begin().cell_overlap().size(i)
                    == MultiYGrid<dim,ctype>::begin().cell_global().origin(i) +
                    MultiYGrid<dim,ctype>::begin().cell_global().size(i)));
        }
      }
      nBSegments = 0;
      for (int k=0; k<dim; k++)
      {
        int offset = 1;
        for (int l=0; l<dim; l++)
        {
          if (l==k) continue;
          offset *= size[l];
        }
        nBSegments += sides[k]*offset;
      }
    }

  public:

    using MultiYGrid<dim,yaspgrid_ctype>::defaultLoadbalancer;

    //! define type used for coordinates in grid module
    typedef yaspgrid_ctype ctype;

    // define the persistent index type
    typedef bigunsignedint<dim*yaspgrid_dim_bits+yaspgrid_level_bits+yaspgrid_codim_bits> PersistentIndexType;

    //! the GridFamily of this grid
    typedef YaspGridFamily<dim,dim> GridFamily;
    // the Traits
    typedef typename YaspGridFamily<dim,dim>::Traits Traits;

    // need for friend declarations in entity
    typedef YaspLevelIndexSet<YaspGrid<dim> > LevelIndexSetType;
    typedef YaspLeafIndexSet<YaspGrid<dim> > LeafIndexSetType;
    typedef YaspGlobalIdSet<YaspGrid<dim> > GlobalIdSetType;

    //! maximum number of levels allowed
    enum { MAXL=64 };

    //! shorthand for base class data types
    typedef MultiYGrid<dim,ctype> YMG;
    typedef typename MultiYGrid<dim,ctype>::YGridLevelIterator YGLI;
    typedef typename SubYGrid<dim,ctype>::TransformingSubIterator TSI;
    typedef typename MultiYGrid<dim,ctype>::Intersection IS;
    typedef typename std::deque<IS>::const_iterator ISIT;

    /*! Constructor for a YaspGrid, they are all forwarded to the base class
       @param comm MPI communicator where this mesh is distributed to
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::MPIHelper::MPICommunicator comm,
              Dune::FieldVector<ctype, dim> L,
              Dune::FieldVector<int, dim> s,
              Dune::FieldVector<bool, dim> periodic, int overlap,
              const YLoadBalance<dim>* lb = defaultLoadbalancer())
#if HAVE_MPI
      : YMG(comm,L,s,periodic,overlap,lb), ccobj(comm),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#else
      : YMG(L,s,periodic,overlap,lb),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#endif
    {
      init();
    }


    /*! Constructor for a sequential YaspGrid, they are all forwarded to the base class.

       Sequential here means that the whole grid is living on one process even if your program is running
       in parallel.
       @see YaspGrid(Dune::MPIHelper::MPICommunicator, Dune::FieldVector<ctype, dim>, Dune::FieldVector<int, dim>,  Dune::FieldVector<bool, dim>, int)
       for constructing one parallel grid decomposed between the processors.
       @param L extension of the domain
       @param s number of cells on coarse mesh in each direction
       @param periodic tells if direction is periodic or not
       @param overlap size of overlap on coarsest grid (same in all directions)
       @param lb pointer to an overloaded YLoadBalance instance
     */
    YaspGrid (Dune::FieldVector<ctype, dim> L,
              Dune::FieldVector<int, dim> s,
              Dune::FieldVector<bool, dim> periodic, int overlap,
              const YLoadBalance<dim>* lb = YMG::defaultLoadbalancer())
#if HAVE_MPI
      : YMG(MPI_COMM_SELF,L,s,periodic,overlap,lb), ccobj(MPI_COMM_SELF),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#else
      : YMG(L,s,periodic,overlap,lb),
        keep_ovlp(true), adaptRefCount(0), adaptActive(false)
#endif
    {
      init();
    }

    ~YaspGrid()
    {
      deallocatePointers(indexsets);
      deallocatePointers(theleafindexset);
      deallocatePointers(theglobalidset);
    }

  private:
    // do not copy this class
    YaspGrid(const YaspGrid&);

  public:

    /*! Return maximum level defined in this grid. Levels are numbered
          0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {return MultiYGrid<dim,ctype>::maxlevel();} // delegate

    //! refine the grid refCount times. What about overlap?
    void globalRefine (int refCount)
    {
      if (refCount < -maxLevel())
        DUNE_THROW(GridError, "Only " << maxLevel() << " levels left. " <<
                   "Coarsening " << -refCount << " levels requested!");
      for (int k=refCount; k<0; k++)
      {
        MultiYGrid<dim,ctype>::coarsen();
        setsizes();
        indexsets.pop_back();
      }
      for (int k=0; k<refCount; k++)
      {
        MultiYGrid<dim,ctype>::refine(keep_ovlp);
        setsizes();
        indexsets.push_back( new YaspLevelIndexSet<const YaspGrid<dim> >(*this,maxLevel()) );
      }
    }

    /**
       \brief set options for refinement
       @param keepPhysicalOverlap [true] keep the physical size of the overlap, [false] keep the number of cells in the overlap.  Default is [true].
     */
    void refineOptions (bool keepPhysicalOverlap)
    {
      keep_ovlp = keepPhysicalOverlap;
    }

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

       \param[in] refCount Number of subdivisions that should be applied. Negative value means coarsening.
       \param[in] e        Entity to Entity that should be refined

       \return true if Entity was marked, false otherwise.

       \note
          -  On yaspgrid marking one element will mark all other elements of the level aswell
          -  If refCount is lower than refCount of a previous mark-call, nothing is changed
     */
    bool mark( int refCount, const typename Traits::template Codim<0>::Entity & e )
    {
      assert(adaptActive == false);
      if (e.level() != maxLevel()) return false;
      adaptRefCount = std::max(adaptRefCount, refCount);
      return true;
    }

    /** \brief returns adaptation mark for given entity

       \param[in] e   Entity for which adaptation mark should be determined

       \return int adaptation mark, here the default value 0 is returned
     */
    int getMark ( const typename Traits::template Codim<0>::Entity &e ) const
    {
      return ( e.level() == maxLevel() ) ? adaptRefCount : 0;
    }

    //! map adapt to global refine
    bool adapt ()
    {
      globalRefine(adaptRefCount);
      return (adaptRefCount > 0);
    }

    //! returns true, if the grid will be coarsened
    bool preAdapt ()
    {
      adaptActive = true;
      adaptRefCount = comm().max(adaptRefCount);
      return (adaptRefCount < 0);
    }

    //! clean up some markers
    void postAdapt()
    {
      adaptActive = false;
      adaptRefCount = 0;
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,pitype>(level);
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      return levelend<cd,pitype>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
      return levelbegin<cd,All_Partition>(level);
    }

    //! version without second template parameter for convenience
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
      return levelend<cd,All_Partition>(level);
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LeafIterator leafend () const
    {
      return levelend<cd,pitype>(maxLevel());
    }

    //! return LeafIterator which points to the first entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
      return levelbegin<cd,All_Partition>(maxLevel());
    }

    //! return LeafIterator which points behind the last entity in maxLevel
    template<int cd>
    typename Traits::template Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
      return levelend<cd,All_Partition>(maxLevel());
    }

    // \brief obtain EntityPointer from EntitySeed. */
    template <typename Seed>
    typename Traits::template Codim<Seed::codimension>::EntityPointer
    entityPointer(const Seed& seed) const
    {
      static const int codim = Seed::codimension;
      YGLI g = MultiYGrid<dim,ctype>::begin(seed.level());
      switch (codim)
      {
      case 0 :
        return YaspEntityPointer<codim,GridImp>(this,g,
                                                TSI(g.cell_overlap(), seed.coord()));
      case dim :
        return YaspEntityPointer<codim,GridImp>(this,g,
                                                TSI(g.vertex_overlap(), seed.coord()));
      default :
        DUNE_THROW(GridError, "YaspEntityPointer: codim not implemented");
      }
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int level, int codim) const
    {
      YGLI g = MultiYGrid<dim,ctype>::begin(level);
      return g.overlap();
    }

    //! return size (= distance in graph) of overlap region
    int overlapSize (int codim) const
    {
      YGLI g = MultiYGrid<dim,ctype>::begin(maxLevel());
      return g.overlap();
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int level, int codim) const
    {
      return 0;
    }

    //! return size (= distance in graph) of ghost region
    int ghostSize (int codim) const
    {
      return 0;
    }

    //! number of entities per level and codim in this process
    int size (int level, int codim) const
    {
      return sizes[level][codim];
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return sizes[maxLevel()][codim];
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return (type.isCube()) ? sizes[level][dim-type.dim()] : 0;
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return size(maxLevel(),type);
    }

    //! \brief returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const
    {
      return nBSegments;
    }

    /*! The new communication interface

       communicate objects for all codims on a given level
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,level);
    }

    /*! The new communication interface

       communicate objects for all codims on the leaf grid
     */
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp,DataType> & data, InterfaceType iftype, CommunicationDirection dir) const
    {
      YaspCommunicateMeta<dim,dim>::comm(*this,data,iftype,dir,this->maxLevel());
    }

    /*! The new communication interface

       communicate objects for one codim
     */
    template<class DataHandle, int codim>
    void communicateCodim (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      // check input
      if (!data.contains(dim,codim)) return; // should have been checked outside

      // data types
      typedef typename DataHandle::DataType DataType;

      // access to grid level
      YGLI g = MultiYGrid<dim,ctype>::begin(level);

      // find send/recv lists or throw error
      const std::deque<IS>* sendlist=0;
      const std::deque<IS>* recvlist=0;
      if (codim==0) // the elements
      {
        if (iftype==InteriorBorder_InteriorBorder_Interface)
          return; // there is nothing to do in this case
        if (iftype==InteriorBorder_All_Interface)
        {
          sendlist = &g.send_cell_interior_overlap();
          recvlist = &g.recv_cell_overlap_interior();
        }
        if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface || iftype==All_All_Interface)
        {
          sendlist = &g.send_cell_overlap_overlap();
          recvlist = &g.recv_cell_overlap_overlap();
        }
      }
      if (codim==dim) // the vertices
      {
        if (iftype==InteriorBorder_InteriorBorder_Interface)
        {
          sendlist = &g.send_vertex_interiorborder_interiorborder();
          recvlist = &g.recv_vertex_interiorborder_interiorborder();
        }

        if (iftype==InteriorBorder_All_Interface)
        {
          sendlist = &g.send_vertex_interiorborder_overlapfront();
          recvlist = &g.recv_vertex_overlapfront_interiorborder();
        }
        if (iftype==Overlap_OverlapFront_Interface || iftype==Overlap_All_Interface)
        {
          sendlist = &g.send_vertex_overlap_overlapfront();
          recvlist = &g.recv_vertex_overlapfront_overlap();
        }
        if (iftype==All_All_Interface)
        {
          sendlist = &g.send_vertex_overlapfront_overlapfront();
          recvlist = &g.recv_vertex_overlapfront_overlapfront();
        }
      }

      // change communication direction?
      if (dir==BackwardCommunication)
        std::swap(sendlist,recvlist);

      int cnt;

      // Size computation (requires communication if variable size)
      std::vector<int> send_size(sendlist->size(),-1);    // map rank to total number of objects (of type DataType) to be sent
      std::vector<int> recv_size(recvlist->size(),-1);    // map rank to total number of objects (of type DataType) to be recvd
      std::vector<size_t*> send_sizes(sendlist->size(),static_cast<size_t*>(0)); // map rank to array giving number of objects per entity to be sent
      std::vector<size_t*> recv_sizes(recvlist->size(),static_cast<size_t*>(0)); // map rank to array giving number of objects per entity to be recvd
      if (data.fixedsize(dim,codim))
      {
        // fixed size: just take a dummy entity, size can be computed without communication
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          send_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          recv_size[cnt] = is->grid.totalsize() * data.size(*it);
          cnt++;
        }
      }
      else
      {
        // variable size case: sender side determines the size
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          // allocate send buffer for sizes per entitiy
          size_t *buf = new size_t[is->grid.totalsize()];
          send_sizes[cnt] = buf;

          // loop over entities and ask for size
          int i=0; size_t n=0;
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
          {
            buf[i] = data.size(*it);
            n += buf[i];
            i++;
          }

          // now we know the size for this rank
          send_size[cnt] = n;

          // hand over send request to torus class
          MultiYGrid<dim,ctype>::torus().send(is->rank,buf,is->grid.totalsize()*sizeof(size_t));
          cnt++;
        }

        // allocate recv buffers for sizes and store receive request
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          // allocate recv buffer
          size_t *buf = new size_t[is->grid.totalsize()];
          recv_sizes[cnt] = buf;

          // hand over recv request to torus class
          MultiYGrid<dim,ctype>::torus().recv(is->rank,buf,is->grid.totalsize()*sizeof(size_t));
          cnt++;
        }

        // exchange all size buffers now
        MultiYGrid<dim,ctype>::torus().exchange();

        // release send size buffers
        cnt=0;
        for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
        {
          delete[] send_sizes[cnt];
          send_sizes[cnt] = 0;
          cnt++;
        }

        // process receive size buffers
        cnt=0;
        for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
        {
          // get recv buffer
          size_t *buf = recv_sizes[cnt];

          // compute total size
          size_t n=0;
          for (int i=0; i<is->grid.totalsize(); ++i)
            n += buf[i];

          // ... and store it
          recv_size[cnt] = n;
          ++cnt;
        }
      }


      // allocate & fill the send buffers & store send request
      std::vector<DataType*> sends(sendlist->size(), static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        //      std::cout << "[" << this->comm().rank() << "] "
        //                << " send " << " dest=" << is->rank
        //                << " size=" << send_size[cnt]
        //                << std::endl;

        // allocate send buffer
        DataType *buf = new DataType[send_size[cnt]];

        // remember send buffer
        sends[cnt] = buf;

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // fill send buffer; iterate over cells in intersection
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
        typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
        tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
        for ( ; it!=tsubend; ++it)
          data.gather(mb,*it);

        // hand over send request to torus class
        MultiYGrid<dim,ctype>::torus().send(is->rank,buf,send_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // allocate recv buffers and store receive request
      std::vector<DataType*> recvs(recvlist->size(),static_cast<DataType*>(0)); // store pointers to send buffers
      cnt=0;
      for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        //      std::cout << "[" << this->comm().rank() << "] "
        //                << " recv " << "  src=" << is->rank
        //                << " size=" << recv_size[cnt]
        //                << std::endl;

        // allocate recv buffer
        DataType *buf = new DataType[recv_size[cnt]];

        // remember recv buffer
        recvs[cnt] = buf;

        // hand over recv request to torus class
        MultiYGrid<dim,ctype>::torus().recv(is->rank,buf,recv_size[cnt]*sizeof(DataType));
        cnt++;
      }

      // exchange all buffers now
      MultiYGrid<dim,ctype>::torus().exchange();

      // release send buffers
      cnt=0;
      for (ISIT is=sendlist->begin(); is!=sendlist->end(); ++is)
      {
        delete[] sends[cnt];
        sends[cnt] = 0;
        cnt++;
      }

      // process receive buffers and delete them
      cnt=0;
      for (ISIT is=recvlist->begin(); is!=recvlist->end(); ++is)
      {
        // get recv buffer
        DataType *buf = recvs[cnt];

        // make a message buffer
        MessageBuffer<DataType> mb(buf);

        // copy data from receive buffer; iterate over cells in intersection
        if (data.fixedsize(dim,codim))
        {
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          size_t n=data.size(*it);
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
            data.scatter(mb,*it,n);
        }
        else
        {
          int i=0;
          size_t *sbuf = recv_sizes[cnt];
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          it(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubbegin()));
          typename Traits::template Codim<codim>::template Partition<All_Partition>::LevelIterator
          tsubend(YaspLevelIterator<codim,All_Partition,GridImp>(this,g,is->grid.tsubend()));
          for ( ; it!=tsubend; ++it)
            data.scatter(mb,*it,sbuf[i++]);
          delete[] sbuf;
        }

        // delete buffer
        delete[] buf; // hier krachts !
        cnt++;
      }
    }

    // The new index sets from DDM 11.07.2005
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return *(theglobalidset[0]);
    }

    const typename Traits::LocalIdSet& localIdSet() const
    {
      return *(theglobalidset[0]);
    }

    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      return *(indexsets[level]);
    }

    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return *(theleafindexset[0]);
    }

#if HAVE_MPI
    /*! @brief return a collective communication object
     */
    const CollectiveCommunication<MPI_Comm>& comm () const
    {
      return ccobj;
    }
#else
    /*! @brief return a collective communication object
     */
    const CollectiveCommunication<YaspGrid>& comm () const
    {
      return ccobj;
    }
#endif

    YaspIntersectionIterator<const YaspGrid<dim> >&
    getRealIntersectionIterator(typename Traits::LevelIntersectionIterator& it)
    {
      return this->getRealImplementation(it);
    }

    const YaspIntersectionIterator<const YaspGrid<dim> >&
    getRealIntersectionIterator(const typename Traits::LevelIntersectionIterator& it) const
    {
      return this->getRealImplementation(it);
    }

  private:

#if HAVE_MPI
    CollectiveCommunication<MPI_Comm> ccobj;
#else
    CollectiveCommunication<YaspGrid> ccobj;
#endif

    std::vector<YaspLevelIndexSet<const YaspGrid<dim> >*> indexsets;
    std::vector<YaspLeafIndexSet<const YaspGrid<dim> >*> theleafindexset;
    std::vector<YaspGlobalIdSet<const YaspGrid<dim> >*> theglobalidset;
    int nBSegments;

    // Index classes need access to the real entity
    friend class Dune::YaspLevelIndexSet<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspLeafIndexSet<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspGlobalIdSet<const Dune::YaspGrid<dim> >;

    friend class Dune::YaspIntersectionIterator<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspIntersection<const Dune::YaspGrid<dim> >;
    friend class Dune::YaspEntity<0, dim, const Dune::YaspGrid<dim> >;

    template<class T>
    void deallocatePointers(T& container)
    {
      typedef typename T::iterator Iterator;

      for(Iterator entry=container.begin(); entry != container.end(); ++entry)
        delete (*entry);
    }

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    template<class DT>
    class MessageBuffer {
    public:
      // Constructor
      MessageBuffer (DT *p)
      {
        a=p;
        i=0;
        j=0;
      }

      // write data to message buffer, acts like a stream !
      template<class Y>
      void write (const Y& data)
      {
        dune_static_assert(( is_same<DT,Y>::value ), "DataType missmatch");
        a[i++] = data;
      }

      // read data from message buffer, acts like a stream !
      template<class Y>
      void read (Y& data) const
      {
        dune_static_assert(( is_same<DT,Y>::value ), "DataType missmatch");
        data = a[j++];
      }

    private:
      DT *a;
      int i;
      mutable int j;
    };

    void setsizes ()
    {
      for (YGLI g=MultiYGrid<dim,ctype>::begin(); g!=MultiYGrid<dim,ctype>::end(); ++g)
      {
        // codim 0 (elements)
        sizes[g.level()][0] = 1;
        for (int i=0; i<dim; ++i)
          sizes[g.level()][0] *= g.cell_overlap().size(i);

        // codim 1 (faces)
        if (dim>1)
        {
          sizes[g.level()][1] = 0;
          for (int i=0; i<dim; ++i)
          {
            int s=g.cell_overlap().size(i)+1;
            for (int j=0; j<dim; ++j)
              if (j!=i)
                s *= g.cell_overlap().size(j);
            sizes[g.level()][1] += s;
          }
        }

        // codim dim-1 (edges)
        if (dim>2)
        {
          sizes[g.level()][dim-1] = 0;
          for (int i=0; i<dim; ++i)
          {
            int s=g.cell_overlap().size(i);
            for (int j=0; j<dim; ++j)
              if (j!=i)
                s *= g.cell_overlap().size(j)+1;
            sizes[g.level()][dim-1] += s;
          }
        }

        // codim dim (vertices)
        sizes[g.level()][dim] = 1;
        for (int i=0; i<dim; ++i)
          sizes[g.level()][dim] *= g.vertex_overlapfront().size(i);
      }
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelbegin (int level) const
    {
      dune_static_assert( cd == dim || cd == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      YGLI g = MultiYGrid<dim,ctype>::begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      if (pitype==Ghost_Partition)
        return levelend <cd, pitype> (level);
      if (cd==0)   // the elements
      {
        if (pitype<=InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_interior().tsubbegin());
        if (pitype<=All_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_overlap().tsubbegin());
      }
      if (cd==dim)   // the vertices
      {
        if (pitype==Interior_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interior().tsubbegin());
        if (pitype==InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interiorborder().tsubbegin());
        if (pitype==Overlap_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlap().tsubbegin());
        if (pitype<=All_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlapfront().tsubbegin());
      }
      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    YaspLevelIterator<cd,pitype,GridImp> levelend (int level) const
    {
      dune_static_assert( cd == dim || cd == 0 ,
                          "YaspGrid only supports Entities with codim=dim and codim=0");
      YGLI g = MultiYGrid<dim,ctype>::begin(level);
      if (level<0 || level>maxLevel()) DUNE_THROW(RangeError, "level out of range");
      if (cd==0)   // the elements
      {
        if (pitype<=InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_interior().tsubend());
        if (pitype<=All_Partition || pitype == Ghost_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.cell_overlap().tsubend());
      }
      if (cd==dim)   // the vertices
      {
        if (pitype==Interior_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interior().tsubend());
        if (pitype==InteriorBorder_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_interiorborder().tsubend());
        if (pitype==Overlap_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlap().tsubend());
        if (pitype<=All_Partition || pitype == Ghost_Partition)
          return YaspLevelIterator<cd,pitype,GridImp>(this,g,g.vertex_overlapfront().tsubend());
      }
      DUNE_THROW(GridError, "YaspLevelIterator with this codim or partition type not implemented");
    }

    int sizes[MAXL][dim+1]; // total number of entities per level and codim
    bool keep_ovlp;
    int adaptRefCount;
    bool adaptActive;
  };

  namespace Capabilities
  {

    /** \struct hasEntity
       \ingroup YaspGrid
     */

    /** \struct hasBackupRestoreFacilities
       \ingroup YaspGrid
     */

    /** \brief YaspGrid has only one geometry type for codim 0 entities
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasSingleGeometryType< YaspGrid<dim> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GenericGeometry :: CubeTopology< dim > :: type :: id ;
    };

    /** \brief YaspGrid is a Cartesian grid
        \ingroup YaspGrid
     */
    template<int dim>
    struct isCartesian< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has codim=0 entities (elements)
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasEntity< YaspGrid<dim>, 0 >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid has codim=dim entities (vertices)
       \ingroup YaspGrid
     */
    template<int dim>
    struct hasEntity< YaspGrid<dim>, dim >
    {
      static const bool v = true;
    };

    template< int dim >
    struct canCommunicate< YaspGrid< dim >, 0 >
    {
      static const bool v = true;
    };

    template< int dim >
    struct canCommunicate< YaspGrid< dim >, dim >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is parallel
       \ingroup YaspGrid
     */
    template<int dim>
    struct isParallel< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is levelwise conforming
       \ingroup YaspGrid
     */
    template<int dim>
    struct isLevelwiseConforming< YaspGrid<dim> >
    {
      static const bool v = true;
    };

    /** \brief YaspGrid is leafwise conforming
       \ingroup YaspGrid
     */
    template<int dim>
    struct isLeafwiseConforming< YaspGrid<dim> >
    {
      static const bool v = true;
    };

  }

} // end namespace


#endif