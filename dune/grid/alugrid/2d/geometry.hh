// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_2D_GEOMETRY_HH
#define DUNE_ALUGRID_2D_GEOMETRY_HH

#include <dune/common/array.hh>
#include <dune/common/misc.hh>
#include <dune/common/nullptr.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/alugrid/2d/alu2dinclude.hh>
#include <dune/grid/alugrid/3d/mappings.hh>
#include <dune/grid/alugrid/common/memory.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int cdim, class GridImp>
  class ALU2dGridGeometry;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;



  // MyALU2dGridGeometryImplBase
  // ---------------------------

  template< int ncorners, class Mapping >
  class MyALU2dGridGeometryImplBase
  {
  private:
    // prohibited due to reference counting
    MyALU2dGridGeometryImplBase( const MyALU2dGridGeometryImplBase& );

  protected:
    //! number of corners
    static const int corners_ = ncorners;

    //! the type of the mapping
    typedef Mapping MappingType;

    typedef typename MappingType::ctype ctype;

    typedef typename MappingType::map_t map_t;
    typedef typename MappingType::world_t world_t;

    typedef typename MappingType::matrix_t matrix_t;
    typedef typename MappingType::inv_t inv_t;

    // get my dimension
    enum { mydim = ncorners < 3 ? ncorners-1 : 2 };
    typedef Dune::ReferenceElement< ctype, mydim > ReferenceElementType;

    //! the mapping
    MappingType mapping_;

    //! reference element
    const ReferenceElementType& referenceElement_;

    //! volume of element
    ctype volume_;

    //! the reference counter
    mutable unsigned int refCount_;

    //! valid flag, true if mapping was built
    bool valid_ ;

    const MappingType& mapping() const
    {
      assert( valid() );
      return mapping_;
    }

  public:
    //! default constructor
    MyALU2dGridGeometryImplBase( const GeometryType type )
      : mapping_(),
        referenceElement_( Dune::ReferenceElements< ctype, mydim >::general( type ) ),
        volume_( 1.0 )
    {
      reset();
    }

    //! reset status and reference count
    void reset ()
    {
      // reset reference counter
      refCount_ = 1;
      // reset status
      valid_ = false;
    }

    //! increase reference count
    void operator++ () { ++refCount_; }

    //! decrease reference count
    void operator-- () { assert( refCount_ > 0 ); --refCount_; }

    //! return true if object has no references anymore
    bool operator ! () const { return (refCount_ == 0); }

    // set status to invalid
    void invalidate () { valid_ = false ; }

    // return true if geometry is valid
    bool valid () const { return valid_; }

    // return volume
    double volume() const { return volume_; }

    // return true if geometry is affine
    bool affine() const
    {
      assert( valid() );
      return mapping().affine();
    }

    // return number of corners
    int corners () const { return corners_ ; }

    // return coordinate of the i-th corner
    world_t corner( int i ) const
    {
      world_t coordinate;
      map2world( referenceElement_.position( i, mydim ), coordinate );
      return coordinate;
    }

    // map from reference to global coordinates
    void map2world ( const map_t &m, world_t &w ) const
    {
      return mapping().map2world( m, w );
    }

    // map from global to reference coordinates
    void world2map ( const world_t &w, map_t &m ) const
    {
      return mapping().world2map( w, m );
    }

    // return jacobian transposed
    const matrix_t &jacobianTransposed ( const map_t &m ) const
    {
      return mapping().jacobianTransposed( m );
    }

    // return jacobian inverse transposed
    const inv_t &jacobianInverseTransposed ( const map_t &m ) const
    {
      return mapping().jacobianInverseTransposed( m );
    }

    // return determinante of mapping
    ctype det ( const map_t &m ) const {  return mapping().det( m );  }
  };



  // MyALU2dGridGeometryImpl
  // -----------------------

  template< int mydim, int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl;

  // geometry implementation for vertices
  template< int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl< 0, cdim, eltype >
    : public MyALU2dGridGeometryImplBase< 1, LinearMapping< cdim, 0 > >
  {
    typedef MyALU2dGridGeometryImplBase< 1, LinearMapping< cdim, 0 > > BaseType;

  protected:
    using BaseType::mapping_;
    using BaseType::valid_;
    using BaseType::volume_;

  public:
    using BaseType :: valid ;
    using BaseType :: invalidate ;
    using BaseType :: corners ;

    typedef typename BaseType :: ctype ctype ;
    typedef typename BaseType :: map_t map_t ;

    // default constructor
    MyALU2dGridGeometryImpl () : BaseType( type() ) {}

    GeometryType type () const
    {
      return GeometryType(
               (eltype == ALU2DSPACE triangle ?
                GenericGeometry :: SimplexTopology< 0 > :: type :: id :
                GenericGeometry :: CubeTopology   < 0 > :: type :: id),
               0 );
    }

    // update geometry coordinates
    template< class CoordVector >
    void update ( const CoordVector &coordVector )
    {
      mapping_.buildMapping( coordVector[ 0 ] );
      volume_ = coordVector.volume();
      valid_ = true;
    }

    // return determinante of mapping
    ctype det ( const map_t &m ) const { return volume_; }
  };

  // geometry implementation for lines
  template< int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl< 1, cdim, eltype >
    : public MyALU2dGridGeometryImplBase< 2, LinearMapping< cdim, 1 > >
  {
    typedef MyALU2dGridGeometryImplBase< 2, LinearMapping< cdim, 1 > > BaseType;
  protected:
    using BaseType :: mapping_ ;
    using BaseType :: valid_ ;
    using BaseType :: volume_ ;
    using BaseType :: corners_ ;

  public:
    using BaseType :: valid ;
    using BaseType :: invalidate ;
    using BaseType :: corners ;

    typedef typename BaseType :: ctype ctype ;
    typedef typename BaseType :: map_t map_t ;

    // default constructor
    MyALU2dGridGeometryImpl () : BaseType( type() ) {}

    GeometryType type () const
    {
      return GeometryType(
               (eltype == ALU2DSPACE triangle ?
                GenericGeometry :: SimplexTopology< 1 > :: type :: id :
                GenericGeometry :: CubeTopology   < 1 > :: type :: id),
               1 );
    }

    // update geometry coordinates
    template< class Vector >
    void update ( const Vector &p0, const Vector &p1, ctype volume )
    {
      mapping_.buildMapping( p0, p1 );
      volume_ = volume;
      valid_ = true;
    }

    template< class CoordVector >
    void update ( const CoordVector &coordVector )
    {
      mapping_.buildMapping( coordVector[ 0 ], coordVector[ 1 ] );
      volume_ = coordVector.volume();
      valid_ = true;
    }

    // return determinant of mapping
    ctype det ( const map_t &m ) const { return volume_; }
  };

  // geometry implementation for triangles
  template< int cdim >
  class MyALU2dGridGeometryImpl< 2, cdim, ALU2DSPACE triangle >
    : public MyALU2dGridGeometryImplBase< 3, LinearMapping< cdim, 2 > >
  {
    typedef MyALU2dGridGeometryImplBase< 3, LinearMapping< cdim, 2 > > BaseType;
  protected:
    using BaseType :: mapping_ ;
    using BaseType :: valid_ ;
    using BaseType :: volume_ ;
    using BaseType :: corners_ ;
    using BaseType :: referenceElement_;

  public:
    using BaseType :: valid ;
    using BaseType :: invalidate ;
    using BaseType :: corners ;

    typedef typename BaseType :: ctype ctype ;
    typedef typename BaseType :: map_t map_t ;

    // default constructor
    MyALU2dGridGeometryImpl () : BaseType( type() )
    {
      // if this assertion fails change factor in method det below
      assert( std::abs( referenceElement_.volume() - 0.5 ) < 1e-10 );
    }

    GeometryType type () const
    {
      return GeometryType( GenericGeometry :: SimplexTopology< 2 > :: type :: id , 2 );
    }

    template< class CoordVector >
    void update ( const CoordVector &coordVector )
    {
      mapping_.buildMapping( coordVector[ 0 ], coordVector[ 1 ], coordVector[ 2 ] );
      volume_ = coordVector.volume();
      valid_ = true;
    }

    // return determinante of mapping
    ctype det ( const map_t &m ) const
    {
      return 2.0 * volume_;
    }
  };

  // geometry implementation for quadrilaterals
  template< int cdim >
  class MyALU2dGridGeometryImpl< 2, cdim, ALU2DSPACE quadrilateral >
    : public MyALU2dGridGeometryImplBase< 4, BilinearMapping< cdim > >
  {
    typedef MyALU2dGridGeometryImplBase< 4, BilinearMapping< cdim > > BaseType;
  protected:
    using BaseType :: mapping_ ;
    using BaseType :: valid_ ;
    using BaseType :: volume_ ;
    using BaseType :: corners_ ;
    using BaseType :: referenceElement_;

  public:
    using BaseType :: valid ;
    using BaseType :: invalidate ;
    using BaseType :: corners ;

    // default constructor
    MyALU2dGridGeometryImpl () : BaseType( type() ) {}

    GeometryType type () const
    {
      return GeometryType( GenericGeometry :: CubeTopology< 2 > :: type :: id, 2 ) ;
    }

    template< class CoordVector >
    void update ( const CoordVector &coordVector )
    {
      mapping_.buildMapping( coordVector[ 0 ], coordVector[ 1 ], coordVector[ 2 ], coordVector[ 3 ] );
      volume_ = coordVector.volume();
      valid_ = true;
    }
  };

  // geometry implementation for triangles
  template< int cdim >
  class MyALU2dGridGeometryImpl< 2, cdim, ALU2DSPACE mixed >
    : public MyALU2dGridGeometryImplBase< 4, BilinearMapping< cdim > >
  {
    typedef MyALU2dGridGeometryImplBase< 4, BilinearMapping< cdim > > BaseType;
  protected:
    typedef typename BaseType :: MappingType BilinearMappingType;
    typedef Dune :: LinearMapping< cdim, 2 >  LinearMappingType;

    using BaseType :: mapping_ ;
    using BaseType :: volume_ ;
    using BaseType :: valid_ ;
    using BaseType :: referenceElement_ ;

    typedef typename LinearMappingType::ctype ctype;

    typedef typename LinearMappingType::map_t map_t;
    typedef typename LinearMappingType::world_t world_t;

    typedef typename LinearMappingType::matrix_t matrix_t;
    typedef typename LinearMappingType::inv_t inv_t;

    typedef typename BaseType :: ReferenceElementType ReferenceElementType;

    //! reference element
    const ReferenceElementType& simplexReferenceElement_ ;

    // for the mixed geometry the number of corners can vary
    int myCorners_;

  public:
    using BaseType :: valid ;
    using BaseType :: invalidate ;

    // default constructor
    MyALU2dGridGeometryImpl ()
      : BaseType( type( 4 ) ),
        simplexReferenceElement_( Dune::ReferenceElements< ctype, 2 >::general( type( 3 ) ) ),
        myCorners_( 0 )
    {
      // make sure that bilinear mapping reserves more memory, othersize change
      assert( sizeof( BilinearMappingType ) >= sizeof( LinearMappingType ) );
    }

  public:
    // return true if current mapping is affine
    bool affine () const
    {
      return (corners() == 3 ? linearMapping().affine() : bilinearMapping().affine());
    }

    // return current number of corners
    int corners () const { return myCorners_; }

    // return coordinate of the i-th corner
    world_t corner( int i ) const
    {
      world_t coordinate;
      if( corners() == 3 )
        linearMapping().map2world( simplexReferenceElement_.position( i, 2 ), coordinate );
      else
        bilinearMapping().map2world( referenceElement_.position( i, 2 ), coordinate );
      return coordinate;
    }

    // return current type of geometry
    GeometryType type () const { return type( corners() ); }

  protected:
    // return type of geometry depending on number of corners
    GeometryType type ( const int corners ) const
    {
      return GeometryType( (corners == 3 ?
                            GenericGeometry :: SimplexTopology< 2 > :: type :: id :
                            GenericGeometry :: CubeTopology   < 2 > :: type :: id), 2);
    }

  public:
    void map2world ( const map_t &m, world_t &w ) const
    {
      if( corners() == 3 )
        linearMapping().map2world( m, w );
      else
        bilinearMapping().map2world( m, w );
    }

    void world2map ( const world_t &w, map_t &m ) const
    {
      if( corners() == 3 )
        linearMapping().world2map( w, m );
      else
        bilinearMapping().world2map( w, m );
    }

    const matrix_t &jacobianTransposed ( const map_t &m ) const
    {
      return (corners() == 3 ? linearMapping().jacobianTransposed( m ) : bilinearMapping().jacobianTransposed( m ));
    }

    const inv_t &jacobianInverseTransposed ( const map_t &m ) const
    {
      return (corners() == 3 ? linearMapping().jacobianInverseTransposed( m ) : bilinearMapping().jacobianInverseTransposed( m ));
    }

    ctype det ( const map_t &m ) const
    {
      return (corners() == 3 ? linearMapping().det( m ) : bilinearMapping().det( m ));
    }

    template< class CoordVector >
    void update ( const CoordVector &coordVector )
    {
      const int corners = coordVector.corners();
      updateMapping( corners );
      if( corners == 3 )
        linearMapping().buildMapping( coordVector[ 0 ], coordVector[ 1 ], coordVector[ 2 ] );
      else
        bilinearMapping().buildMapping( coordVector[ 0 ], coordVector[ 1 ], coordVector[ 2 ], coordVector[ 3 ] );
      volume_ = coordVector.volume();
      valid_ = true;
    }

  private:
    MyALU2dGridGeometryImpl &operator= ( const MyALU2dGridGeometryImpl &other );

    const LinearMappingType &linearMapping () const
    {
      assert( valid() );
      return static_cast< const LinearMappingType * >( &mapping_ );
    }

    LinearMappingType &linearMapping ()
    {
      assert( valid() );
      return static_cast< LinearMappingType * >( &mapping_ );
    }

    const BilinearMappingType &bilinearMapping () const
    {
      assert( valid() );
      return static_cast< const BilinearMappingType * >( &mapping_ );
    }

    BilinearMappingType &bilinearMapping ()
    {
      assert( valid() );
      return static_cast< BilinearMappingType * >( &mapping_ );
    }

    void updateMapping ( const int corners )
    {
      assert( (corners == 3) || (corners == 4) );
      if( corners != myCorners_ )
      {
        destroyMapping();
        corners = myCorners_;
        if( corners == 3 )
          new( &mapping_ )LinearMappingType;
        else
          new( &mapping_ )BilinearMappingType;
      }
    }

    void destroyMapping ()
    {
      if( corners() == 3 )
        linearMapping().~LinearMappingType();
      else if( corners() == 4 )
        bilinearMapping().~BilinearMappingType();
    }
  };



  // ALU2dGridGeometry
  // -----------------

  template< int mydim, int cdim, class GridImp >
  class ALU2dGridGeometry
    : public GeometryDefaultImplementation< mydim, cdim, GridImp, ALU2dGridGeometry >
  {
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    //! type of our Geometry interface
    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    //! type of our Geometry implementation
    typedef ALU2dGridGeometry<mydim,cdim,GridImp> GeometryImp;

    typedef typename ALU2dImplTraits< GridImp::dimensionworld, eltype >::HElementType HElementType ;
    typedef typename ALU2dImplInterface< 0, GridImp::dimensionworld, eltype >::Type VertexType;

    // type of specialized geometry implementation
    typedef MyALU2dGridGeometryImpl< mydim, cdim, eltype > GeometryImplType;

  public:
    typedef FieldVector< alu2d_ctype, cdim > GlobalCoordinate;
    typedef FieldVector< alu2d_ctype, mydim > LocalCoordinate;

  public:
    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    ALU2dGridGeometry ();

    //! copy constructor copying pointer and increasing reference counter
    ALU2dGridGeometry ( const ALU2dGridGeometry & );

    //! destructor releasing object
    ~ALU2dGridGeometry ();

    //! assigment operator
    const ALU2dGridGeometry &operator= ( const ALU2dGridGeometry & );

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    const GeometryType type () const { return geoImpl().type(); }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const { return geoImpl().corners(); }

    //! access to coordinates of corners. Index is the number of the corner
    GlobalCoordinate corner ( int i ) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    GlobalCoordinate global ( const LocalCoordinate& local ) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    LocalCoordinate local (const GlobalCoordinate& global) const;

    //! A(l) , see grid.hh
    alu2d_ctype integrationElement (const LocalCoordinate& local) const;

    //! return volume of geometry
    alu2d_ctype volume () const;

    //! return true if geometry has affine mapping
    bool affine() const { return geoImpl().affine(); }

    //! jacobian inverse transposed
    const FieldMatrix<alu2d_ctype,cdim,mydim>& jacobianInverseTransposed (const LocalCoordinate& local) const;

    //! jacobian transposed
    const FieldMatrix<alu2d_ctype,mydim,cdim>& jacobianTransposed (const LocalCoordinate& local) const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for out of given ALU2dGridElement
    template< class CoordVector >
    bool buildGeom ( const CoordVector &coordVector );

    // returns true if geometry information is valid
    inline bool valid() const { return geoImpl().valid(); }

    // invalidate geometry implementation
    void invalidate () ;

  protected:
    //! get a new pointer object
    void getObject();

    // type of object provider
    typedef ALUMemoryProvider< GeometryImplType > GeometryProviderType ;

    //! return storage provider for geometry objects
    static GeometryProviderType& geoProvider()
    {
      static GeometryProviderType storage;
      return storage;
    }

    // return reference to geometry implementation
    GeometryImplType &geoImpl () const
    {
      assert( geoImpl_ );
      return *geoImpl_;
    }

    // implementation of the coordinates and mapping
    GeometryImplType* geoImpl_;
  };

  namespace FacadeOptions
  {
    //! geometry can be stored as an object
    template< int mydim, int cdim, class GridImp >
    struct StoreGeometryReference< mydim, cdim, GridImp, ALU2dGridGeometry >
    {
      static const bool v = false;
    };
  }



  // Implementation of ALU2dGridGeometry
  // -----------------------------------

  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >::ALU2dGridGeometry ()
  {
    getObject();
  }


  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >::ALU2dGridGeometry ( const ALU2dGridGeometry &other )
    : geoImpl_( other.geoImpl_ )
  {
    ++geoImpl();
  }


  template< int mydim, int cdim, class GridImp >
  inline ALU2dGridGeometry< mydim, cdim, GridImp >::~ALU2dGridGeometry ()
  {
    --geoImpl();
    if( !geoImpl() )
      geoProvider().freeObject( geoImpl_ );
  }


  template< int mydim, int cdim, class GridImp >
  inline const ALU2dGridGeometry< mydim, cdim, GridImp > &
  ALU2dGridGeometry< mydim, cdim, GridImp >::operator= ( const ALU2dGridGeometry &other )
  {
    ++other.geoImpl();
    --geoImpl();
    if( !geoImpl() )
      geoProvider().freeObject( geoImpl_ );
    geoImpl_ = other.geoImpl_;
    return *this;
  }


  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::getObject ()
  {
    geoImpl_ = geoProvider().getEmptyObject();
    geoImpl().reset();
  }


  template< int mydim, int cdim, class GridImp >
  inline void ALU2dGridGeometry< mydim, cdim, GridImp >::invalidate ()
  {
    --geoImpl();
    if( !geoImpl() )
    {
      ++geoImpl();
      geoImpl().invalidate();
    }
    else
      getObject();
  }


  //! access to coordinates of corners. Index is the number of the corner
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::corner ( int i ) const
  {
    return geoImpl().corner( i );
  }


  //! maps a local coordinate within reference element to
  //! global coordinate in element
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::GlobalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::global ( const LocalCoordinate &local ) const
  {
    GlobalCoordinate global;
    geoImpl().map2world( local, global );
    return global;
  }


  //! maps a global coordinate within the element to a
  //! local coordinate in its reference element
  template< int mydim, int cdim, class GridImp >
  inline typename ALU2dGridGeometry< mydim, cdim, GridImp >::LocalCoordinate
  ALU2dGridGeometry< mydim, cdim, GridImp >::local ( const GlobalCoordinate &global ) const
  {
    if( mydim == 0 )
      return LocalCoordinate( 1 );

    LocalCoordinate local;
    geoImpl().world2map( global, local );
    return local;
  }


  template< int mydim, int cdim, class GridImp >
  inline alu2d_ctype
  ALU2dGridGeometry< mydim, cdim, GridImp >::integrationElement ( const LocalCoordinate &local ) const
  {
    assert( geoImpl().valid() );
    return geoImpl().det( local );
  }


  template< int mydim, int cdim, class GridImp >
  inline alu2d_ctype ALU2dGridGeometry< mydim, cdim, GridImp >::volume () const
  {
    assert( geoImpl().valid() );
    return geoImpl().volume();
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix<alu2d_ctype,mydim,cdim>&
  ALU2dGridGeometry< mydim, cdim, GridImp>::jacobianTransposed ( const LocalCoordinate &local ) const
  {
    return geoImpl().jacobianTransposed( local );
  }


  template< int mydim, int cdim, class GridImp >
  inline const FieldMatrix<alu2d_ctype,cdim,mydim>&
  ALU2dGridGeometry< mydim, cdim, GridImp >::jacobianInverseTransposed ( const LocalCoordinate &local ) const
  {
    return geoImpl().jacobianInverseTransposed( local );
  }


  //! built Geometry for triangles
  template< int mydim, int cdim, class GridImp >
  template< class CoordVector >
  inline bool
  ALU2dGridGeometry< mydim, cdim, GridImp >::buildGeom ( const CoordVector &coordVector )
  {
    geoImpl().update( coordVector );
    return true;
  }



  // ALU2dGridCoordVectorBase
  // ------------------------

  template< class Grid >
  struct ALU2dGridCoordVectorBase
  {
    static const int dimensionworld = remove_const< Grid >::type::dimensionworld;
    static const ALU2DSPACE ElementType elementType = remove_const< Grid >::type::elementType;

    typedef typename remove_const< Grid >::type::ctype ctype;
    typedef FieldVector< ctype, dimensionworld > GlobalCoordinate;

  protected:
    template< class Vector >
    static GlobalCoordinate convert ( const Vector &x )
    {
      GlobalCoordinate y;
      for( int i = 0; i < dimensionworld; ++i )
        y[ i ] = x[ i ];
      return y;
    }
  };



  // ALU2dGridCoordVector
  // --------------------

  template< int mydim, class Grid >
  class ALU2dGridCoordVector;

  template< class Grid >
  class ALU2dGridCoordVector< 2, Grid >
    : public ALU2dGridCoordVectorBase< Grid >
  {
    typedef ALU2dGridCoordVectorBase< Grid > Base;

  public:
    typedef typename ALU2dImplInterface< 2, Base::dimensionworld, Base::elementType >::Type Item;

    ALU2dGridCoordVector ( const Item &item )
      : item_( item )
    {}

    int corners () const
    {
      switch( Base::elementType )
      {
      case ALU2DSPACE triangle :
        return 3;

      case ALU2DSPACE quadrilateral :
        return 4;

      case ALU2DSPACE mixed :
        return item_.numvertices();
      }
    }

    typename Base::GlobalCoordinate operator[] ( int i ) const
    {
      assert( (i >= 0) && (i < corners()) );
      // for cubes: exchange corners 2 and 3
      const int j = (corners() == 4 ? i^(i >> 1) : i);
      return Base::convert( item_.getVertex( j )->coord() );
    }

    typename Base::ctype volume () const { return item_.area(); }

  private:
    const Item &item_;
  };

  template< class Grid >
  class ALU2dGridCoordVector< 1, Grid >
    : public ALU2dGridCoordVectorBase< Grid >
  {
    typedef ALU2dGridCoordVectorBase< Grid > Base;

  public:
    typedef typename ALU2dImplInterface< 1, Base::dimensionworld, Base::elementType >::Type Item;

    ALU2dGridCoordVector ( const Item &item, int aluFace )
      : item_( item ),
        aluFace_( aluFace ),
        numFaces_( item_.numfaces() )
    {}

    int corners () const { return 2; }

    typename Base::GlobalCoordinate operator[] ( int i ) const
    {
      assert( (i >= 0) && (i < corners()) );

      // for triangles face 1 is twisted and for quatrilaterals faces 1 and 2
      const int twist = (numFaces_ == 3 ? (aluFace_ % 2) : (aluFace_ >> 1)^(aluFace_ & 1));
      const int j = (i + twist) & 1;

      return Base::convert( item_.getVertex( (aluFace_ + 1 + j) % numFaces_ )->coord() );
    }

    typename Base::ctype volume () const { return item_.sidelength( aluFace_ ); }

  private:
    const Item &item_;
    int aluFace_, numFaces_;
  };

  template< class Grid >
  class ALU2dGridCoordVector< 0, Grid >
    : public ALU2dGridCoordVectorBase< Grid >
  {
    typedef ALU2dGridCoordVectorBase< Grid > Base;

  public:
    typedef typename ALU2dImplInterface< 0, Base::dimensionworld, Base::elementType >::Type Item;

    ALU2dGridCoordVector ( const Item &item, int )
      : item_( item )
    {}

    int corners () const { return 1; }

    typename Base::GlobalCoordinate operator[] ( int i ) const
    {
      assert( (i >= 0) && (i < corners()) );
      return Base::convert( item_.coord() );
    }

    typename Base::ctype volume () const { return 1; }

  private:
    const Item &item_;
  };

} // namespace Dune

#endif // #ifndef DUNE_ALUGRID_2D_GEOMETRY_HH
