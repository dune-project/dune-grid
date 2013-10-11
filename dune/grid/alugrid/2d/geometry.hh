// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGEOMETRY_HH
#define DUNE_ALU2DGRIDGEOMETRY_HH

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/alugrid/2d/alu2dinclude.hh>
#include <dune/grid/alugrid/3d/mappings.hh>
#include <dune/grid/alugrid/common/memory.hh>

namespace Dune
{

  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int cdim, class GridImp>
  class ALU2dGridGeometry;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;


  template< int mydim, int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl;

  template <int ncorners, class Mapping>
  class MyALU2dGridGeometryImplBase
  {
  private:
    // prohibited due to reference counting
    MyALU2dGridGeometryImplBase( const MyALU2dGridGeometryImplBase& );

  protected:
    //! number of corners
    static const int corners_ = ncorners ;

    //! the type of the mapping
    typedef Mapping MappingType;

    typedef typename MappingType::ctype ctype;

    typedef typename MappingType::map_t map_t;
    typedef typename MappingType::world_t world_t;

    typedef typename MappingType::matrix_t matrix_t;
    typedef typename MappingType::inv_t inv_t;

    // get my dimension
    enum { mydim = ncorners < 3 ? ncorners-1 : 2 };
    typedef ReferenceElement< ctype, mydim > ReferenceElementType ;

    //! the mapping
    MappingType mapping_;

    //! reference element
    const ReferenceElementType& referenceElement_ ;

    //! volume of element
    double volume_ ;

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
        referenceElement_( ReferenceElements< ctype, mydim >::general( type ) ),
        volume_( 1.0 )
    {
      reset();
    }

    //! reset status and reference count
    void reset()
    {
      // reset reference counter
      refCount_ = 1;
      // reset status
      valid_ = false ;
    }

    //! increase reference count
    void operator ++ () { ++ refCount_; }

    //! decrease reference count
    void operator -- () { assert( refCount_ > 0 ); --refCount_; }

    //! return true if object has no references anymore
    bool operator ! () const { return refCount_ == 0; }

    //! return true if there exists more then on reference
    bool stillUsed () const { return refCount_ > 1 ; }

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

  // geometry implementation for vertices
  template< int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl< 0, cdim, eltype >
    : public MyALU2dGridGeometryImplBase< 1, LinearMapping< cdim, 0 > >
  {
    typedef MyALU2dGridGeometryImplBase< 1, LinearMapping< cdim, 0 > > BaseType;
  protected:
    using BaseType :: mapping_ ;
    using BaseType :: valid_ ;

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
    template< class Vector >
    void update ( const Vector &p0 )
    {
      mapping_.buildMapping( p0 );
      valid_ = true ;
    }

    // return determinante of mapping
    ctype det ( const map_t &m ) const {  return 1.0;  }
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

    // update geometry in father coordinates
    template< class Geo, class LocalGeo >
    void updateLocal ( const Geo &geo, const LocalGeo &localGeo )
    {
      assert( localGeo.corners() == corners_ );
      // compute the local coordinates in father refelem
      FieldMatrix< alu2d_ctype, corners_, cdim > coord;
      for( int i = 0; i < corners_; ++i )
      {
        // calculate coordinate
        coord[ i ] = geo.local( localGeo.corner( i ) );
        // to avoid rounding errors
        for( int j = 0; j < cdim; ++j )
          coord[ i ][ j ] = (coord[ i ][ j ] < 1e-14 ? 0 : coord[ i ][ j ]);
      }
      mapping_.buildMapping( coord[ 0 ], coord[ 1 ] );
      volume_ = mapping_.det( map_t(0.25) );
      valid_ = true ;
    }

    // update geometry coordinates
    template< class Vector >
    void update ( const Vector &p0, const Vector &p1, const double volume )
    {
      mapping_.buildMapping( p0, p1 );
      volume_ = volume ;
      valid_ = true ;
    }

    // return determinante of mapping
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

    // update geometry in father coordinates
    template< class Geo, class LocalGeo >
    void updateLocal ( const Geo &geo, const LocalGeo &localGeo )
    {
      assert( localGeo.corners() == corners_ );
      // compute the local coordinates in father refelem
      FieldMatrix< alu2d_ctype, corners_, cdim > coord;
      for( int i = 0; i < corners_; ++i )
      {
        // calculate coordinate
        coord[ i ] = geo.local( localGeo.corner( i ) );
        // to avoid rounding errors
        for( int j = 0; j < cdim; ++j )
          coord[ i ][ j ] = (coord[ i ][ j ] < 1e-14 ? 0 : coord[ i ][ j ]);
      }
      mapping_.buildMapping( coord[ 0 ], coord[ 1 ], coord[ 2 ] );
      // set volume to a fraction of the reference volume which is 0.5
      volume_ = referenceElement_.volume() * ( localGeo.volume() / geo.volume() );
      assert( (volume_ > 0.0) && (volume_ < referenceElement_.volume() ) );
      valid_  = true ;
    }

    template< class HElement >
    void update ( const HElement &item )
    {
      mapping_.buildMapping( item.getVertex( 0 )->coord(), item.getVertex( 1 )->coord(),
                             item.getVertex( 2 )->coord() );
      volume_ = item.area();
      valid_ = true ;
    }

    // return determinante of mapping
    ctype det ( const map_t &m ) const
    {
      return 2.0 * volume_ ;
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

    // update geometry in father coordinates
    template< class Geo, class LocalGeo >
    void updateLocal ( const Geo &geo, const LocalGeo &localGeo )
    {
      assert( localGeo.corners() == corners_ );
      // compute the local coordinates in father refelem
      FieldMatrix< alu2d_ctype, corners_, cdim > coord;
      for( int i = 0; i < corners_; ++i )
      {
        // calculate coordinate
        coord[ i ] = geo.local( localGeo.corner( i ) );
        // to avoid rounding errors
        for( int j = 0; j < cdim; ++j )
          coord[ i ][ j ] = (coord[ i ][ j ] < 1e-14 ? 0 : coord[ i ][ j ]);
      }
      mapping_.buildMapping( coord[ 0 ], coord[ 1 ], coord[ 2 ], coord[ 3 ] );
      volume_ = referenceElement_.volume() * ( localGeo.volume() / geo.volume() );
      assert( (volume_ > 0.0) && (volume_ < referenceElement_.volume() ) );
      valid_ = true ;
    }

    template< class HElement >
    void update ( const HElement &item )
    {
      mapping_.buildMapping( item.getVertex( 0 )->coord(), item.getVertex( 1 )->coord(),
                             item.getVertex( 3 )->coord(), item.getVertex( 2 )->coord() );
      volume_ = item.area();
      valid_ = true ;
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
        simplexReferenceElement_( ReferenceElements< ctype, 2 >::general( type( 3 ) ) ),
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

    // update geometry in father coordinates
    template< class Geo, class LocalGeo >
    void updateLocal ( const Geo &geo, const LocalGeo &localGeo )
    {
      const int corners = localGeo.corners();

      // compute the local coordinates in father refelem
      FieldMatrix< alu2d_ctype, 4, cdim > coord;
      for( int i = 0; i < corners; ++i )
      {
        // calculate coordinate
        coord[ i ] = geo.local( localGeo.corner( i ) );
        // to avoid rounding errors
        for( int j = 0; j < cdim; ++j )
          coord[ i ][ j ] = (coord[ i ][ j ] < 1e-14 ? 0 : coord[ i ][ j ]);
      }

      updateMapping( corners );
      if( corners == 3 )
      {
        linearMapping().buildMapping( coord[ 0 ], coord[ 1 ], coord[ 2 ] );
        volume_ = simplexReferenceElement_.volume() * ( localGeo.volume() / geo.volume() );
      }
      else
      {
        bilinearMapping().buildMapping( coord[ 0 ], coord[ 1 ], coord[ 2 ], coord[ 3 ] );
        volume_ = referenceElement_.volume() * ( localGeo.volume() / geo.volume() );
      }

      assert( (volume_ > 0.0) && (volume_ < 1.0) );
      valid_ = true ;
    }

    template< class HElement >
    void update ( const HElement &item )
    {
      const int corners = item.numvertices();
      updateMapping( corners );
      if( corners == 3 )
        linearMapping().buildMapping( item.getVertex( 0 )->coord(), item.getVertex( 1 )->coord(),
                                      item.getVertex( 2 )->coord() );
      else
        bilinearMapping().buildMapping( item.getVertex( 0 )->coord(), item.getVertex( 1 )->coord(),
                                        item.getVertex( 3 )->coord(), item.getVertex( 2 )->coord() );

      volume_ = item.area();
      valid_  = true ;
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


  //**********************************************************************
  //
  // --ALU2dGridGeometry
  // --Geometry
  //**********************************************************************
  /*!
     Defines the geometry part of a mesh entity. Works for all dimensions, element types and dimensions
     of world. Provides reference element and mapping between local and global coordinates.
     The element may have different implementations because the mapping can be
     done more efficient for structured meshes than for unstructured meshes.

     dim: An element is a polygonal in a hyperplane of dimension dim. 0 <= dim <= 2 is typically
     dim=0 is a point.

     dimworld: Each corner is a point with dimworld coordinates.
   */

  //! ALU2dGridGeometry
  //! Empty definition, needs to be specialized for element type
  template< int mydim, int cdim, class GridImp >
  class ALU2dGridGeometry
    : public GeometryDefaultImplementation< mydim, cdim, GridImp, ALU2dGridGeometry >
  {
    static const ALU2DSPACE ElementType eltype = GridImp::elementType;

    //! type of our Geometry interface
    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    //! type of our Geometry implementation
    typedef ALU2dGridGeometry<mydim,cdim,GridImp> GeometryImp;
    //! know dimension of barycentric coordinates
    enum { dimbary=mydim+1};

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
    ALU2dGridGeometry();

    //! copy constructor copying pointer and increasing reference counter
    ALU2dGridGeometry( const ALU2dGridGeometry& );

    //! destructor releasing object
    ~ALU2dGridGeometry() ;

    //! assigment operator
    ALU2dGridGeometry& operator =( const ALU2dGridGeometry& );

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
    // method for elements
    bool buildGeom(const HElementType & item);
    // method for edges
    bool buildGeom(const HElementType & item, const int aluFace);
    // method for vertices
    bool buildGeom(const VertexType & item, const int );

    //! build geometry for intersectionSelfLocal and
    //! intersectionNeighborLocal
    template <class GeometryType, class LocalGeomType >
    bool buildLocalGeom(const GeometryType & geo , const LocalGeomType & lg);

    //! build local geometry given local face number
    bool buildLocalGeometry(const int faceNumber, const int twist,const int coorns);

    //! return non-const reference to coord vecs
    GlobalCoordinate& getCoordVec (int i);

    //! print internal data
    void print (std::ostream& ss) const;

    //! build geometry with local coords of child in reference element
    inline bool buildGeomInFather(const Geometry &fatherGeom,
                                  const Geometry & myGeom );

    // returns true if geometry information is valid
    inline bool valid() const { return geoImpl().valid(); }

    // invalidate geometry implementation
    void invalidate () ;

  protected:
    // return reference coordinates of the alu triangle
    static std::pair< FieldMatrix< alu2d_ctype, 4, 2 >, FieldVector< alu2d_ctype, 4 > >
    calculateReferenceCoords ( const int corners );
  protected:
    //! assign pointer
    void assign( const ALU2dGridGeometry& other );
    //! remove pointer object
    void removeObj();
    //! get a new pointer object
    void getObject();

    // type of object provider
    typedef ALUMemoryProvider< GeometryImplType > GeometryProviderType ;

    //! return storage provider for geometry objects
    static GeometryProviderType& geoProvider()
    {
#ifdef USE_SMP_PARALLEL
      typedef ALUGridObjectFactory< GridImp >  GridObjectFactoryType;
      static std::vector< GeometryProviderType > storage( GridObjectFactoryType :: maxThreads() );
      return storage[ GridObjectFactoryType :: threadNumber () ];
#else
      static GeometryProviderType storage;
      return storage;
#endif
    }

    // return reference to geometry implementation
    GeometryImplType& geoImpl() const
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
      //! Whether to store by reference.
      static const bool v = false;
    };
  }

} // end namespace Dune

#include "geometry_imp.cc"

#endif
