// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDGEOMETRY_HH
#define DUNE_ALU3DGRIDGEOMETRY_HH

// System includes

// Dune includes
#include <dune/common/power.hh>
#include <dune/grid/common/grid.hh>

// Local includes
#include "alu3dinclude.hh"
#include "topology.hh"
#include "mappings.hh"
#include <dune/grid/alugrid/common/objectfactory.hh>

namespace Dune
{

  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU3dGridEntity;
  template<int cd, class GridImp >
  class ALU3dGridEntityPointer;
  template<int mydim, int coorddim, class GridImp>
  class ALU3dGridGeometry;
  template< ALU3dGridElementType, class >
  class ALU3dGrid;
  class BilinearSurfaceMapping;
  class TrilinearMapping;

  template< class GridImp >
  class ALU3dGridIntersectionIterator;

  template <int cdim>
  class MyALUGridGeometryImplementation
  {
  public:
    typedef FieldVector<alu3d_ctype, cdim> CoordinateVectorType;

    static const signed char invalid      = -1; // means geometry is not meaningful
    static const signed char updated      =  0; // means the point values have been set
    static const signed char buildmapping =  1; // means updated and mapping was build

    template <int dim, int corners, class Mapping>
    class GeometryImplBase
    {
    private:
      // prohibited due to reference counting
      GeometryImplBase( const GeometryImplBase& );

    protected:
      //! number of corners
      static const int corners_ = corners ;

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners , cdim>  CoordinateMatrixType;

      template <int dummy, int dimused>
      struct CoordTypeExtractorType
      {
        typedef CoordinateMatrixType Type;
      };

      template <int dummy>
      struct CoordTypeExtractorType< dummy, 3 >
      {
        typedef CoordinateMatrixType* Type;
      };

      typedef typename CoordTypeExtractorType< 0, dim > :: Type CoordinateStorageType ;

      //! the type of the mapping
      typedef Mapping MappingType;

      //! to coordinates
      CoordinateStorageType coord_ ;

      //! the mapping
      MappingType map_;

      //! volume of element
      double volume_ ;

      //! the reference counter
      mutable unsigned int refCount_;

      //! the status (see different status above)
      signed char status_ ;
    public:
      //! default constructor
      GeometryImplBase()
        : coord_( 0 ),
          map_(),
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
        status_   = invalid ;
      }

      //! increase reference count
      void operator ++ () { ++ refCount_; }

      //! decrease reference count
      void operator -- () { assert( refCount_ > 0 ); --refCount_; }

      //! return true if object has no references anymore
      bool operator ! () const { return refCount_ == 0; }

      //! return true if there exists more then on reference
      bool stillUsed () const { return refCount_ > 1 ; }

      // copy coordinate vector from field vector or alu3d_ctype[cdim]
      template <class CoordPtrType>
      static inline void copy(const CoordPtrType& p,
                              CoordinateVectorType& c)
      {
        assert( cdim == 3 );
        c[0] = p[0];
        c[1] = p[1];
        c[2] = p[2];
      }

      template <class CoordPtrType>
      void update(const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType& ) const
      {
        DUNE_THROW(InvalidStateException,"This method should not be called!");
      }

      template <class CoordPtrType>
      void update(const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType& ) const
      {
        DUNE_THROW(InvalidStateException,"This method should not be called!");
      }

      template <class CoordPtrType>
      void update(const CoordPtrType&,
                  const CoordPtrType&,
                  const CoordPtrType& ) const
      {
        DUNE_THROW(InvalidStateException,"This method should not be called!");
      }

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }

      // set volume
      void setVolume( const double volume ) { volume_ = volume ; }

      // return volume
      double volume() const { return volume_; }
    };

    //! general type of geometry implementation
    template <int dummy, int dim,
        ALU3dGridElementType eltype> class GeometryImpl;
  public:
    // geometry implementation for edges and vertices
    template <int dummy, int dim, ALU3dGridElementType eltype>
    class GeometryImpl : public GeometryImplBase< dim, dim+1, LinearMapping<cdim, dim> >
    {
      typedef GeometryImplBase< dim, dim+1, LinearMapping<cdim, dim> > BaseType;

      using BaseType :: corners_ ;
      using BaseType :: copy ;
      using BaseType :: coord_ ;
      using BaseType :: map_ ;
      using BaseType :: status_ ;

      typedef typename BaseType :: MappingType MappingType ;
    public:
      using BaseType :: update ;
      using BaseType :: valid ;

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( valid() );
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      inline MappingType& mapping()
      {
        assert( valid() );
        if( status_ == buildmapping ) return map_;

        map_.buildMapping( coord_[0] );
        status_ = buildmapping ;
        return map_;
      }

      // update vertex
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0)
      {
        assert( corners_ == 1 );
        copy( p0, coord_[0] );
        // we need to update the mapping
        status_ = updated ;
      }
    };

    // geometry implementation for edges and vertices
    template <int dummy, ALU3dGridElementType eltype>
    class GeometryImpl<dummy,1,eltype>
      : public GeometryImplBase< 1, 2, LinearMapping<cdim, 1> >
    {
      enum { dim = 1 };
      typedef GeometryImplBase< dim, dim+1, LinearMapping<cdim, dim> > BaseType;

      using BaseType :: corners_ ;
      using BaseType :: copy ;
      using BaseType :: coord_ ;
      using BaseType :: map_ ;
      using BaseType :: status_ ;

      typedef typename BaseType :: MappingType MappingType;
    public:
      using BaseType :: update ;
      using BaseType :: valid ;

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( valid() );
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      inline MappingType& mapping()
      {
        assert( valid() );
        if( status_ == buildmapping ) return map_;

        map_.buildMapping( coord_[0], coord_[1] );
        status_ = buildmapping ;
        return map_;
      }

      // update edge
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1)
      {
        assert( corners_ == 2 );
        copy( p0, coord_[0] );
        copy( p1, coord_[1] );
        status_ = updated;
      }
    };

    // geom impl for simplex faces (triangles)
    template <int dummy>
    class GeometryImpl<dummy, 2, tetra>
      : public GeometryImplBase< 2, 3, LinearMapping<cdim, 2> >
    {
      // dim = 2, corners = 3
      typedef GeometryImplBase< 2, 3, LinearMapping<cdim, 2> > BaseType;

      using BaseType :: corners_ ;
      using BaseType :: copy ;
      using BaseType :: coord_ ;
      using BaseType :: map_ ;
      using BaseType :: status_ ;

      typedef typename BaseType :: MappingType MappingType ;
    public:
      using BaseType :: update ;
      using BaseType :: valid ;

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( valid() );
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      // update geometry coordinates
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1,
                         const CoordPtrType& p2)
      {
        copy(p0, coord_[0] );
        copy(p1, coord_[1] );
        copy(p2, coord_[2] );
        status_ = updated;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        assert( valid() );
        if( status_ == buildmapping ) return map_;

        map_.buildMapping( coord_[0], coord_[1], coord_[2] );
        status_ = buildmapping ;
        return map_;
      }
    };

    ///////////////////////////////////////////////////////////////
    //
    //  hexa specializations
    //
    ///////////////////////////////////////////////////////////////

    // geom impl for cube faces (quadrilaterals)
    template <int dummy>
    class GeometryImpl<dummy, 2, hexa>
      : public GeometryImplBase< 2, 4, BilinearSurfaceMapping >
    {
      // dim = 2, corners = 4
      typedef GeometryImplBase< 2, 4, BilinearSurfaceMapping > BaseType;

      using BaseType :: corners_ ;
      using BaseType :: copy ;
      using BaseType :: coord_ ;
      using BaseType :: map_ ;
      using BaseType :: status_ ;

      typedef typename BaseType :: MappingType MappingType ;
    public:
      using BaseType :: update ;
      using BaseType :: valid ;

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( valid() );
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      // update geometry coordinates
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1,
                         const CoordPtrType& p2,
                         const CoordPtrType& p3)
      {
        copy(p0, coord_[0] );
        copy(p1, coord_[1] );
        copy(p2, coord_[2] );
        copy(p3, coord_[3] );
        status_ = updated;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        assert( valid() );
        if( status_ == buildmapping ) return map_;

        map_.buildMapping( coord_[0], coord_[1], coord_[2], coord_[3] );
        status_ = buildmapping ;
        return map_;
      }
    };

    // geometry impl for hexahedrons
    template <int dummy>
    class GeometryImpl<dummy,3, hexa>
      : public GeometryImplBase< 3, 8, TrilinearMapping >
    {
      // dim = 3, corners = 8
      typedef GeometryImplBase< 3, 8, TrilinearMapping > BaseType;

      using BaseType :: corners_ ;
      using BaseType :: copy ;
      using BaseType :: coord_ ;
      using BaseType :: map_ ;
      using BaseType :: status_ ;

      typedef typename BaseType :: MappingType MappingType ;
      typedef typename BaseType :: CoordinateMatrixType CoordinateMatrixType;

      typedef alu3d_ctype CoordPtrType[cdim];

      // coordinate pointer vector
      const alu3d_ctype* coordPtr_[ corners_ ];
    public:
      using BaseType :: update ;
      using BaseType :: valid ;

      //! constructor creating geo impl
      GeometryImpl() : BaseType()
      {
        // set initialize coord pointers
        for( int i=0; i<corners_; ++i )
          coordPtr_[ i ] = 0;
      }

      // desctructor
      ~GeometryImpl()
      {
        if( coord_ ) delete coord_;
      }

      const alu3d_ctype* point( const int i ) const
      {
        assert( valid() );
        assert( i>=0 && i<corners_ );
        assert( coordPtr_[i] );
        return coordPtr_[ i ];
      }

      // return coordinates
      inline CoordinateVectorType operator [] (const int i) const
      {
        CoordinateVectorType coord ;
        copy( point( i ), coord );
        return coord ;
      }

      // update geometry coordinates
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1,
                         const CoordPtrType& p2,
                         const CoordPtrType& p3,
                         const CoordPtrType& p4,
                         const CoordPtrType& p5,
                         const CoordPtrType& p6,
                         const CoordPtrType& p7)
      {
        coordPtr_[0] = &p0[ 0 ];
        coordPtr_[1] = &p1[ 0 ];
        coordPtr_[2] = &p2[ 0 ];
        coordPtr_[3] = &p3[ 0 ];
        coordPtr_[4] = &p4[ 0 ];
        coordPtr_[5] = &p5[ 0 ];
        coordPtr_[6] = &p6[ 0 ];
        coordPtr_[7] = &p7[ 0 ];
        status_ = updated;
      }

      // update geometry in father coordinates
      template <class GeometryImp>
      inline void updateInFather(const GeometryImp &fatherGeom ,
                                 const GeometryImp &myGeom)
      {
        if( coord_ == 0 )
        {
          coord_ = new CoordinateMatrixType();
        }

        CoordinateMatrixType& coord = *coord_;
        // compute the local coordinates in father refelem
        for(int i=0; i < myGeom.corners() ; ++i)
        {
          // calculate coordinate
          coord[i] = fatherGeom.local( myGeom.corner( i ) );

          // set pointer
          coordPtr_[i] = (&(coord[i][0]));

          // to avoid rounding errors
          for(int j=0; j<cdim; ++j)
          {
            if ( coord[i][j] < 1e-16) coord[i][j] = 0.0;
          }
        }

        status_ = updated ;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        assert( valid() );
        if( status_ == buildmapping ) return map_;

        map_.buildMapping( point( 0 ), point( 1 ), point( 2 ), point( 3 ),
                           point( 4 ), point( 5 ), point( 6 ), point( 7 ) );

        status_ = buildmapping;
        return map_;
      }

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
    };


    // geometry impl for hexahedrons
    template <int dummy>
    class GeometryImpl<dummy,3, tetra>
      : public GeometryImplBase< 3, 4, LinearMapping<cdim, cdim> >
    {
      // dim = 3, corners = 8
      typedef GeometryImplBase< 3, 4, LinearMapping<cdim, cdim> > BaseType;

      using BaseType :: corners_ ;
      using BaseType :: copy ;
      using BaseType :: coord_ ;
      using BaseType :: map_ ;
      using BaseType :: status_ ;

      typedef typename BaseType :: MappingType MappingType ;
      typedef typename BaseType :: CoordinateMatrixType CoordinateMatrixType;

      typedef alu3d_ctype CoordPtrType[cdim];

      // coordinate pointer vector
      const alu3d_ctype* coordPtr_[ corners_ ];
    public:
      using BaseType :: update ;
      using BaseType :: valid ;

      // default constructor
      GeometryImpl() : BaseType()
      {
        // set initialize coord pointers
        for( int i=0; i<corners_; ++i )
          coordPtr_[ i ] = 0;
      }

      // destructor
      ~GeometryImpl()
      {
        if( coord_ ) delete coord_;
      }

      const alu3d_ctype* point( const int i ) const
      {
        assert( valid() );
        assert( i>=0 && i<corners_ );
        assert( coordPtr_[ i ] );
        return coordPtr_[ i ];
      }

      // return coordinate vector
      inline CoordinateVectorType operator [] (const int i) const
      {
        CoordinateVectorType coord ;
        copy( point( i ), coord );
        return coord ;
      }

      // update geometry coordinates
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1,
                         const CoordPtrType& p2,
                         const CoordPtrType& p3)
      {
        coordPtr_[0] = &p0[ 0 ];
        coordPtr_[1] = &p1[ 0 ];
        coordPtr_[2] = &p2[ 0 ];
        coordPtr_[3] = &p3[ 0 ];
        status_ = updated;
      }

      // update geometry in father coordinates
      template <class GeometryImp>
      inline void updateInFather(const GeometryImp &fatherGeom ,
                                 const GeometryImp & myGeom)
      {
        if( coord_ == 0 )
        {
          coord_ = new CoordinateMatrixType();
        }

        CoordinateMatrixType& coord = *coord_;
        // compute the local coordinates in father refelem
        for(int i=0; i < myGeom.corners() ; ++i)
        {
          // calculate coordinate
          coord[i] = fatherGeom.local( myGeom.corner( i ) );

          // set pointer
          coordPtr_[i] = (&(coord[i][0]));

          // to avoid rounding errors
          for(int j=0; j<cdim; ++j)
          {
            if ( coord[i][j] < 1e-16) coord[i][j] = 0.0;
          }
        }

        status_ = updated;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        assert( valid() );
        if( status_ == buildmapping ) return map_;

        map_.buildMapping( point( 0 ), point( 1 ), point( 2 ), point( 3 ) );

        status_ = buildmapping;
        return map_;
      }
    };
  }; // end of class ALUGridGeometryImplementation

  template <int mydim, int cdim, class GridImp>
  class ALU3dGridGeometry :
    public GeometryDefaultImplementation<mydim, cdim, GridImp, ALU3dGridGeometry>
  {
    static const ALU3dGridElementType elementType = GridImp::elementType;

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridIntersectionIterator<GridImp>;

    typedef typename ALU3dImplTraits< elementType, Comm >::IMPLElementType IMPLElementType;
    typedef typename ALU3dImplTraits< elementType, Comm >::GEOFaceType GEOFaceType;
    typedef typename ALU3dImplTraits< elementType, Comm >::GEOEdgeType GEOEdgeType;
    typedef typename ALU3dImplTraits< elementType, Comm >::GEOVertexType GEOVertexType;

    // interface types
    typedef typename ALU3dImplTraits< elementType, Comm >::HFaceType HFaceType;
    typedef typename ALU3dImplTraits< elementType, Comm >::HEdgeType HEdgeType;
    typedef typename ALU3dImplTraits< elementType, Comm >::VertexType VertexType;

    typedef ElementTopologyMapping<elementType> ElementTopo;
    typedef FaceTopologyMapping<elementType> FaceTopo;

    enum { corners_      = (elementType == hexa) ? StaticPower<2,mydim>::power : mydim+1 };

    // type of specialized geometry implementation
    typedef typename MyALUGridGeometryImplementation<cdim> ::
    template GeometryImpl<0, mydim, elementType > GeometryImplType;

  public:
    typedef typename GridImp :: ctype ctype;

    //! type of local coordinates
    typedef FieldVector<ctype, mydim> LocalCoordinate;

    //! type of the global coordinates
    typedef FieldVector<ctype, cdim > GlobalCoordinate;

    //! type of jacobian inverse transposed
    typedef FieldMatrix<ctype,cdim,mydim> JacobianInverseTransposed;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

    // type of coordinate matrix for faces
    typedef FieldMatrix<ctype,
        EntityCount< elementType > :: numVerticesPerFace , 3> FaceCoordinatesType;

    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    ALU3dGridGeometry();

    //! copy constructor copying pointer and increasing reference count
    ALU3dGridGeometry( const ALU3dGridGeometry& );

    //! copy constructor copying pointer and increasing reference count
    ALU3dGridGeometry& operator = ( const ALU3dGridGeometry& );

    //! destructor decreasing reference count and freeing object
    ~ALU3dGridGeometry( );

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    GeometryType type () const;

    //! return the number of corners of this element. Corners are numbered 0..n-1
    int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    GlobalCoordinate corner (int i) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    GlobalCoordinate global (const LocalCoordinate& local) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    LocalCoordinate local (const GlobalCoordinate& global) const;

    //! A(l) , see grid.hh
    ctype integrationElement (const LocalCoordinate& local) const;

    //! can only be called for dim=dimworld! (Trivially true, since there is no
    //! other specialization...)
    const JacobianInverseTransposed &jacobianInverseTransposed (const LocalCoordinate& local) const;

    //! jacobian transposed
    const JacobianTransposed& jacobianTransposed (const LocalCoordinate& local) const;

    //! returns true if mapping is affine
    inline bool affine () const;

    //! returns volume of geometry
    ctype volume () const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry out of a given ALU3dGridElement
    bool buildGeom(const IMPLElementType & item);
    bool buildGeom(const HFaceType & item, int twist, int faceNum);
    bool buildGeom(const HEdgeType & item, int twist, int);
    bool buildGeom(const VertexType & item, int twist, int);

    // this method is used by the intersection iterator
    bool buildGeom(const FaceCoordinatesType& coords);

    // this method is used by the intersection iterator
    template <class coord_t>
    bool buildGeom(const coord_t& p0,
                   const coord_t& p1,
                   const coord_t& p2,
                   const coord_t& p3);

    // this method is used by the intersection iterator
    template <class coord_t>
    bool buildGeom(const coord_t& p0,
                   const coord_t& p1,
                   const coord_t& p2);

    //! build geometry of local coordinates relative to father
    template <class GeometryType>
    bool buildGeomInFather(const GeometryType &fatherGeom , const GeometryType & myGeom);

    //! print internal data
    //! no interface method
    void print (std::ostream& ss) const;

    //! invalidate geometry implementation to avoid errors
    void invalidate () ;

    //! invalidate geometry implementation to avoid errors
    bool valid () const ;

  protected:
    //! assign pointer
    void assign( const ALU3dGridGeometry& other );
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
    struct StoreGeometryReference< mydim, cdim, GridImp, ALU3dGridGeometry >
    {
      //! Whether to store by reference.
      static const bool v = false;
    };
  }

} // end namespace Dune

#include "geometry_imp.cc"

#endif
