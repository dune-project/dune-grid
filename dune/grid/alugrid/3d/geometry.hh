// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDGEOMETRY_HH
#define DUNE_ALU3DGRIDGEOMETRY_HH

// System includes

// Dune includes
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

// Local includes
#include "alu3dinclude.hh"
#include "topology.hh"
#include "mappings.hh"

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

    static const signed char invalid      = -1; // means geometry is not meaningfull
    static const signed char updated      =  0; // means the point values have been set
    static const signed char buildmapping =  1; // means updated and mapping was build

    struct CoordVecCopy
    {

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
    };

    //! general type of geometry implementation
    template <int dummy, int dim,
        ALU3dGridElementType eltype> class GeometryImpl;
  public:
    // geometry implementation for edges and vertices
    template <int dummy, int dim, ALU3dGridElementType eltype>
    class GeometryImpl : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { corners_ = dim+1 };
      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      typedef LinearMapping<cdim, dim> MappingType;

      CoordinateMatrixType coord_;
      MappingType liMap_;
      signed char status_;

    public:
      using CoordVecCopy :: update ;

      GeometryImpl() : coord_() , liMap_() , status_( invalid ) {}
      GeometryImpl(const GeometryImpl& other)
        : coord_(other.coord_),
          liMap_(other.liMap_),
          status_(other.status_)
      {}

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
        if( status_ == buildmapping ) return liMap_;

        liMap_.buildMapping( coord_[0] );
        status_ = buildmapping ;
        return liMap_;
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

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
    };

    // geometry implementation for edges and vertices
    template <int dummy, ALU3dGridElementType eltype>
    class GeometryImpl<dummy,1,eltype> : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { dim = 1 };
      enum { corners_ = dim+1 };
      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      typedef LinearMapping<cdim, dim> MappingType;

      // for edges use LinearMapping<cdim, 1> here that has all features
      // implemented

      CoordinateMatrixType coord_;
      MappingType liMap_;
      signed char status_;

    public:
      using CoordVecCopy :: update ;

      GeometryImpl() : coord_() , liMap_() , status_( invalid ) {}
      GeometryImpl(const GeometryImpl& other)
        : coord_(other.coord_),
          liMap_(other.liMap_),
          status_(other.status_)
      {}

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
        if( status_ == buildmapping ) return liMap_;

        liMap_.buildMapping( coord_[0], coord_[1] );
        status_ = buildmapping ;
        return liMap_;
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

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
    };

    // geom impl for simplex faces (triangles)
    template <int dummy>
    class GeometryImpl<dummy,2, tetra> : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { corners_ = 3 };

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;
      typedef LinearMapping<cdim, 2>   MappingType;

      //! the coordinates (for dim = 3 only a pointer)
      CoordinateMatrixType coord_;

      MappingType liMap_;
      signed char status_;

    public:
      using CoordVecCopy :: update ;

      // constructor
      GeometryImpl() : coord_(), liMap_(), status_( invalid ) {}
      // copy constructor
      GeometryImpl(const GeometryImpl& other) :
        coord_(other.coord_),
        liMap_(other.liMap_),
        status_(other.status_)
      {}

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
        if( status_ == buildmapping ) return liMap_;

        liMap_.buildMapping( coord_[0], coord_[1], coord_[2] );
        status_ = buildmapping ;
        return liMap_;
      }

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
    };

    ///////////////////////////////////////////////////////////////
    //
    //  hexa specializations
    //
    ///////////////////////////////////////////////////////////////

    // geom impl for cube faces (quadrilaterals)
    template <int dummy>
    class GeometryImpl<dummy,2, hexa> : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { corners_ = 4 };

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;
      typedef BilinearSurfaceMapping MappingType;

      //! the coordinates (for dim = 3 only a pointer)
      CoordinateMatrixType coord_;

      MappingType biMap_;
      signed char status_;

    public:
      using CoordVecCopy :: update ;

      // constructor
      GeometryImpl() : coord_(), biMap_(), status_( invalid ) {}
      // copy constructor
      GeometryImpl(const GeometryImpl& other) :
        coord_(other.coord_),
        biMap_(other.biMap_),
        status_(other.status_)
      {}

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
        if( status_ == buildmapping ) return biMap_;

        biMap_.buildMapping( coord_[0], coord_[1], coord_[2], coord_[3] );
        status_ = buildmapping ;
        return biMap_;
      }

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
    };

    // geometry impl for hexahedrons
    template <int dummy>
    class GeometryImpl<dummy,3, hexa> : public CoordVecCopy
    {
      enum { corners_ = 8 };

      //! the vertex coordinates
      typedef alu3d_ctype CoordPtrType[cdim];

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      typedef TrilinearMapping MappingType;

      // coordinate pointer vector
      const alu3d_ctype* coordPtr_[ corners_ ];
      MappingType triMap_;
      CoordinateMatrixType* fatherCoord_;
      signed char status_;

    public:
      using CoordVecCopy :: update ;
      using CoordVecCopy :: copy ;

      //! constructor creating geo impl
      GeometryImpl() : triMap_(),
                       fatherCoord_(0),
                       status_( invalid )
      {
        // set initialize coord pointers
        for( int i=0; i<corners_; ++i )
          coordPtr_[ i ] = 0;
      }


      //! conpy constructor
      GeometryImpl(const GeometryImpl& other) :
        triMap_(other.triMap_),
        fatherCoord_(0),
        status_(other.status_)
      {
        // copy coord pointers
        for( int i=0; i<corners_; ++i )
          coordPtr_[ i ] = other.coordPtr_[ i ];

        // if father coords are set, then reset coordPtr
        if( other.fatherCoord_ )
        {
          fatherCoord_ = new CoordinateMatrixType(*other.fatherCoord_);
          CoordinateMatrixType& coord = *fatherCoord_;
          for(int i=0; i<corners_; ++i)
          {
            coordPtr_[i] = (&(coord[i][0]));
          }
        }
      }

      // desctructor
      ~GeometryImpl()
      {
        if( fatherCoord_ ) delete fatherCoord_;
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
        if( fatherCoord_ == 0 )
        {
          fatherCoord_ = new CoordinateMatrixType();
        }

        CoordinateMatrixType& coord = *fatherCoord_;
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
        if( status_ == buildmapping ) return triMap_;

        triMap_.buildMapping( point( 0 ), point( 1 ), point( 2 ), point( 3 ),
                              point( 4 ), point( 5 ), point( 6 ), point( 7 ) );

        status_ = buildmapping;
        return triMap_;
      }

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
    };


    // geometry impl for hexahedrons
    template <int dummy>
    class GeometryImpl<dummy,3, tetra> : public CoordVecCopy
    {
      enum { corners_ = 4 };

      //! the vertex coordinates
      typedef alu3d_ctype CoordPtrType[cdim];

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      typedef LinearMapping<cdim, cdim> MappingType;

      // coordinate pointer vector
      const alu3d_ctype* coordPtr_[ corners_ ];
      MappingType liMap_;
      CoordinateMatrixType* fatherCoord_;
      signed char status_;

    public:
      using CoordVecCopy :: update ;
      using CoordVecCopy :: copy ;

      // default constructor
      GeometryImpl() : liMap_(),
                       fatherCoord_(0),
                       status_( invalid )
      {
        // set initialize coord pointers
        for( int i=0; i<corners_; ++i )
          coordPtr_[ i ] = 0;
      }

      // copy constructor
      GeometryImpl(const GeometryImpl& other) :
        liMap_(other.liMap_),
        fatherCoord_(0),
        status_(other.status_)
      {
        // copy coord pointers
        for( int i=0; i<corners_; ++i )
          coordPtr_[ i ] = other.coordPtr_[ i ];

        // if father coords are set, then reset coordPtr
        if( other.fatherCoord_ )
        {
          fatherCoord_ = new CoordinateMatrixType(*other.fatherCoord_);
          CoordinateMatrixType& coord = *fatherCoord_;
          for(int i=0; i<corners_; ++i)
          {
            coordPtr_[i] = (&(coord[i][0]));
          }
        }
      }

      // destructor
      ~GeometryImpl()
      {
        if( fatherCoord_ ) delete fatherCoord_;
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
        if( fatherCoord_ == 0 )
        {
          fatherCoord_ = new CoordinateMatrixType();
        }

        CoordinateMatrixType& coord = *fatherCoord_;
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
        if( status_ == buildmapping ) return liMap_;

        liMap_.buildMapping( point( 0 ), point( 1 ), point( 2 ), point( 3 ) );

        status_ = buildmapping;
        return liMap_;
      }

      // set status to invalid
      void invalidate () { status_ = invalid ; }

      // return true if geometry is valid
      bool valid () const { return status_ != invalid ; }
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

    enum { corners_      = (elementType == hexa) ? Power_m_p<2,mydim>::power : mydim+1 };

    // type of specialized geometry implementation
    typedef typename MyALUGridGeometryImplementation<cdim> ::
    template GeometryImpl<0, mydim, elementType > GeometryImplType;

  public:
    typedef typename GridImp :: ctype ctype;

    //! type of local coordinates
    typedef FieldVector<ctype, mydim> LocalCoordinate;

    //! type of the global coordinates
    typedef FieldVector<ctype, cdim > GlobalCoordinate;

    //! type of jacobian (also of jacobian inverse transposed)
    typedef FieldMatrix<ctype,cdim,mydim> Jacobian;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

    // type of coordinate matrix for faces
    typedef FieldMatrix<ctype,
        EntityCount< elementType > :: numVerticesPerFace , 3> FaceCoordinatesType;

    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    ALU3dGridGeometry();

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
    const Jacobian& jacobianInverseTransposed (const LocalCoordinate& local) const;

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

  private:
    // implementation of coord and mapping
    mutable GeometryImplType geoImpl_;
    // volume
    mutable ctype volume_;
  };

} // end namespace Dune

#include "geometry_imp.cc"

#endif
