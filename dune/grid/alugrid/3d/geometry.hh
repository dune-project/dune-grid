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
  template<int dim, int dimworld, ALU3dGridElementType elType>
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
      bool builtMapping_;

    public:
      GeometryImpl() : coord_() , liMap_() , builtMapping_(false) {}
      GeometryImpl(const GeometryImpl& other)
        : coord_(other.coord_),
          liMap_(other.liMap_),
          builtMapping_(other.builtMapping_)
      {}

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      inline MappingType& mapping()
      {
        if( builtMapping_ ) return liMap_;

        liMap_.buildMapping( coord_[0] );
        builtMapping_ = true ;
        return liMap_;
      }

      // update vertex
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0)
      {
        assert( corners_ == 1 );
        copy( p0, coord_[0] );
      }
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
      bool builtMapping_;

    public:
      GeometryImpl() : coord_() , liMap_() , builtMapping_(false) {}
      GeometryImpl(const GeometryImpl& other)
        : coord_(other.coord_),
          liMap_(other.liMap_),
          builtMapping_(other.builtMapping_)
      {}

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      inline MappingType& mapping()
      {
        if( builtMapping_ ) return liMap_;

        liMap_.buildMapping( coord_[0], coord_[1] );
        builtMapping_ = true ;
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
        builtMapping_ = false;
      }

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
      bool builtMapping_;

    public:
      // constructor
      GeometryImpl() : coord_(), liMap_(), builtMapping_(false) {}
      // copy constructor
      GeometryImpl(const GeometryImpl& other) :
        coord_(other.coord_),
        liMap_(other.liMap_),
        builtMapping_(other.builtMapping_)
      {}

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
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
        builtMapping_ = false;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        if( builtMapping_ ) return liMap_;

        liMap_.buildMapping( coord_[0], coord_[1], coord_[2] );
        builtMapping_ = true ;
        return liMap_;
      }
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
      bool builtMapping_;

    public:
      // constructor
      GeometryImpl() : coord_(), biMap_(), builtMapping_(false) {}
      // copy constructor
      GeometryImpl(const GeometryImpl& other) :
        coord_(other.coord_),
        biMap_(other.biMap_),
        builtMapping_(other.builtMapping_)
      {}

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
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
        builtMapping_ = false;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        if( builtMapping_ ) return biMap_;

        biMap_.buildMapping( coord_[0], coord_[1], coord_[2], coord_[3] );
        builtMapping_ = true ;
        return biMap_;
      }
    };

    // geometry impl for hexahedrons
    template <int dummy>
    class GeometryImpl<dummy,3, hexa>
    {
      enum { corners_ = 8 };

      //! the vertex coordinates
      typedef alu3d_ctype CoordPtrType[cdim];

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      typedef TrilinearMapping MappingType;

      // coordinate pointer vector
      FieldVector<const CoordPtrType*, corners_ > coordPtr_;
      MappingType triMap_;
      CoordinateMatrixType* fatherCoord_;
      bool builtMapping_;

    public:
      GeometryImpl() : coordPtr_((CoordPtrType*) 0), triMap_(),
                       fatherCoord_(0),
                       builtMapping_(false)
      {}


      GeometryImpl(const GeometryImpl& other) :
        coordPtr_(other.coordPtr_),
        triMap_(other.triMap_),
        fatherCoord_(0),
        builtMapping_(other.builtMapping_)
      {
        // if father coords are set, then reset coordPtr
        if( other.fatherCoord_ )
        {
          fatherCoord_ = new CoordinateMatrixType(*other.fatherCoord_);
          CoordinateMatrixType& coord = *fatherCoord_;
          for(int i=0; i<corners_; ++i)
          {
            coordPtr_[i] = reinterpret_cast<const CoordPtrType*> (&(coord[i][0]));
          }
        }
      }

      // desctructor
      ~GeometryImpl()
      {
        if( fatherCoord_ ) delete fatherCoord_;
      }

      // return coordinates
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( i>=0 && i<corners_ );
        return reinterpret_cast<const CoordinateVectorType&> (*(coordPtr_[i]));
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
        coordPtr_[0] = &p0;
        coordPtr_[1] = &p1;
        coordPtr_[2] = &p2;
        coordPtr_[3] = &p3;
        coordPtr_[4] = &p4;
        coordPtr_[5] = &p5;
        coordPtr_[6] = &p6;
        coordPtr_[7] = &p7;
        builtMapping_ = false;
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
          coordPtr_[i] = reinterpret_cast<const CoordPtrType*> (&(coord[i][0]));

          // to avoid rounding errors
          for(int j=0; j<cdim; ++j)
          {
            if ( coord[i][j] < 1e-16) coord[i][j] = 0.0;
          }
        }

        builtMapping_ = false ;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        if( builtMapping_ ) return triMap_;

        triMap_.buildMapping( (*this)[0], (*this)[1], (*this)[2], (*this)[3],
                              (*this)[4], (*this)[5], (*this)[6], (*this)[7] );

        builtMapping_ = true;
        return triMap_;
      }
    };


    // geometry impl for hexahedrons
    template <int dummy>
    class GeometryImpl<dummy,3, tetra>
    {
      enum { corners_ = 4 };

      //! the vertex coordinates
      typedef alu3d_ctype CoordPtrType[cdim];

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      typedef LinearMapping<cdim, cdim> MappingType;

      // coordinate pointer vector
      FieldVector<const CoordPtrType*, corners_ > coordPtr_;
      MappingType liMap_;
      CoordinateMatrixType* fatherCoord_;
      bool builtMapping_;

    public:
      // default constructor
      GeometryImpl() : coordPtr_((CoordPtrType*) 0), liMap_(),
                       fatherCoord_(0),
                       builtMapping_(false)
      {}

      // copy constructor
      GeometryImpl(const GeometryImpl& other) :
        coordPtr_(other.coordPtr_),
        liMap_(other.liMap_),
        fatherCoord_(0),
        builtMapping_(other.builtMapping_)
      {
        // if father coords are set, then reset coordPtr
        if( other.fatherCoord_ )
        {
          fatherCoord_ = new CoordinateMatrixType(*other.fatherCoord_);
          CoordinateMatrixType& coord = *fatherCoord_;
          for(int i=0; i<corners_; ++i)
          {
            coordPtr_[i] = reinterpret_cast<const CoordPtrType*> (&(coord[i][0]));
          }
        }
      }

      // destructor
      ~GeometryImpl()
      {
        if( fatherCoord_ ) delete fatherCoord_;
      }

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( i>=0 && i<corners_ );
        return reinterpret_cast<const CoordinateVectorType&> (*(coordPtr_[i]));
      }

      // update geometry coordinates
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1,
                         const CoordPtrType& p2,
                         const CoordPtrType& p3)
      {
        coordPtr_[0] = &p0;
        coordPtr_[1] = &p1;
        coordPtr_[2] = &p2;
        coordPtr_[3] = &p3;
        builtMapping_ = false;
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
          coordPtr_[i] = reinterpret_cast<const CoordPtrType*> (&(coord[i][0]));

          // to avoid rounding errors
          for(int j=0; j<cdim; ++j)
          {
            if ( coord[i][j] < 1e-16) coord[i][j] = 0.0;
          }
        }

        builtMapping_ = false ;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        if( builtMapping_ ) return liMap_;

        liMap_.buildMapping( (*this)[0], (*this)[1], (*this)[2], (*this)[3] );

        builtMapping_ = true;
        return liMap_;
      }
    };
  }; // end of class ALUGridGeometryImplementation

  template <int mydim, int cdim, class GridImp>
  class ALU3dGridGeometry :
    public GeometryDefaultImplementation<mydim, cdim, GridImp, ALU3dGridGeometry>
  {
    static const ALU3dGridElementType elementType = GridImp :: elementType ;

    friend class ALU3dGridIntersectionIterator<GridImp>;

    typedef typename ALU3dImplTraits<elementType>::IMPLElementType IMPLElementType;
    typedef typename ALU3dImplTraits<elementType>::PLLBndFaceType PLLBndFaceType;
    typedef typename ALU3dImplTraits<elementType>::GEOFaceType GEOFaceType;
    typedef typename ALU3dImplTraits<elementType>::GEOEdgeType GEOEdgeType;
    typedef typename ALU3dImplTraits<elementType>::GEOVertexType GEOVertexType;

    // interface types
    typedef typename ALU3dImplTraits<elementType>::HFaceType HFaceType;
    typedef typename ALU3dImplTraits<elementType>::HEdgeType HEdgeType;
    typedef typename ALU3dImplTraits<elementType>::VertexType VertexType;

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
    const GlobalCoordinate& operator[] (int i) const;

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

  private:
    // implementation of coord and mapping
    mutable GeometryImplType geoImpl_;
    // volume
    mutable ctype volume_;
  };

} // end namespace Dune

#include "geometry_imp.cc"

#endif
