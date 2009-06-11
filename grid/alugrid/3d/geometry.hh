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



  //**********************************************************************
  //
  // --ALU3dGridGeometry
  // --Geometry
  /*!
     Defines the geometry part of a mesh entity. Works for all dimensions, element types and dimensions
     of world. Provides reference element and mapping between local and global coordinates.
     The element may have different implementations because the mapping can be
     done more efficient for structured meshes than for unstructured meshes.

     dim: An element is a polygonal in a hyperplane of dimension dim. 0 <= dim <= 3 is typically
     dim=0 is a point.

     dimworld: Each corner is a point with dimworld coordinates.
   */
  //! ALU3dGridGeometry
  //! Empty definition, needs to be specialized for element type
  template <int mydim, int cdim, class GridImp>
  class ALU3dGridGeometry :
    public GeometryDefaultImplementation <mydim,cdim,GridImp,ALU3dGridGeometry> {};

  //! Specialisation for tetrahedra
  template <int mydim, int cdim>
  class ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, tetra> > :
    public GeometryDefaultImplementation<mydim, cdim, const ALU3dGrid<3, 3, tetra>,
        ALU3dGridGeometry> {
    typedef const ALU3dGrid<3, 3, tetra> GridImp;

    typedef typename ALU3dImplTraits<tetra>::IMPLElementType IMPLElementType;
    typedef typename ALU3dImplTraits<tetra>::PLLBndFaceType PLLBndFaceType;
    typedef typename ALU3dImplTraits<tetra>::GEOFaceType GEOFaceType;
    typedef typename ALU3dImplTraits<tetra>::GEOEdgeType GEOEdgeType;
    typedef typename ALU3dImplTraits<tetra>::GEOVertexType GEOVertexType;

    typedef ElementTopologyMapping<tetra> ElementTopo;
    typedef FaceTopologyMapping<tetra> FaceTopo;
    //! know number of vertices
    enum { corners_ = mydim+1 };
  public:
    typedef FieldMatrix<alu3d_ctype, 3, 3> FaceCoordinatesType;
    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    inline ALU3dGridGeometry();

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    inline GeometryType type () const;

    //! return the number of corners of this element. Corners are numbered 0...n-1
    inline int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    inline const FieldVector<alu3d_ctype, cdim>& operator[] (int i) const;

    //! access to coordinates of corners. Index is the number of the corner
    inline FieldVector<alu3d_ctype, cdim> corner (int i) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    inline FieldVector<alu3d_ctype, cdim> global (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    inline FieldVector<alu3d_ctype,  mydim> local (const FieldVector<alu3d_ctype, cdim>& global) const;

#if 0
    //! returns true if the point in local coordinates is inside reference element
    inline bool checkInside(const FieldVector<alu3d_ctype, mydim>& local) const;
#endif

    //! A(l) , see grid.hh
    inline alu3d_ctype integrationElement (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! can only be called for dim=dimworld!
    inline const FieldMatrix<alu3d_ctype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! returns true if mapping is affine
    inline bool affine () const;

    //! returns volume of geometry
    inline alu3d_ctype volume () const;
    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for out of given ALU3dGridElement
    inline bool buildGeom(const IMPLElementType & item);
    inline bool buildGeom(const ALU3DSPACE HFaceType & item, int twist, int face );
    inline bool buildGeom(const ALU3DSPACE HEdgeType & item, int twist, int );
    inline bool buildGeom(const ALU3DSPACE VertexType & item, int twist, int);

    // this method is used by the intersection iterator
    inline bool buildGeom(const FaceCoordinatesType& coords);
    // this method is used by the intersection iterator
    template <class coord_t>
    bool buildGeom(const coord_t& p0,
                   const coord_t& p1,
                   const coord_t& p2);

    //! build geometry of local coordinates relative to father
    template <class GeometryImpType>
    inline bool buildGeomInFather(const GeometryImpType &fatherGeom , const GeometryImpType & myGeom);

    //! print internal data
    //! no interface method
    inline void print (std::ostream& ss) const;

  private:
    // calculates determinant and volume
    void calculateDeterminant () const ;

    //! calculates the vertex index in the reference element out of a face index
    //! and a local vertex index
    inline int faceIndex(int faceIdx, int vtxIdx) const;

    // generate transposed Jacobian Inverse and calculate integration_element
    inline void buildJacobianInverseTransposed() const;

    // calculates the element matrix for calculation of the jacobian inverse
    inline void calcElMatrix () const;

    // copies the values of point to the values of coord
    template <class coord_t>
    inline void copyCoordVec(const coord_t& point,
                             FieldVector<alu3d_ctype,cdim> & coord ) const;

    //! the vertex coordinates
    mutable FieldMatrix<alu3d_ctype, corners_, cdim> coord_;

    mutable FieldMatrix<alu3d_ctype,cdim,mydim> Jinv_; //!< storage for inverse of jacobian
    mutable alu3d_ctype detDF_;              //!< storage of integration_element
    mutable alu3d_ctype volume_;             //!< storage of volume
    mutable FieldMatrix<alu3d_ctype, cdim , mydim> A_;    //!< transformation matrix (transposed)
    mutable FieldMatrix<alu3d_ctype, mydim, mydim> AT_A_;    //!< transformation matrix (transposed)

    mutable FieldVector<alu3d_ctype, mydim> AT_x_;
    mutable FieldVector<alu3d_ctype, mydim> localCoord_;
    mutable FieldVector<alu3d_ctype, cdim>  globalCoord_;

    mutable FieldVector<alu3d_ctype,cdim> tmpV_; //! temporary memory
    mutable FieldVector<alu3d_ctype,cdim> tmpU_; //! temporary memory

    const GeometryType myGeomType_;

    //! is true if Jinv_, A and detDF_ is calced
    mutable bool builtinverse_;
    mutable bool builtA_;
    mutable bool builtDetDF_;
  };

  //! Specialisation for hexahedra
  template <int mydim, int cdim>
  class ALU3dGridGeometry<mydim, cdim, const ALU3dGrid<3, 3, hexa> > :
    public GeometryDefaultImplementation<mydim, cdim, const ALU3dGrid<3, 3, hexa>,
        ALU3dGridGeometry> {
    typedef const ALU3dGrid<3, 3, hexa> GridImp;
    friend class ALU3dGridIntersectionIterator<GridImp>;

    typedef typename ALU3dImplTraits<hexa>::IMPLElementType IMPLElementType;
    typedef typename ALU3dImplTraits<hexa>::PLLBndFaceType PLLBndFaceType;
    typedef typename ALU3dImplTraits<hexa>::GEOFaceType GEOFaceType;
    typedef typename ALU3dImplTraits<hexa>::GEOEdgeType GEOEdgeType;
    typedef typename ALU3dImplTraits<hexa>::GEOVertexType GEOVertexType;

    typedef ElementTopologyMapping<hexa> ElementTopo;
    typedef FaceTopologyMapping<hexa> FaceTopo;

    enum { corners_ = Power_m_p<2,mydim>::power };

    // geometry implementation for edges and vertices
    template <int dummy, int dim>
    class GeometryImpl
    {
      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;
      typedef FieldVector<alu3d_ctype, cdim> CoordinateVectorType;

      CoordinateMatrixType coord_;

    public:
      GeometryImpl() : coord_() {}
      GeometryImpl(const GeometryImpl& other) : coord_(other.coord_) {}

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      // update edge
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1)
      {
        assert( corners_ == 2 );
        copyCoordVec( p0, coord_[0] );
        copyCoordVec( p1, coord_[1] );
      }

      // update vertex
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0)
      {
        assert( corners_ == 1 );
        copyCoordVec( p0, coord_[0] );
      }

    protected:
      template <class CoordPtrType>
      inline void copyCoordVec(const CoordPtrType& p,
                               CoordinateVectorType& c)
      {
        assert( cdim == 3 );
        c[0] = p[0];
        c[1] = p[1];
        c[2] = p[2];
      }
    };


    // geom impl for faces
    template <int dummy>
    class GeometryImpl<dummy,2>
    {
      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;
      typedef FieldVector<alu3d_ctype, cdim> CoordinateVectorType;
      typedef BilinearSurfaceMapping MappingType;

      //! the coordinates (for dim = 3 only a pointer)
      CoordinateMatrixType coord_;

      MappingType biMap_;
      bool builtMapping_;

    public:
      GeometryImpl() : coord_(), biMap_(), builtMapping_(false) {}
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
        copyCoordVec(p0, coord_[0] );
        copyCoordVec(p1, coord_[1] );
        copyCoordVec(p2, coord_[2] );
        copyCoordVec(p3, coord_[3] );
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

    protected:
      template <class CoordPtrType>
      inline void copyCoordVec(const CoordPtrType& p,
                               CoordinateVectorType& c)
      {
        assert( cdim == 3 );
        c[0] = p[0];
        c[1] = p[1];
        c[2] = p[2];
      }

    };

    // geometry impl for elements
    template <int dummy>
    class GeometryImpl<dummy,3>
    {
      //! the vertex coordinates
      typedef double CoordPtrType[cdim];

      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;
      typedef FieldVector<alu3d_ctype, cdim> CoordinateVectorType;

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

      ~GeometryImpl()
      {
        if( fatherCoord_ ) delete fatherCoord_;
      }

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

    // type of specialized geometry implementation
    typedef GeometryImpl<0, mydim> GeometryImplType;
  public:
    typedef FieldMatrix<alu3d_ctype, 4, 3> FaceCoordinatesType;

    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    ALU3dGridGeometry();

    //! Destructor
    ~ALU3dGridGeometry();

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    GeometryType type () const;

    //! return the number of corners of this element. Corners are numbered 0..n-1
    int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<alu3d_ctype, cdim>& operator[] (int i) const;

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector<alu3d_ctype, cdim> corner (int i) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    FieldVector<alu3d_ctype, cdim> global (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    FieldVector<alu3d_ctype,  mydim> local (const FieldVector<alu3d_ctype, cdim>& global) const;

#if 0
    //! returns true if the point in local coordinates is inside reference
    //! element
    bool checkInside(const FieldVector<alu3d_ctype, mydim>& local) const;
#endif

    //! A(l) , see grid.hh
    alu3d_ctype integrationElement (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! can only be called for dim=dimworld! (Trivially true, since there is no
    //! other specialization...)
    const FieldMatrix<alu3d_ctype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<alu3d_ctype, mydim>& local) const;

    //! returns true if mapping is affine
    inline bool affine () const;

    //! returns volume of geometry
    alu3d_ctype volume () const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry out of a given ALU3dGridElement
    bool buildGeom(const IMPLElementType & item);
    bool buildGeom(const ALU3DSPACE HFaceType & item, int twist, int faceNum);
    bool buildGeom(const ALU3DSPACE HEdgeType & item, int twist, int);
    bool buildGeom(const ALU3DSPACE VertexType & item, int twist, int);

    // this method is used by the intersection iterator
    bool buildGeom(const FaceCoordinatesType& coords);
    // this method is used by the intersection iterator
    template <class coord_t>
    bool buildGeom(const coord_t& p0,
                   const coord_t& p1,
                   const coord_t& p2,
                   const coord_t& p3);

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
    mutable alu3d_ctype volume_;
  };

} // end namespace Dune

#include "geometry_imp.cc"

#endif
