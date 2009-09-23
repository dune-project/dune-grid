// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGEOMETRY_HH
#define DUNE_ALU2DGRIDGEOMETRY_HH

// Dune includes
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

// alugrid mappings
#include <dune/grid/alugrid/3d/mappings.hh>

namespace Dune
{

  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int mydim, int cdim, class GridImp>
  class ALU2dGridGeometry;
  template<int dim, int dimworld>
  class ALU2dGrid;

  template <int cdim>
  class MyALU2dGridGeometryImplementation
  {
  public:
    typedef FieldVector<alu2d_ctype, cdim> CoordinateVectorType;

    struct CoordVecCopy
    {
      // copy coordinate vector from field vector or double[cdim]
      template <class CoordPtrType>
      static inline void copy(const CoordPtrType& p,
                              CoordinateVectorType& c)
      {
        assert( cdim == 2 );
        c[0] = p[0];
        c[1] = p[1];
      }
    };

    //! general type of geometry implementation
    template <int dummy, int dim> class GeometryImpl;
  public:
    // geometry implementation for edges and vertices
    template <int dummy, int dim>
    class GeometryImpl : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { corners_ = dim+1 };
      //! the vertex coordinates
      typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

      // for edges use LinearMapping<cdim, 1> here that has all features
      // implemented

      CoordinateMatrixType coord_;

    public:
      GeometryImpl() : coord_() {}
      GeometryImpl(const GeometryImpl& other) : coord_(other.coord_) {}

      // return true since we have affine mapping
      bool affine() const { return true; }

      // return coordinate vector
      inline const CoordinateVectorType& operator [] (const int i) const
      {
        assert( i>=0 && i<corners_ );
        return coord_[i];
      }

      // update vertex
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0)
      {
        assert( corners_ == 1 );
        copy( p0, coord_[0] );
      }
    };

    // geom impl for edges
    template <int dummy>
    class GeometryImpl<dummy,1> : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { corners_ = 2 };

      //! the vertex coordinates
      typedef FieldMatrix<alu2d_ctype, corners_ , cdim>  CoordinateMatrixType;
      typedef LinearMapping<cdim, 1>   MappingType;

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

      // return true since we have affine mapping
      bool affine() const { return liMap_.affine(); }

      // update geometry in father coordinates
      template <class GeometryImp, class LocalGeoImp>
      inline void updateLocal(const GeometryImp& globalGeom,
                              const LocalGeoImp& localGeom)
      {
        // compute the local coordinates in father refelem
        for(int i=0; i < localGeom.corners() ; ++i)
        {
          // calculate coordinate
          coord_[i] = globalGeom.local( localGeom.corner( i ) );

          // to avoid rounding errors
          for(int j=0; j<cdim; ++j)
          {
            if ( coord_[i][j] < 1e-14) coord_[i][j] = 0.0;
          }
        }

        builtMapping_ = false ;
      }

      // update geometry coordinates
      template <class CoordPtrType>
      inline void update(const CoordPtrType& p0,
                         const CoordPtrType& p1)
      {
        copy(p0, coord_[0] );
        copy(p1, coord_[1] );
        builtMapping_ = false;
      }

      // return mapping (always up2date)
      inline MappingType& mapping()
      {
        if( builtMapping_ ) return liMap_;

        liMap_.buildMapping( coord_[0], coord_[1] );
        builtMapping_ = true ;
        return liMap_;
      }
    };

    // geom impl for simplices
    template <int dummy>
    class GeometryImpl<dummy,2> : public CoordVecCopy
    {
      using CoordVecCopy :: copy ;

      enum { corners_ = 3 };

      //! the vertex coordinates
      typedef FieldMatrix<alu2d_ctype, corners_ , cdim>  CoordinateMatrixType;
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

      // return true since we have affine mapping
      bool affine() const { return liMap_.affine(); }

      // update geometry in father coordinates
      template <class GeometryImp, class LocalGeoImp>
      inline void updateLocal(const GeometryImp& globalGeom,
                              const LocalGeoImp& localGeom)
      {
        // compute the local coordinates in father refelem
        for(int i=0; i < localGeom.corners() ; ++i)
        {
          // calculate coordinate
          coord_[i] = globalGeom.local( localGeom.corner( i ) );

          // to avoid rounding errors
          for(int j=0; j<cdim; ++j)
          {
            if ( coord_[i][j] < 1e-14) coord_[i][j] = 0.0;
          }
        }
        builtMapping_ = false ;
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
  }; // end of class ALUGridGeometryImplementation

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
  template <int mydim, int cdim, class GridImp>
  class ALU2dGridGeometry :
    public GeometryDefaultImplementation <mydim,cdim,GridImp,ALU2dGridGeometry>
  {

    //! type of our Geometry interface
    typedef typename GridImp::template Codim<0>::Geometry Geometry;
    //! type of our Geometry implementation
    typedef ALU2dGridGeometry<mydim,cdim,GridImp> GeometryImp;
    //! know dimension of barycentric coordinates
    enum { dimbary=mydim+1};

    typedef typename ALU2DSPACE Hmesh_basic::helement_t HElementType ;

    // type of specialized geometry implementation
    typedef typename MyALU2dGridGeometryImplementation<cdim> ::
    template GeometryImpl<0, mydim> GeometryImplType;

  public:
    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    ALU2dGridGeometry();

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    const GeometryType type () const;

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<alu2d_ctype, cdim>& operator[] (int i) const;

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector<alu2d_ctype, cdim> corner (int i) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    FieldVector<alu2d_ctype, cdim> global (const FieldVector<alu2d_ctype, mydim>& local) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    FieldVector<alu2d_ctype,  mydim> local (const FieldVector<alu2d_ctype, cdim>& global) const;

#if 0
    //! returns true if the point in local coordinates is inside reference element
    bool checkInside(const FieldVector<alu2d_ctype, mydim>& local) const;
#endif

    //! A(l) , see grid.hh
    alu2d_ctype integrationElement (const FieldVector<alu2d_ctype, mydim>& local) const;

    //! return volume of geometry
    alu2d_ctype volume () const;

    //! return true if geometry has affine mapping
    bool affine() const { return geoImpl_.affine(); }

    //! can only be called for dim=dimworld!
    const FieldMatrix<alu2d_ctype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<alu2d_ctype, mydim>& local) const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for out of given ALU2dGridElement
    // method for elements
    bool buildGeom(const HElementType & item);
    // method for edges
    bool buildGeom(const HElementType & item, const int face);
    // method for vertices
    bool buildGeom(const ALU2DSPACE Vertex & item, const int );

    //! build geometry for intersectionSelfLocal and
    //! intersectionNeighborLocal
    template <class GeometryType, class LocalGeomType >
    bool buildLocalGeom(const GeometryType & geo , const LocalGeomType & lg);

    //! build local geometry given local face number
    bool buildLocalGeom(const int faceNumber, const int twist);

    //! return non-const reference to coord vecs
    FieldVector<alu2d_ctype, cdim>& getCoordVec (int i);

    //! print internal data
    void print (std::ostream& ss) const;

    //! build geometry with local coords of child in reference element
    inline bool buildGeomInFather(const Geometry &fatherGeom , const Geometry & myGeom);

    // returns true if geometry is up-2-date
    inline bool up2Date() const { return up2Date_; }

    // set up2Date flag to false
    inline void unsetUp2Date() const { up2Date_ = false; }

  protected:
    // implementation of coord and mapping
    mutable GeometryImplType geoImpl_;
    // determinant
    mutable alu3d_ctype det_;

    //! is true if geom is up2date
    mutable bool up2Date_;
  };

  template <class GeometryImp, int nChild>
  class ALU2DLocalGeometryStorage {

    // array with pointers to the geometries
    std::vector < GeometryImp * > geoms_;
    // count local geometry creation
    int count_;
  public:
    // create empty storage
    ALU2DLocalGeometryStorage () : geoms_ (nChild) , count_ (0)
    {
      for(size_t i=0 ; i<geoms_.size(); ++i) geoms_[i] = 0;
    }

    // desctructor deleteing geometries
    ~ALU2DLocalGeometryStorage ()
    {
      for(size_t i=0 ; i<geoms_.size(); ++i)
        if(geoms_[i]) delete geoms_[i];
    }

    // check if geometry has been created
    bool geomCreated(int child) const { return geoms_[child] != 0; }

    // create local geometry
    template <class GridImp, class Geometry>
    void create (const GridImp & grid,
                 const Geometry & father,
                 const Geometry & son, const int child)
    {
      assert( !geomCreated(child) );
      assert( child >=0 && child < nChild );

      assert( count_ < nChild );
      ++count_;

      typedef typename GeometryImp :: ImplementationType ImplType;
      GeometryImp * g = new GeometryImp(ImplType());
      geoms_[child] = g;
      GeometryImp & geo = *g;
      grid.getRealImplementation(geo).buildGeomInFather( father, son );
    }
    // return reference to local geometry
    const GeometryImp & operator [] (int child) const
    {
      assert( geomCreated(child) );
      return *(geoms_[child]);
    }
  };

} // end namespace Dune

#include "geometry_imp.cc"

#endif
