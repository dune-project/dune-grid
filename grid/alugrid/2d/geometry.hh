// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGEOMETRY_HH
#define DUNE_ALU2DGRIDGEOMETRY_HH

// System includes

// Dune includes
#include <dune/common/misc.hh>
#include <dune/grid/common/grid.hh>

// Local includes

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

    template <class GeomImp, int gdim>
    struct CopyCoordinates;

    template <class GeomImp>
    struct CopyCoordinates<GeomImp,2>
    {
      template <class VecType>
      static inline void copyData(const double p[cdim], VecType& vec)
      {
        assert( cdim == 2 );
        vec[0] = p[0];
        vec[1] = p[1];
      }

      // copy coordinates and return determinant
      template<class CoordMatrixType>
      static inline double copy(const HElementType& item, CoordMatrixType& coord, const int )
      {
        assert( cdim == 2 );
        // copy coordinates
        copyData(item.getVertex(0)->coord(),coord[0]);
        copyData(item.getVertex(1)->coord(),coord[1]);
        copyData(item.getVertex(2)->coord(),coord[2]);
        return 2.0*item.area();
      }
    };

    template <class GeomImp>
    struct CopyCoordinates<GeomImp,1>
    {
      template <class VecType>
      static inline void copyData(const double p[cdim], VecType& vec)
      {
        assert( cdim == 2 );
        vec[0] = p[0];
        vec[1] = p[1];
      }

      // copy coordinates and return determinant
      template<class CoordMatrixType>
      static inline double copy(const HElementType& item, CoordMatrixType& coord, const int face)
      {
        // copy coordinates
        for(int i=0; i<2; ++i)
        {
          const int vx = (i+face+1)%3;
          copyData(item.getVertex(vx)->coord(),coord[i]);
        }
        return item.sidelength((face)%3);
      }
    };

  public:
    //! for makeRefGeometry == true a Geometry with the coordinates of the
    //! reference element is made
    ALU2dGridGeometry();

    //! constructor building geometry in father
    ALU2dGridGeometry(const int child, const int orientation );

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

    //! can only be called for dim=dimworld!
    const FieldMatrix<alu2d_ctype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<alu2d_ctype, mydim>& local) const;

    //***********************************************************************
    //!  Methods that not belong to the Interface, but have to be public
    //***********************************************************************
    //! generate the geometry for out of given ALU2dGridElement
    bool builtGeom(const ALU2DSPACE Vertex & item, int );
    bool builtGeom(const HElementType & item, int face);

    //! build geometry for intersectionSelfLocal and
    //! intersectionNeighborLocal
    template <class GeometryType, class LocalGeomType >
    bool builtLocalGeom(const GeometryType & geo , const LocalGeomType & lg);

    //! build local geometry given local face number
    bool builtLocalGeom(const int faceNumber, const int twist);

    //! return non-const reference to coord vecs
    FieldVector<alu2d_ctype, cdim>& getCoordVec (int i);

    //! print internal data
    void print (std::ostream& ss) const;

    //! build geometry with local coords of child in reference element
    //void buildGeomInFather(const int child, const int orientation);
    inline bool buildGeomInFather(const Geometry &fatherGeom , const Geometry & myGeom);

    inline bool up2Date() const { return up2Date_; }
    inline void unsetUp2Date() const { up2Date_ = false; }
  private:
    // set ref coords
    void setReferenceCoordinates();

    //! build the transposed of the jacobian inverse and store the volume
    void buildJacobianInverseTransposed () const;

    //! calculates the element matrix for calculation of the jacobian inverse
    void calcElMatrix () const;

    //! calculates the volume of the element
    alu2d_ctype elDeterminant () const;

    //! temporary need vector
    mutable FieldVector<alu2d_ctype, mydim+1> tmpVec_;

    //! the vertex coordinates
    mutable FieldMatrix<alu2d_ctype, mydim+1, cdim> coord_;

    //! the reference element coordinates, last is the length of face
    FieldMatrix<alu2d_ctype, 3 , 3> refCoord_;

    //! storage for global coords
    mutable FieldVector<alu2d_ctype, cdim> globalCoord_;

    //! storage for local coords
    mutable FieldVector<alu2d_ctype, mydim> localCoord_;

    enum { matdim = (mydim > 0) ? mydim : 1 };
    mutable FieldMatrix<alu2d_ctype,cdim,matdim> Jinv_; //!< storage for inverse of jacobian
    mutable FieldMatrix<alu2d_ctype,matdim,matdim> Mtmp_;    //!< storage for inverse of jacobian

    mutable FieldMatrix<alu2d_ctype,cdim,mydim> elMat_; //!< storage for mapping matrix
    mutable FieldMatrix<alu2d_ctype,matdim,matdim> elMatT_elMat_; //!< storage for mapping matrix

    //! is true if elMat_ was calced
    mutable bool builtElMat_;
    //! is true if Jinv_ and volume_ is calced
    mutable bool builtinverse_;

#ifndef NDEBUG
    mutable bool calcedDet_; //! true if determinant was calculated
#endif
    mutable alu2d_ctype elDet_;                           //!< storage of integration_element

    //! is true if geom is up2date
    mutable bool up2Date_;

    // temporary mem for integrationElement with mydim < cdim
    mutable FieldVector<alu2d_ctype,cdim> tmpV_; //! temporary memory
    mutable FieldVector<alu2d_ctype,cdim> tmpU_; //! temporary memory
    mutable FieldVector<alu2d_ctype,cdim> tmpZ_;

    mutable FieldVector<alu2d_ctype,mydim> AT_x_;

  };

  template <class GeometryImp, int nChild>
  class ALU2DLocalGeometryStorage {

    // array with pointers to the geometries
    std::vector < GeometryImp * > geoms_;
    // count local geometry creation
    int count_;
  public:
    // create empty storage
    ALU2DLocalGeometryStorage () : geoms_ (nChild) , count_ (0) {
      for(size_t i=0 ; i<geoms_.size(); ++i) geoms_[i] = 0;
    }

    // desctructor deleteing geometries
    ~ALU2DLocalGeometryStorage () {
      for(size_t i=0 ; i<geoms_.size(); ++i)
        if(geoms_[i]) delete geoms_[i];
    }

    // check if geometry has been created
    bool geomCreated(int child) const { return geoms_[child] != 0; }

    // create local geometry
    template <class GridImp, class Geometry>
    void create (const GridImp & grid, const Geometry & father, const Geometry & son, const int child)
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
