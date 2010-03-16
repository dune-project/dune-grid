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
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;


  template< int mydim, int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl;

  // geometry implementation for vertices
  template< int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl< 0, cdim, eltype >
  {
    enum { corners_ = 1 };

    //! the vertex coordinates
    typedef FieldMatrix<alu3d_ctype, corners_ , cdim>  CoordinateMatrixType;

    typedef LinearMapping< cdim, 0 > MappingType;

    MappingType liMap_;

  public:
    // return true since we have affine mapping
    bool affine() const { return true; }

    const MappingType &mapping() const { return liMap_; }

    // update vertex
    template< class CoordPtrType >
    inline void update( const CoordPtrType &p0 )
    {
      liMap_.buildMapping( p0 );
    }
  };

  template< int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl< 1, cdim, eltype >
  {
    enum { corners_ = 2 };

    //! the vertex coordinates
    typedef FieldMatrix<alu2d_ctype, corners_ , cdim>  CoordinateMatrixType;
    typedef LinearMapping< cdim, 1 >   MappingType;

    MappingType liMap_;

  public:
    // return true since we have affine mapping
    bool affine() const { return liMap_.affine(); }

    // update geometry in father coordinates
    template <class GeometryImp, class LocalGeoImp>
    inline void updateLocal(const GeometryImp& globalGeom,
                            const LocalGeoImp& localGeom)
    {
      CoordinateMatrixType coord;
      // compute the local coordinates in father refelem
      for(int i=0; i < localGeom.corners() ; ++i)
      {
        // calculate coordinate
        coord[i] = globalGeom.local( localGeom.corner( i ) );

        // to avoid rounding errors
        for(int j=0; j<cdim; ++j)
        {
          if ( coord[i][j] < 1e-14) coord[i][j] = 0.0;
        }
      }
      liMap_.buildMapping( coord[ 0 ], coord[ 1 ] );
    }

    // update geometry coordinates
    template <class CoordPtrType>
    inline void update(const CoordPtrType& p0,
                       const CoordPtrType& p1)
    {
      liMap_.buildMapping( p0, p1 );
    }

    // return mapping (always up2date)
    inline MappingType &mapping()
    {
      return liMap_;
    }
  };

  // geom impl for elements
  template< int cdim, ALU2DSPACE ElementType eltype >
  class MyALU2dGridGeometryImpl< 2, cdim, eltype >
  {
    enum { corners_ = (eltype == ALU2DSPACE triangle ? 3 : 4) };

    //! the vertex coordinates
    typedef FieldMatrix<alu2d_ctype, corners_ , cdim>  CoordinateMatrixType;

    typedef LinearMapping< cdim, 2 > LinearMappingType;
    typedef BilinearMapping< cdim > BilinearMappingType;

    typedef LinearMappingType MappingType;

    MappingType liMap_;

  public:
    // return true since we have affine mapping
    bool affine() const { return liMap_.affine(); }

    // update geometry in father coordinates
    template <class GeometryImp, class LocalGeoImp>
    inline void updateLocal(const GeometryImp& globalGeom,
                            const LocalGeoImp& localGeom)
    {
      CoordinateMatrixType coord;
      // compute the local coordinates in father refelem
      for(int i=0; i < localGeom.corners() ; ++i)
      {
        // calculate coordinate
        coord[i] = globalGeom.local( localGeom.corner( i ) );

        // to avoid rounding errors
        for(int j=0; j<cdim; ++j)
        {
          if ( coord[i][j] < 1e-14) coord[i][j] = 0.0;
        }
      }
      liMap_.buildMapping( coord[ 0 ], coord[ 1 ], coord[ 2 ] );
    }

    // update geometry coordinates
    template <class CoordPtrType>
    inline void update(const CoordPtrType& p0,
                       const CoordPtrType& p1,
                       const CoordPtrType& p2)
    {
      liMap_.buildMapping( p0, p1, p2 );
    }

    // return mapping (always up2date)
    inline MappingType& mapping()
    {
      return liMap_;
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

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    const GeometryType type () const;

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const;

    //! access to coordinates of corners. Index is the number of the corner
    const GlobalCoordinate &operator[] ( int i ) const;

    //! access to coordinates of corners. Index is the number of the corner
    GlobalCoordinate corner ( int i ) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    GlobalCoordinate global ( const LocalCoordinate &local ) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    FieldVector<alu2d_ctype,  mydim> local (const FieldVector<alu2d_ctype, cdim>& global) const;

    //! A(l) , see grid.hh
    alu2d_ctype integrationElement (const FieldVector<alu2d_ctype, mydim>& local) const;

    //! return volume of geometry
    alu2d_ctype volume () const;

    //! return true if geometry has affine mapping
    bool affine() const { return geoImpl_.affine(); }

    //! jacobian inverse transposed
    const FieldMatrix<alu2d_ctype,cdim,mydim>& jacobianInverseTransposed (const FieldVector<alu2d_ctype, mydim>& local) const;

    //! jacobian transposed
    const FieldMatrix<alu2d_ctype,mydim,cdim>& jacobianTransposed (const FieldVector<alu2d_ctype, mydim>& local) const;

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
    bool buildLocalGeometry(const int faceNumber, const int twist);

    //! return non-const reference to coord vecs
    FieldVector<alu2d_ctype, cdim>& getCoordVec (int i);

    //! print internal data
    void print (std::ostream& ss) const;

    //! build geometry with local coords of child in reference element
    inline bool buildGeomInFather(const Geometry &fatherGeom ,
                                  const Geometry & myGeom,
                                  const bool hasBndProjection = false );

    // returns true if geometry is up-2-date
    inline bool up2Date() const { return up2Date_; }

    // set up2Date flag to false
    inline void unsetUp2Date() const { up2Date_ = false; }

  protected:
    // return reference coordinates of the alu triangle
    FieldMatrix<alu2d_ctype, 3, 3> calculateReferenceCoords() const;

    // implementation of coord and mapping
    mutable GeometryImplType geoImpl_;

    // determinant
    mutable alu2d_ctype det_;

    //! is true if geom is up2date
    mutable bool up2Date_;

#ifndef NDEBUG
    //! true if boundary projection is set
    mutable bool haveProjection_;
#endif
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
      grid.getRealImplementation(geo).
      buildGeomInFather( father, son, grid.hasBoundaryProjection() );
    }
    // return reference to local geometry
    const GeometryImp & operator [] (int child) const
    {
      assert( geomCreated(child) );
      return *(geoms_[child]);
    }
  };

  template <class LocalGeometry, class LocalGeometryImp>
  class ALU2DIntersectionGeometryStorage
  {
    // one geometry for each face and twist 0 and 1
    LocalGeometry* geoms_[ 3 ][ 2 ];

  protected:
    // create empty storage
    ALU2DIntersectionGeometryStorage ()
    {
      for( int i=0; i<3; ++i)
      {
        for( int j=0; j<2; ++j)
        {
          LocalGeometryImp geo;
          // build geometry
          geo.buildLocalGeometry( i, j );
          // create dune geoemtry
          geoms_[ i ][ j ] = new LocalGeometry( geo );
        }
      }
    }

  public:
    // destructor
    ~ALU2DIntersectionGeometryStorage()
    {
      for( size_t i=0; i<3; ++i)
        for( size_t j=0; j<2; ++j)
          delete geoms_[ i ][ j ];
    }

    typedef ALU2DIntersectionGeometryStorage<LocalGeometry, LocalGeometryImp> ThisType;

    // return reference to local geometry
    const LocalGeometry& localGeom(const int aluFace, const int twist) const
    {
      assert( 0 <= aluFace && aluFace < 3 );
      assert( twist == 0 || twist == 1 );
      assert( geoms_[ aluFace ][ twist ] );
      return *geoms_[ aluFace ][ twist ];
    }

    // return static instance
    static const ThisType& instance()
    {
      static const ThisType geomStorage;
      return geomStorage;
    }
  };

} // end namespace Dune

#include "geometry_imp.cc"

#endif
