// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRY_HH
#define DUNE_ALBERTA_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>

#include <dune/grid/genericgeometry/geometry.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/elementinfo.hh>

// set to 1 to use generic geometries in AlbertaGrid
#ifndef USE_GENERICGEOMETRY
#define USE_GENERICGEOMETRY 0
#endif

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



  // AlbertaGridCoordinateReader
  // ---------------------------

  template< int codim, class GridImp >
  struct AlbertaGridCoordinateReader
  {
    typedef typename remove_const< GridImp >::type Grid;

    static const int dimension = Grid::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;
    static const int coorddimension = Grid::dimensionworld;

    typedef Alberta::Real ctype;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef FieldVector< ctype, coorddimension > Coordinate;

  private:
    const Grid &grid_;
    const ElementInfo &elementInfo_;
    const int subEntity_;

  public:
    AlbertaGridCoordinateReader ( const GridImp &grid,
                                  const ElementInfo &elementInfo,
                                  int subEntity )
      : grid_( grid ),
        elementInfo_( elementInfo ),
        subEntity_( subEntity )
    {}

    void coordinate ( int i, Coordinate &x ) const
    {
      assert( !elementInfo_ == false );
      assert( (i >= 0) && (i <= mydimension) );

      const int k = mapVertices( subEntity_, i );
      const Alberta::GlobalVector &coord = grid_.getCoord( elementInfo_, k );
      for( int j = 0; j < coorddimension; ++j )
        x[ j ] = coord[ j ];
    }

    bool hasDeterminant () const
    {
      return ((codimension == 0) && elementInfo_.isLeaf());
    }

    ctype determinant () const
    {
      assert( hasDeterminant() );

      const Alberta::Element *el = elementInfo_.el();
      typedef typename Grid::LeafDataType::Data LeafData;
      LeafData *leafdata = (LeafData *)el->child[ 1 ];
      assert( leafdata != NULL );
      return leafdata->determinant;
    }

  private:
    static int mapVertices ( int subEntity, int i )
    {
      return Alberta::MapVertices< dimension, codimension >::apply( subEntity, i );
    }
  };



  // AlbertaGridCoordStorage
  // -----------------------

  template< class CoordTraits, class Topology, unsigned int dimW >
  class AlbertaGridCornerStorage
  {
    typedef AlbertaGridCornerStorage< CoordTraits, Topology, dimW > This;

  public:
    static const unsigned int size = Topology::numCorners;

    static const unsigned int dimWorld = dimW;

    typedef typename CoordTraits::template Vector< dimWorld >::type
    GlobalCoordinate;

    template< class SubTopology >
    struct SubStorage
    {
      typedef AlbertaGridCornerStorage< CoordTraits, SubTopology, dimWorld > type;
    };

  private:
    GlobalCoordinate coords_[ size ];

  public:
    template< class CoordReader >
    explicit AlbertaGridCornerStorage ( const CoordReader &coordReader )
    {
      for( unsigned int i = 0; i < size; ++i )
        coordReader.coordinate( i, coords_[ i ] );
    }

    template< class Mapping, unsigned int codim >
    explicit AlbertaGridCornerStorage ( const GenericGeometry::SubMappingCoords< Mapping, codim > &coords )
    {
      for( unsigned int i = 0; i < size; ++i )
        coords_[ i ] = coords[ i ];
    }

    const GlobalCoordinate &operator[] ( unsigned int i ) const
    {
      return coords_[ i ];
    }
  };



  // template< int dim, int dimworld, int cdim >
  template <class GridImp,int cdim>
  struct AlbertaGridGeometryTraits
  {
    // typedef AlbertaGrid< dim, dimworld > Grid;
    typedef typename remove_const<GridImp>::type Grid;

    typedef GenericGeometry::DuneCoordTraits< Alberta::Real > CoordTraits;

    static const int dimGrid = Grid::dimension;
    static const int dimWorld = cdim;

    static const bool hybrid = false;
    static const GeometryType::BasicType dunetype = GeometryType::simplex;

    static const GeometryType::BasicType linetype = GeometryType::simplex;

    template< class Topology >
    struct Mapping
    {
      typedef AlbertaGridCornerStorage< CoordTraits, Topology, dimWorld > CornerStorage;
      typedef GenericGeometry::CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage > type;
    };

    struct Caching
    {
      static const GenericGeometry::EvaluationType evaluateJacobianTransposed = GenericGeometry::ComputeOnDemand;
      static const GenericGeometry::EvaluationType evaluateJacobianInverseTransposed = GenericGeometry::ComputeOnDemand;
      static const GenericGeometry::EvaluationType evaluateIntegrationElement = GenericGeometry::ComputeOnDemand;
      static const GenericGeometry::EvaluationType evaluateNormal = GenericGeometry::ComputeOnDemand;
    };
  };



  // AlbertaGridGeometry
  // -------------------

#if USE_GENERICGEOMETRY
  template< int mydim, int cdim, class GridImp >
  class AlbertaGridGeometry
    : public GenericGeometry::BasicGeometry
      < mydim, AlbertaGridGeometryTraits< GridImp, cdim > >
  {
    typedef AlbertaGridGeometry< mydim, cdim, GridImp > This;
    typedef GenericGeometry::BasicGeometry
    < mydim, AlbertaGridGeometryTraits< GridImp, cdim > >
    Base;

  public:
    //! Default constructor
    AlbertaGridGeometry ()
      : Base ()
    {}

    AlbertaGridGeometry ( const This &other )
      : Base ( other )
    {}

    template< class CoordReader >
    AlbertaGridGeometry ( const CoordReader &coordReader )
      : Base( GeometryType( GeometryType::simplex, mydim ), coordReader )
    {}

    template< class CoordReader >
    void build ( const CoordReader &coordReader )
    {
      (*this) = AlbertaGridGeometry( coordReader );
    }
  };
#endif // #if USE_GENERICGEOMETRY

#if !USE_GENERICGEOMETRY
  /** \class AlbertaGridGeometry
   *  \brief geometry implementation for AlbertaGrid
   *
   *  Defines the geometry part of a mesh entity. Works for all dimensions,
   *  element types and dim of world. Provides reference element and mapping
   *  between local and global coordinates.
   *
   *  \tparam  mydim    dimension of the element (0 <= dim <= 3)
   *  \tparam  cdim     dimension of global coordinates
   *  \tparam  GridImp  grid implementation
   *                    (always const AlbertaGrid< dim, dimworld >)
   */
  template< int mydim, int cdim, class GridImp >
  class AlbertaGridGeometry
  {
    typedef AlbertaGridGeometry< mydim, cdim, GridImp > This;

    // remember type of the grid
    typedef GridImp Grid;

    // dimension of barycentric coordinates
    static const int dimbary = mydim + 1;

  public:
    //! type of coordinates
    typedef Alberta::Real ctype;

    static const int dimension = Grid :: dimension;
    static const int mydimension = mydim;
    static const int codimension = dimension - mydimension;
    static const int coorddimension = cdim;

    typedef FieldVector< ctype, mydimension > LocalVector;
    typedef FieldVector< ctype, coorddimension > GlobalVector;

    typedef FieldMatrix< ctype, mydimension, coorddimension >
    JacobianTransposed;
    typedef FieldMatrix< ctype, coorddimension, mydimension >
    JacobianInverseTransposed;

  private:
    static const int numCorners = mydimension + 1;

    typedef FieldMatrix< ctype, numCorners, coorddimension > CoordMatrix;

  public:
    //! Default constructor
    AlbertaGridGeometry ();

    AlbertaGridGeometry ( const This &other );

    template< class CoordReader >
    AlbertaGridGeometry ( const CoordReader &coordReader );

    //! return the element type identifier
    //! line , triangle or tetrahedron, depends on dim
    GeometryType type () const;

    /** \brief obtain the number of corners of this element */
    int corners () const;

    GlobalVector corner ( const int i ) const;

    //! access to coordinates of corners. Index is the number of the corner
    const GlobalVector &operator[] ( const int i ) const;

    //! maps a local coordinate within reference element to
    //! global coordinate in element
    GlobalVector global ( const LocalVector &local ) const;

    //! maps a global coordinate within the element to a
    //! local coordinate in its reference element
    LocalVector local ( const GlobalVector &global ) const;

#if 0
    //! returns true if the point in local coordinates is inside reference element
    bool checkInside( const LocalVector &local ) const;
#endif

    /*!
       Copy from sgrid.hh:

       Integration over a general element is done by integrating over the reference element
       and using the transformation from the reference element to the global element as follows:
       \f[\int\limits_{\Omega_e} f(x) dx = \int\limits_{\Omega_{ref}} f(g(l)) A(l) dl \f] where
       \f$g\f$ is the local to global mapping and \f$A(l)\f$ is the integration element.

       For a general map \f$g(l)\f$ involves partial derivatives of the map (surface element of
       the first kind if \f$d=2,w=3\f$, determinant of the Jacobian of the transformation for
       \f$d=w\f$, \f$\|dg/dl\|\f$ for \f$d=1\f$).

       For linear elements, the derivatives of the map with respect to local coordinates
       do not depend on the local coordinates and are the same over the whole element.

       For a structured mesh where all edges are parallel to the coordinate axes, the
       computation is the length, area or volume of the element is very simple to compute.

       Each grid module implements the integration element with optimal efficieny. This
       will directly translate in substantial savings in the computation of finite element
       stiffness matrices.
     */

    // A(l)
    ctype integrationElement ( const LocalVector &local ) const;

    // volume if geometry
    ctype volume () const;

    const JacobianTransposed &jacobianTransposed () const;

    const JacobianTransposed &
    jacobianTransposed ( const LocalVector &local ) const
    {
      return jacobianTransposed();
    }

    const JacobianInverseTransposed &jacobianInverseTransposed () const;

    const JacobianInverseTransposed &
    jacobianInverseTransposed ( const LocalVector &local ) const
    {
      return jacobianInverseTransposed();
    }

    //***********************************************************************
    //  Methods that not belong to the Interface, but have to be public
    //***********************************************************************

    void invalidate ();

    template< class CoordReader >
    void build ( const CoordReader &coordReader );

    //! print internal data
    //! no interface method
    void print ( std::ostream &out ) const;

  private:
    // calculates the volume of the element
    ctype elDeterminant () const;

    //! the vertex coordinates
    CoordMatrix coord_;

    // storage for the transposed of the jacobian
    mutable JacobianTransposed jT_;

    // storage for the transposed inverse of the jacboian
    mutable JacobianInverseTransposed jTInv_;

    // has jT_ been computed, yet?
    mutable bool builtJT_;
    // has jTInv_ been computed, yet?
    mutable bool builtJTInv_;

    mutable bool calcedDet_; //!< true if determinant was calculated
    mutable ctype elDet_; //!< storage of element determinant
  };
#endif // #if !USE_GENERICGEOMETRY



  // AlbertaGridLocalGeometryProvider
  // --------------------------------

  template< class Grid >
  class AlbertaGridLocalGeometryProvider
  {
    typedef AlbertaGridLocalGeometryProvider< Grid > This;

  public:
    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;

    template< int codim >
    struct Codim
    {
      typedef Geometry< dimension-codim, dimension, Grid, AlbertaGridGeometry >
      LocalGeometry;
    };

    typedef typename Codim< 0 >::LocalGeometry LocalElementGeometry;
    typedef typename Codim< 1 >::LocalGeometry LocalFaceGeometry;

    static const int numChildren = 2;
    static const int numFaces = dimension + 1;

    static const int minFaceTwist = Alberta::Twist< dimension, dimension-1 >::minTwist;
    static const int maxFaceTwist = Alberta::Twist< dimension, dimension-1 >::maxTwist;
    static const int numFaceTwists = maxFaceTwist - minFaceTwist + 1;

  private:
    struct GeoInFatherCoordReader;
    struct FaceCoordReader;

    AlbertaGridLocalGeometryProvider ()
    {
      buildGeometryInFather();
      buildFaceGeometry();
    }

    ~AlbertaGridLocalGeometryProvider ()
    {
      for( int child = 0; child < numChildren; ++child )
      {
        delete geometryInFather_[ child ][ 0 ];
        delete geometryInFather_[ child ][ 1 ];
      }

      for( int i = 0; i < numFaces; ++i )
      {
        for( int j = 0; j < numFaceTwists; ++j )
          delete faceGeometry_[ i ][ j ];
      }
    }

    void buildGeometryInFather();
    void buildFaceGeometry();

  public:
    const LocalElementGeometry &
    geometryInFather ( int child, const int orientation = 1 ) const
    {
      assert( (child >= 0) && (child < numChildren) );
      assert( (orientation == 1) || (orientation == -1) );
      return *geometryInFather_[ child ][ (orientation + 1) / 2 ];
    }

    const LocalFaceGeometry &
    faceGeometry ( int face, int twist = 0 ) const
    {
      assert( (face >= 0) && (face < numFaces) );
      assert( (twist >= minFaceTwist) && (twist <= maxFaceTwist) );
      return *faceGeometry_[ face ][ twist - minFaceTwist ];
    }

    static const This &instance ()
    {
      static This theInstance;
      return theInstance;
    }

  private:
    template< int codim >
    static int mapVertices ( int subEntity, int i )
    {
      return Alberta::MapVertices< dimension, codim >::apply( subEntity, i );
    }

    const LocalElementGeometry *geometryInFather_[ numChildren ][ 2 ];
    const LocalFaceGeometry *faceGeometry_[ numFaces ][ numFaceTwists ];
  };

}

#endif // #ifndef DUNE_ALBERTA_GEOMETRY_HH
