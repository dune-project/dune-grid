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
      return false;
    }

    ctype determinant () const
    {
      assert( hasDeterminant() );
      return ctype( 0 );
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
    AlbertaGridGeometry ();

    AlbertaGridGeometry ( const This &other );

    template< class CoordReader >
    AlbertaGridGeometry ( const CoordReader &coordReader );

    /** \brief obtain the type of reference element */
    GeometryType type () const;

    /** \brief number of corner the geometry */
    int corners () const;

    /** \brief obtain the i-th corner of this geometry */
    GlobalVector corner ( const int i ) const;

    /** \brief deprecated way of obtaining the i-th corner */
    const GlobalVector &operator[] ( const int i ) const;

    /** \brief map a point from the refence element to the geometry */
    GlobalVector global ( const LocalVector &local ) const;

    /** \brief map a point from the geometry to the reference element */
    LocalVector local ( const GlobalVector &global ) const;

    /** \brief integration element of the geometry mapping
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    ctype integrationElement () const;

    /** \brief integration element of the geometry mapping */
    ctype integrationElement ( const LocalVector &local ) const
    {
      return integrationElement();
    }

    /** \brief volume of geometry */
    ctype volume () const;

    /** \brief transposed of the geometry mapping's Jacobian
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    const JacobianTransposed &jacobianTransposed () const;

    /** \brief transposed of the geometry mapping's Jacobian */
    const JacobianTransposed &
    jacobianTransposed ( const LocalVector &local ) const
    {
      return jacobianTransposed();
    }

    /** \brief transposed inverse of the geometry mapping's Jacobian
     *
     *  \note This method is not part of the geometry interface. It is used
     *        internally only.
     */
    const JacobianInverseTransposed &jacobianInverseTransposed () const;

    /** \brief transposed inverse of the geometry mapping's Jacobian */
    const JacobianInverseTransposed &
    jacobianInverseTransposed ( const LocalVector &local ) const
    {
      return jacobianInverseTransposed();
    }

    //***********************************************************************
    //  Methods that not belong to the Interface, but have to be public
    //***********************************************************************

    /** \brief invalidate the geometry */
    void invalidate ();

    /** \brief build the geometry from a coordinate reader */
    template< class CoordReader >
    void build ( const CoordReader &coordReader );

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
