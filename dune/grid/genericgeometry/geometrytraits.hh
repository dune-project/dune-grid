// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH
#define DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH

#include <dune/common/geometrytype.hh>
#include <dune/common/polyallocator.hh>

#include <dune/grid/genericgeometry/matrixhelper.hh>
#include <dune/grid/genericgeometry/cornermapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // DuneCoordTraits
    // ---------------

    template< class ct >
    struct DuneCoordTraits
    {
      typedef ct ctype;

      template< int dim >
      struct Vector
      {
        typedef FieldVector< ctype, dim > type;
      };

      template< int rows, int cols >
      struct Matrix
      {
        typedef FieldMatrix< ctype, rows, cols > type;
      };

      // This limit is, e.g., used in the termination criterion of the Newton
      // scheme within the generic implementation of the method local
      static const ctype epsilon ()
      {
        return 1e-6;
      }
    };



    // MappingTraits
    // -------------
    /** \class MappingTraits
     *  \ingroup GenericGeometry
     *  \brief Default mapping traits using Dune::FieldVector and
     *  Dune::FieldMatrix
     */
    template< class CT, unsigned int dim, unsigned int dimW >
    struct MappingTraits
    {
      typedef CT CoordTraits;

      static const unsigned int dimension = dim;
      static const unsigned int dimWorld = dimW;

      typedef typename CoordTraits :: ctype FieldType;
      typedef typename CoordTraits :: template Vector< dimension > :: type LocalCoordinate;
      typedef typename CoordTraits :: template Vector< dimWorld > :: type GlobalCoordinate;

      typedef typename CoordTraits :: template Matrix< dimWorld, dimension > :: type
      JacobianType;
      typedef typename CoordTraits :: template Matrix< dimension, dimWorld > :: type
      JacobianTransposedType;

      typedef GenericGeometry :: MatrixHelper< CoordTraits > MatrixHelper;

      template< unsigned int codim >
      struct Codim
      {
        typedef GenericGeometry :: MappingTraits< CoordTraits, dimension - codim, dimWorld >
        MappingTraits;
      };
    };



    // If not affine only volume is cached (based on intElCompute)
    // otherwise all quantities can be cached using:
    //   ComputeOnDemand:    assign if method called using barycenter
    //   PreCompute:         assign in constructor using barycenter
    enum EvaluationType
    {
      //! compute on demand
      ComputeOnDemand,
      //! compute in constructor
      PreCompute
    };



    // DefaultGeometryTraits
    // ---------------------

    /** \class   DefaultGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   default settings for BasicGeometry
     *
     *  The class BasicGeometry requires a template argument <em>Traits</em>.
     *  These traits specify which reference mapping shall be used by the
     *  geometry and tweaks some performance settings.
     *
     *  This default implementation serves two purposed. Firstly, it documents
     *  the exprected parameters. Secondly, the user of BasicGeometry can
     *  derive his traits class from DefaultGeometryTraits. Then, only the
     *  non-default settings have to be specified. Moreover, deriving from
     *  DefaultGeometryTraits makes the user code more robust to changes in
     *  the generic geometries.
     *
     *  \note DefaultGeometryTraits can directly be used for the
     *        <em>Traits</em> argument of BasicGeometry.
     */
    template< class ctype, int dimG, int dimW, bool alwaysAffine = false >
    struct DefaultGeometryTraits
    {
      //! types needed in matrix-vector operations
      typedef DuneCoordTraits< ctype > CoordTraits;

      //! dimension of the grid
      static const int dimGrid = dimG;
      //! dimension of the world
      static const int dimWorld = dimW;

      /** \brief may the grid contain elements of different type?
       *
       *  If the elements (entities of codimension 0) may differ in topology
       *  type, the grid is called hybrid (and this parameter must be set
       *  to true). In this case, all methods of the geometry implementation
       *  are virtual (but no other branching for topology type is used).
       *
       *  If the grid is non-hybrid, <em>hybrid</em> can be set to false.
       *  In this case, virtual methods are not necessary and, hence, the
       *  geometries are a little faster.
       *
       *  If <em>hybrid</em> is set to false, an additional paramter
       *  <em>topologyId</em> is required.
       *  It specifies the topological type of all elements in the grid.
       *  Here's an example:
       *  \code
       *  static const unsigned int topologyId = SimplexTopology< dimGrid >::type::id;
       *  \endcode
       */
      static const bool hybrid = true;
      //   topologyId [ for Codim 0, needed for (hybrid=false) ]
      // static const unsigned int topologyId = SimlexTopology< dimGrid >::type::id;

      /** \brief specifies the reference mapping to be used
       *
       *  \tparam  Topology  type of topology for which the mapping
       *                     implementation is specified
       *
       *  This sturcture contains a single tydedef <em>type</em> specifying
       *  the implementation of the reference mapping. Basically, it looks like
       *  \code
       *  typedef CornerMapping< ... > type;
       *  \endcode
       */
      template< class Topology >
      struct Mapping
      {
        typedef CoordStorage< CoordTraits, Topology, dimWorld > CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage, alwaysAffine > type;
      };

      /** \brief specifies how constant values are to be cached
       *
       *  This structure contains 4 parameters of type
       *  GenericGeometry::EvaluationType:
       *  - evaluateJacobianTransposed
       *  - evaluateJacobianInverseTransposed
       *  - evaluateIntegrationElement
       *  .
       *  These parameters control how eager these evaluations shall be
       *  performed in the case of an affine mapping.
       */
      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
        static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
        static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
      };

      typedef PolyAllocator Allocator;
    };



    /** \struct  GlobalGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by GenericGeometry::Geometry
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via subclassing) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimG, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    typedef DuneCoordTraits< ctype > CoordTraits;
     *
     *    static const int dimGrid = dimG;
     *    static const int dimWorld = dimW;
     *
     *    //   hybrid   [ true if Codim 0 is hybrid ]
     *    static const bool hybrid = true;
     *    //   topologyId [ for Codim 0, needed for (hybrid=false) ]
     *    // static const unsigned int topologyId = SimlexTopology< dimGrid >::type::id;
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordinate >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *    };
     *  };
     *  \endcode
     *
     *  This implementation specifies the information used by
     *  GenericGeometry::Geometry.
     *
     *  \tparam  Grid  type of the grid, this traits class applies to
     */
    template< class Grid >
    struct GlobalGeometryTraits;

    template< class Grid >
    struct GlobalGeometryTraits< const Grid >
      : public GlobalGeometryTraits< Grid >
    {};



    /** \struct  LocalGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by GenericGeometry::LocalGeometry
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via subclassing) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimG, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    typedef DuneCoordTraits< ctype > CoordTraits;
     *
     *    static const int dimGrid = dimG;
     *    static const int dimWorld = dimW;
     *
     *    //   hybrid   [ true if Codim 0 is hybrid ]
     *    static const bool hybrid = true;
     *    //   topologyId [ for Codim 0, needed for (hybrid=false) ]
     *    // static const unsigned int topologyId = SimlexTopology< dimGrid >::type::id;
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordinate >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *    };
     *  };
     *  \endcode
     *
     *  This implementation specifies the information used by
     *  GenericGeometry::LocalGeometry.
     *
     *  \tparam  Grid  type of the grid, this traits class applies to
     */
    template< class Grid >
    struct LocalGeometryTraits;

    template< class Grid >
    struct LocalGeometryTraits< const Grid >
      : public LocalGeometryTraits< Grid >
    {};
  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH
