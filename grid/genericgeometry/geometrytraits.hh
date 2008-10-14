// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH
#define DUNE_GENERICGEOMETRY_GEOMETRYTRAITS_HH

#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/matrix.hh>
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
    };



    // MappingTraits
    // -------------

    template< class CT, unsigned int dim, unsigned int dimW >
    struct MappingTraits
    {
      typedef CT CoordTraits;

      static const unsigned int dimension = dim;
      static const unsigned int dimWorld = dimW;

      typedef typename CoordTraits :: ctype FieldType;
      typedef typename CoordTraits :: template Vector< dimension > :: type LocalCoordType;
      typedef typename CoordTraits :: template Vector< dimWorld > :: type GlobalCoordType;

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
      PreCompute,
    };



    // DefaultGeometryTraits
    // ---------------------

    template< class ctype, int dimG, int dimW >
    struct DefaultGeometryTraits
    {
      typedef DuneCoordTraits< ctype > CoordTraits;

      static const int dimGrid = dimG;
      static const int dimWorld = dimW;

      //   hybrid   [ true if Codim 0 is hybrid ]
      static const bool hybrid = true;
      //   dunetype [ for Codim 0, needed for (hybrid=false) ]
      // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;

      // what basic geometry type shall the line be considered?
      static const GeometryType :: BasicType linetype = GeometryType :: simplex;

      template< class Topology >
      struct Mapping
      {
        typedef CoordPointerStorage< CoordTraits, Topology, dimWorld > CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage, false > type;
      };

      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
        static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
        static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
        static const EvaluationType evaluateNormal = ComputeOnDemand;
      };
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
     *  can be used (via derivation) to provide the necessary information. It
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
     *    //   dunetype [ for Codim 0, needed for (hybrid=false) ]
     *    // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;
     *
     *    // what basic geometry type shall the line be considered?
     *    static const GeometryType :: BasicType linetype = GeometryType :: simplex;
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordType >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *      static const EvaluationType evaluateNormal = ComputeOnDemand;
     *    };
     *
     *    template< class Traits >
     *    struct Caching
     *    : public ComputeAll< Traits >
     *    {};
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
     *  can be used (via derivation) to provide the necessary information. It
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
     *    //   dunetype [ for Codim 0, needed for (hybrid=false) ]
     *    // static const GeometryType :: BasicType dunetype = GeometryType :: simplex;
     *
     *    // what basic geometry type shall the line be considered?
     *    static const GeometryType :: BasicType linetype = GeometryType :: simplex;
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordType >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *      static const EvaluationType evaluateNormal = ComputeOnDemand;
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

#endif
