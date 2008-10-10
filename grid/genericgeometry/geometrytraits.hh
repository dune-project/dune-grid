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

    template< class CT, unsigned int DimG, unsigned int DimW >
    struct MappingTraits
    {
      typedef CT CoordTraits;

      static const unsigned int dimG = DimG;
      static const unsigned int dimW = DimW;

      typedef typename CoordTraits :: ctype FieldType;
      typedef typename CoordTraits :: template Vector< dimG > :: type LocalCoordType;
      typedef typename CoordTraits :: template Vector< dimW > :: type GlobalCoordType;

      typedef typename CoordTraits :: template Matrix< dimW, dimG > :: type
      JacobianType;
      typedef typename CoordTraits :: template Matrix< dimG, dimW > :: type
      JacobianTransposedType;

      typedef GenericGeometry :: MatrixHelper< CoordTraits > MatrixHelper;

      template< unsigned int codim >
      struct Codim
      {
        typedef GenericGeometry :: MappingTraits< CoordTraits, dimG - codim, dimW >
        MappingTraits;
      };
    };



    // If not affine only volume is cached (based on intElCompute)
    // otherwise all quantities can be cached using:
    //   geoCompute:    assign if method called using barycenter
    //   geoPreCompute: assign in constructor using barycenter
    //   geoIsComputed: assign in constructor using barycenter using callback
    enum EvaluationType
    {
      //! compute on demand
      ComputeOnDemand,
      //! compute in constructor
      PreCompute,
      //! assign in constructor using callback
      IsComputed
    };

    template< class Traits >
    struct ComputeAll
    {
      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
      static const EvaluationType evaluateNormal = ComputeOnDemand;

      void jacobianT ( typename Traits :: JacobianTransposedType &jT ) const
      {}

      void integrationElement ( typename Traits :: FieldType &intEl ) const
      {}

      void jacobianInverseTransposed ( typename Traits :: JacobianType &jTInv ) const
      {}

      void normal ( int face, typename Traits :: GlobalCoordType &n ) const
      {}
    };

    template< class Traits >
    struct PreComputeAll
    {
      static const EvaluationType evaluateJacobianTransposed = PreCompute;
      static const EvaluationType evaluateJacobianInverseTransposed = PreCompute;
      static const EvaluationType evaluateIntegrationElement = PreCompute;
      static const EvaluationType evaluateNormal = PreCompute;

      void jacobianT ( typename Traits :: JacobianTransposedType &jT ) const
      {}

      void integrationElement ( typename Traits :: FieldType &intEl ) const
      {}

      void jacobianInverseTransposed ( typename Traits :: JacobianType &jTInv ) const
      {}

      void normal ( int face, typename Traits :: GlobalCoordType &n ) const
      {}
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
        typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
        typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordType >
        CornerStorage;
        typedef CornerMapping< Topology, Traits, CornerStorage > type;
      };

      template< class Traits >
      struct Caching
        : public ComputeAll< Traits >
      {};
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
     *    template< class Traits >
     *    struct Caching
     *    : public ComputeAll< Traits >
     *    {};
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
