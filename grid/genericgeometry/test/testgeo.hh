// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef TESTGEO_HH
#define TESTGEO_HH

#include <dune/grid/alugrid/3d/topology.hh>
#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune
{

  template <int d1,int d2> class YaspGrid;
  template <int d1,int d2> class AlbertaGrid;
  template< int dim, int dimworld > class ALUSimplexGrid;
  template< int dim, int dimworld, ALU3dGridElementType elType > class ALU3dGrid;
  template <int dim> class UGGrid;



  // Topology
  // --------

  template< class Grid >
  struct Topology;

  template< int d >
  struct Topology< YaspGrid< d, d > >
  {
    static const GeometryType :: BasicType basicType = GeometryType :: cube;
    static const int dimension = d;
    static const bool correctJacobian = false;

    typedef GenericGeometry :: Convert< basicType, dimension > Convert;
    typedef typename Convert :: type Type;
  };

  template< int d >
  struct Topology< UGGrid< d > >
  {
    static const GeometryType :: BasicType basicType = GeometryType :: cube;
    static const int dimension = d;
    static const bool correctJacobian = false;

    typedef GenericGeometry :: Convert< basicType, dimension > Convert;
    typedef typename Convert :: type Type;
  };

  template< int d >
  struct Topology< AlbertaGrid< d, d > >
  {
    static const GeometryType :: BasicType basicType = GeometryType :: simplex;
    static const int dimension = d;
    static const bool correctJacobian = false;

    typedef GenericGeometry :: Convert< basicType, dimension > Convert;
    typedef typename Convert :: type Type;
  };

  template<>
  struct Topology< ALU3dGrid< 3, 3, tetra > >
  {
    static const GeometryType :: BasicType basicType = GeometryType :: simplex;
    static const int dimension = 3;
    static const bool correctJacobian = true;

    typedef GenericGeometry :: Convert< basicType, dimension > Convert;
    typedef Convert :: type Type;
  };

  template<>
  struct Topology< ALUSimplexGrid< 3, 3 > >
  {
    static const GeometryType :: BasicType basicType = GeometryType :: simplex;
    static const int dimension = 3;
    static const bool correctJacobian = true;

    typedef GenericGeometry :: Convert< basicType, dimension > Convert;
    typedef Convert :: type Type;
  };

  template<>
  struct Topology< ALUSimplexGrid< 2, 2 > >
  {
    static const GeometryType :: BasicType basicType = GeometryType :: simplex;
    static const int dimension = 2;
    static const bool correctJacobian = true;

    typedef GenericGeometry :: Convert< basicType, dimension > Convert;
    typedef Convert :: type Type;
  };



  // GeometryTraits
  // --------------

  namespace GenericGeometry
  {

    template< class Traits, class DuneGeometry >
    struct DuneCache
      : public GenericGeometry :: ComputeAll< Traits >
    {
      DuneCache ()
      {}

      template< class GeometryType >
      DuneCache( const GeometryType & )
      {}
    };

    template< class CoordTraits, int dimG, int dimW, class DuneGeometry >
    struct DuneCache< MappingTraits< CoordTraits, dimG, dimW >, DuneGeometry >
    {
      typedef MappingTraits< CoordTraits, dimG, dimW > Traits;
      typedef DuneGeometry DuneGeometryType;

      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
      static const EvaluationType evaluateNormal = ComputeOnDemand;

    private:
      const DuneGeometryType &geo_;

    public:
      DuneCache ( const DuneGeometryType &geo )
        : geo_( geo )
      {}

      void jacobianT ( typename Traits :: JacobianTransposedType &jT ) const
      {}

      void integrationElement( typename Traits :: FieldType &intEl ) const
      {
        FieldVector< double, DuneGeometryType :: mydimension > x( 0 );
        intEl = geo_.integrationElement( x );
      }

      void jacobianInverseTransposed ( typename Traits :: JacobianType &jTInv ) const
      {}

      void normal ( int face, typename Traits :: GlobalCoordType &normal ) const
      {}
    };



    template< class Grid >
    struct GlobalGeometryTraits
    {
      typedef typename Grid :: template Codim< 0 > :: Geometry DuneGeometryType;

      typedef DuneCoordTraits< typename DuneGeometryType :: ctype > CoordTraits;

      static const int dimGrid = DuneGeometryType :: dimension;
      static const int dimWorld = DuneGeometryType :: coorddimension;

      static const bool hybrid = true;

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
        : public DuneCache< Traits, DuneGeometryType >
      {};
    };

  }

}

#endif
