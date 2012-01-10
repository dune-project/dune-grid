// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_HOSTCORNERS_HH
#define DUNE_GEOGRID_HOSTCORNERS_HH

#include <dune/geometry/type.hh>

#include <dune/grid/alugrid/3d/topology.hh>
#include <dune/grid/common/entity.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, int, class >
  class ALU3dGridEntity;

  template< ALU3dGridElementType, class >
  struct ALU3dImplTraits;



  namespace GeoGrid
  {

    // HostCorners
    // -----------

    template< class HostEntity >
    class HostCorners
    {
      typedef typename HostEntity::Geometry HostGeometry;

    public:
      typedef typename HostGeometry::GlobalCoordinate Coordinate;

      explicit HostCorners ( const HostEntity &hostEntity )
        : hostGeometry_( hostEntity.geometry() )
      {}

      GeometryType type () const
      {
        return hostGeometry_.type();
      }

      Coordinate corner ( const int i ) const
      {
        return hostGeometry_.corner( i );
      }

      unsigned int numCorners () const
      {
        return hostGeometry_.corners();
      }

    private:
      HostGeometry hostGeometry_;
    };



    // HostCorners for ALU3dGrid
    // -------------------------

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    template< int dim, class Grid >
    class HostCorners< Dune::Entity< 0, dim, Grid, ALU3dGridEntity > >
    {
      typedef Dune::Entity< 0, dim, Grid, ALU3dGridEntity > HostEntity;

      typedef double ALUCoordinate[ 3 ];

      static const ALU3dGridElementType elementType = remove_const< Grid >::type::elementType;
      typedef Dune::ElementTopologyMapping< elementType > ElementTopologyMapping;

      typedef typename remove_const< Grid >::type::MPICommunicatorType Comm;
      typedef ALU3dImplTraits< elementType, Comm > ImplTraits;

    public:
      typedef FieldVector< double, 3 > Coordinate;

      explicit HostCorners ( const HostEntity &hostEntity )
        : item_( hostEntity.impl().getItem() )
      {}

      GeometryType type () const
      {
        if( elementType == tetra )
          return GeometryType( GenericGeometry::SimplexTopology< dim >::type::id, dim );
        else
          return GeometryType( GenericGeometry::CubeTopology< dim >::type::id, dim );
      }

      Coordinate corner ( const int i ) const
      {
        const int j = ElementTopologyMapping::dune2aluVertex( i );
        const ALUCoordinate &point = item_.myvertex( j )->Point();

        Coordinate corner;
        for( int k = 0; k < 3; ++k )
          corner[ k ] = point[ k ];
        return corner;
      }

      unsigned int numCorners () const
      {
        return (elementType == tetra ? dim+1 : (1 << dim));
      }

    private:
      const typename ImplTraits::IMPLElementType &item_;
    };
#endif // #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_HOSTCORNERS_HH
