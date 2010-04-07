// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH
#define DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/cachedmapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // CachedMappingFactory
    // --------------------

    template< class Topology, class GeometryTraits >
    class CachedMappingFactory
    {
      typedef CachedMappingFactory< Topology, GeometryTraits > This;

    public:
      typedef CachedMapping< Topology, GeometryTraits > Mapping;
      typedef typename GeometryTraits::Allocator Allocator;

      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, Allocator &allocator )
      {
        assert( type.dim() == Mapping::dimension );
        Mapping *mapping = allocator.template allocate< Mapping >();
        allocator.construct( mapping, Mapping( coords ) );
        return mapping;
      }
    };



    // VirtualMappingFactory
    // ---------------------

    template< unsigned int dim, class GeometryTraits >
    class VirtualMappingFactory
    {
      typedef VirtualMappingFactory< dim, GeometryTraits > This;

    public:
      typedef HybridMapping< dim, GeometryTraits > Mapping;
      typedef typename GeometryTraits::Allocator Allocator;

    private:
      template< bool > struct AllTypes;
      template< bool > struct OnlySimplexCube;

      template< GeometryType::BasicType type, class CoordVector >
      static Mapping *virtualMapping ( const CoordVector &coords, Allocator &allocator )
      {
        // Note: for dim < 3, Convert may only be instantiated for cube and simplex
        typedef typename Convert< type, dim >::type Topology;
        typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits > VirtualMapping;

        VirtualMapping *mapping = allocator.template allocate< VirtualMapping >();
        allocator.construct( mapping, VirtualMapping( coords ) );
        return mapping;
      }

    public:
      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, Allocator &allocator )
      {
        assert( type.dim() == Mapping::dimension );
        // This construction instantiates either AllTypes or OnlySimplexCube,
        // depending on Mapping::dimension.
        typedef typename SelectType< (Mapping::dimension >= 3), AllTypes<true>, OnlySimplexCube<false> >::Type Switch;
        return Switch::mapping( type.basicType(), coords, allocator );
      }
    };



    // VirtualMappingFactory::AllTypes
    //
    // Allow construction of VirtualMappings for all geometry types
    // (i.e., simplex, cube, prism and pyramid). This construction
    // method is used if dim >= 3.
    //
    // Note: The boolean template argument is a dummy that prevents
    //       instatiation unless this template is used.

    template< unsigned int dim, class GeometryTraits >
    template< bool >
    struct VirtualMappingFactory< dim, GeometryTraits > :: AllTypes
    {
      template< class CoordVector >
      static Mapping *
      mapping ( GeometryType::BasicType type, const CoordVector &coords, Allocator &allocator )
      {
        switch( type )
        {
        case GeometryType::simplex :
          return virtualMapping< GeometryType::simplex, CoordVector >( coords, allocator );

        case GeometryType::cube :
          return virtualMapping< GeometryType::cube, CoordVector >( coords, allocator );

        case GeometryType::prism :
          return virtualMapping< GeometryType::prism, CoordVector >( coords, allocator );

        case GeometryType::pyramid :
          return virtualMapping< GeometryType::pyramid, CoordVector >( coords, allocator );

        default :
          DUNE_THROW( RangeError, "Unknown basic geometry type: " << type );
        }
      }
    };


    // VirtualMappingFactory::OnlySimplexCube
    //
    // Allow construction of VirtualMappings only for simplices and cube. This
    // construction method is used for dim < 3.
    //
    // Note: The boolean template argument is a dummy that prevents
    //       instatiation unless this template is used.

    template< unsigned int dim, class GeometryTraits >
    template< bool >
    struct VirtualMappingFactory< dim, GeometryTraits > :: OnlySimplexCube
    {
      template< class CoordVector >
      static Mapping *
      mapping ( GeometryType::BasicType type, const CoordVector &coords, Allocator &allocator )
      {
        switch( type )
        {
        case GeometryType::simplex :
          return virtualMapping< GeometryType::simplex, CoordVector >( coords, allocator );

        case GeometryType::cube :
          return virtualMapping< GeometryType::cube, CoordVector >( coords, allocator );

        default :
          DUNE_THROW( RangeError, "Unknown basic geometry type: " << type );
        }
      }
    };




    // MappingProvider
    // ---------------

    template< class ElementMapping, unsigned int codim >
    class MappingProvider;


    template< unsigned int dim, class GeometryTraits, unsigned int codim >
    class MappingProvider< HybridMapping< dim, GeometryTraits >, codim >
    {
      typedef MappingProvider< HybridMapping< dim, GeometryTraits >, codim > This;

    public:
      static const unsigned int dimension = dim;
      static const unsigned int codimension = codim;
      static const unsigned int mydimension = dimension - codimension;

      typedef typename GeometryTraits::Allocator Allocator;

    private:
      typedef VirtualMappingFactory< mydimension, GeometryTraits > Factory;

    public:
      typedef typename Factory :: Mapping Mapping;

      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, Allocator &allocator )
      {
        return Factory::mapping( type, coords, allocator );
      }
    };


    template< class Topology, class GeometryTraits, unsigned int codim >
    class MappingProvider< CachedMapping< Topology, GeometryTraits >, codim >
    {
      typedef MappingProvider< CachedMapping< Topology, GeometryTraits >, codim > This;

    public:
      static const unsigned int dimension = Topology :: dimension;
      static const unsigned int codimension = codim;
      static const unsigned int mydimension = dimension - codimension;

      static const bool hybrid = IsCodimHybrid< Topology, codim > :: value;

      typedef typename GeometryTraits::Allocator Allocator;

    private:
      template< bool >
      struct HybridFactory
        : public VirtualMappingFactory< mydimension, GeometryTraits >
      {};

      template< bool >
      struct NonHybridFactory
        : public CachedMappingFactory
          < typename SubTopology< Topology, codim, 0 > :: type, GeometryTraits >
      {};

      typedef typename SelectType< hybrid, HybridFactory<true>, NonHybridFactory<false> >::Type Factory;

    public:
      typedef typename Factory :: Mapping Mapping;

      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, Allocator &allocator )
      {
        return Factory::mapping( type, coords, allocator );
      }
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH
