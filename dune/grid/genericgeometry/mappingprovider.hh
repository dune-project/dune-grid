// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH
#define DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/genericgeometry/maximum.hh>
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

      static const unsigned int mappingSize = sizeof( Mapping );

      template< class CoordVector >
      static Mapping *
      mapping ( const unsigned int topologyId, const CoordVector &coords, Allocator &allocator )
      {
        assert( topologyId == Topology::id );
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

      static const unsigned int numTopologies = (1 << dim);

      template< class CoordVector >
      class ConstructorTable;

      template< int topologyId >
      struct MappingSize
      {
        typedef typename GenericGeometry::Topology< (unsigned int) topologyId, dim >::type Topology;
        static const int v = sizeof( VirtualMapping< Topology, GeometryTraits > );
      };

    public:
      typedef HybridMapping< dim, GeometryTraits > Mapping;
      typedef typename GeometryTraits::Allocator Allocator;

      static const unsigned int mappingSize = Maximum< MappingSize, 0, numTopologies-1 >::v;

      template< class CoordVector >
      static Mapping *
      mapping ( const unsigned int topologyId, const CoordVector &coords, Allocator &allocator )
      {
        static ConstructorTable< CoordVector > construct;
        return construct[ topologyId ]( coords, allocator );
      }
    };


    // VirtualMappingFactory::ConstructorTable
    // ---------------------------------------

    template< unsigned int dim, class GeometryTraits >
    template< class CoordVector >
    class VirtualMappingFactory< dim, GeometryTraits >::ConstructorTable
    {
      typedef Mapping *(*Construct)( const CoordVector &coords, Allocator &allocator );

      template< int i >
      struct Builder;

    public:
      ConstructorTable ()
      {
        ForLoop< Builder, 0, numTopologies-1 >::apply( construct_ );
      }

      Construct operator[] ( const unsigned int topologyId )
      {
        assert( topologyId < numTopologies );
        return construct_[ topologyId ];
      }

    private:
      template< class Topology >
      static Mapping *construct ( const CoordVector &coords, Allocator &allocator )
      {
        typedef VirtualMapping< Topology, GeometryTraits > VMapping;
        VMapping *mapping = allocator.template allocate< VMapping >();
        allocator.construct( mapping, VMapping( coords ) );
        return mapping;
      }

      Construct construct_[ numTopologies ];
    };


    // VirtualMappingFactory::ConstructorTable::Builder
    // ------------------------------------------------

    template< unsigned int dim, class GeometryTraits >
    template< class CoordVector >
    template< int topologyId >
    struct VirtualMappingFactory< dim, GeometryTraits >::ConstructorTable< CoordVector >::Builder
    {
      static void apply ( Construct (&construct)[ numTopologies ] )
      {
        typedef typename GenericGeometry::Topology< (unsigned int) topologyId, dim >::type Topology;
        construct[ topologyId ] = ConstructorTable< CoordVector >::template construct< Topology >;
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
      // Maximal amount of memory required to store a mapping
      static const unsigned int mappingSize = Factory::mappingSize;

      typedef typename Factory::Mapping Mapping;

      template< class CoordVector >
      static Mapping *
      mapping ( const unsigned int topologyId, const CoordVector &coords, Allocator &allocator )
      {
        return Factory::mapping( topologyId, coords, allocator );
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
          < typename SubTopology< Topology, codim, 0 >::type, GeometryTraits >
      {};

      typedef typename SelectType< hybrid, HybridFactory<true>, NonHybridFactory<false> >::Type Factory;

    public:
      // Maximal amount of memory required to store a mapping
      static const unsigned int mappingSize = Factory::mappingSize;

      typedef typename Factory::Mapping Mapping;

      template< class CoordVector >
      static Mapping *
      mapping ( const unsigned int topologyId, const CoordVector &coords, Allocator &allocator )
      {
        return Factory::mapping( topologyId, coords, allocator );
      }
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH
