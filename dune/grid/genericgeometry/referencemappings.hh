// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEMAPPINGS_HH
#define DUNE_GENERICGEOMETRY_REFERENCEMAPPINGS_HH

#include <dune/grid/genericgeometry/referencetopologies.hh>
#include <dune/grid/genericgeometry/referencedomain.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // ReferenceMappingsContainer
    // --------------------------

    template< class ctype, unsigned int dim >
    struct ReferenceMappingsContainer
    {
      static const unsigned int dimension = dim;

      typedef GenericGeometry::ReferenceTopology< dimension > ReferenceTopology;

    private:
      template< class Topology > class CornerStorage;
      template< class Topology > struct Initialize;
      template< int codim > struct Finalize;

      struct GeometryTraits;

    public:
      template< int codim >
      struct Codim
      {
        typedef HybridMapping< dimension-codim, GeometryTraits > Mapping;
      };

      ~ReferenceMappingsContainer ()
      {
        ForLoop< Finalize, 0, dimension >::apply( mappings_, allocator_ );
      }

      const ReferenceTopology &topology () const
      {
        return *topology_;
      }

      template< int codim >
      typename Codim< codim >::Mapping &mapping ( const unsigned int i ) const
      {
        Int2Type< codim > codimVariable;
        return *(mappings_[ codimVariable ][ i ]);
      }

      template< class Topology >
      void initialize ()
      {
        typedef Initialize< Topology > Init;
        typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits > VirtualMapping;

        Int2Type< 0 > codim0Variable;
        topology_ = &(ReferenceTopologies< dimension >::get( Topology::id ));

        VirtualMapping *virtualMapping = allocator_.template allocate< VirtualMapping >();
        allocator_.construct( virtualMapping, VirtualMapping( codim0Variable ) );
        virtualMapping->referenceCount = 1;

        mappings_[ codim0Variable ].resize( 1 );
        mappings_[ codim0Variable ][ 0 ] = virtualMapping;

        ForLoop< Init::template Codim, 1, dimension >::apply( mappings_, allocator_ );
      }

    private:
      template< int codim >
      class MappingArray
        : public std::vector< typename Codim< codim >::Mapping * >
      {};

      typedef CodimTable< MappingArray, dimension > MappingsTable;

      typename GeometryTraits::Allocator allocator_;
      const ReferenceTopology *topology_;
      MappingsTable mappings_;
    };


    template< class ctype, unsigned int dim >
    template< class Topology >
    class ReferenceMappingsContainer< ctype, dim >::CornerStorage
    {
      typedef ReferenceDomain< Topology > RefDomain;

    public:
      static const unsigned int size = Topology::numCorners;

      template< class SubTopology >
      struct SubStorage
      {
        typedef CornerStorage< SubTopology > type;
      };

      explicit CornerStorage ( const Int2Type< 0 > & )
      {
        for( unsigned int i = 0; i < size; ++i )
          RefDomain::template corner< ctype >( i, coords_[ i ] );
      }

      template< class Mapping, unsigned int codim >
      explicit
      CornerStorage ( const SubMappingCoords< Mapping, codim > &coords )
      {
        for( unsigned int i = 0; i < size; ++i )
          coords_[ i ] = coords[ i ];
      }

      const FieldVector< ctype, dim > &operator[] ( unsigned int i ) const
      {
        return coords_[ i ];
      }

    private:
      FieldVector< ctype, dim > coords_[ size ];
    };


    template< class ctype, unsigned int dim >
    template< class Topology >
    struct ReferenceMappingsContainer< ctype, dim >::Initialize
    {
      typedef typename ReferenceMappingsContainer< ctype, dim >::template Codim< 0 >::Mapping ReferenceMapping;

      template< int codim >
      struct Codim
      {
        static void apply ( MappingsTable &mappings, typename GeometryTraits::Allocator &allocator )
        {
          const unsigned int size = GenericGeometry::Size< Topology, codim >::value;

          Int2Type< 0 > codim0Variable;
          const ReferenceMapping &refMapping = *(mappings[ codim0Variable ][ 0 ]);

          Int2Type< codim > codimVariable;
          mappings[ codimVariable ].resize( size );
          for( unsigned int i = 0; i < size; ++i )
          {
            mappings[ codimVariable ][ i ] = refMapping.template trace< codim >( i, allocator );
            mappings[ codimVariable ][ i ]->referenceCount = 1;
          }
        }
      };
    };


    template< class ctype, unsigned int dim >
    template< int codim >
    struct ReferenceMappingsContainer< ctype, dim >::Finalize
    {
      static void apply ( MappingsTable &mappings, typename GeometryTraits::Allocator &allocator )
      {
        Int2Type< codim > codimVariable;
        const unsigned int size = mappings[ codimVariable ].size();
        for( unsigned int i = 0; i < size; ++i )
        {
          allocator.destroy( mappings[ codimVariable ][ i ] );
          allocator.deallocate( mappings[ codimVariable ][ i ] );
        }
      };
    };


    template< class ctype, unsigned int dim >
    struct ReferenceMappingsContainer< ctype, dim >::GeometryTraits
      : public GenericGeometry::DefaultGeometryTraits< ctype, dimension, dimension >
    {
      typedef GenericGeometry::DefaultGeometryTraits< ctype, dimension, dimension > Base;

      typedef typename Base::CoordTraits CoordTraits;

      template< class Topology >
      struct Mapping
      {
        typedef GenericGeometry::CornerMapping< CoordTraits, Topology, dimension, CornerStorage< Topology >, true > type;
      };

      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = PreCompute;
        static const EvaluationType evaluateJacobianInverseTransposed = PreCompute;
        static const EvaluationType evaluateIntegrationElement = PreCompute;
        static const EvaluationType evaluateNormal = PreCompute;
      };
    };



    // ReferenceMappings
    // -----------------

    template< class ctype, unsigned int dim >
    struct ReferenceMappings
    {
      static const unsigned int dimension = dim;
      static const unsigned int numTopologies = (1 << dimension);

      typedef ReferenceMappingsContainer< ctype, dimension > Container;

      static const Container &container ( const unsigned int topologyId )
      {
        static ReferenceMappings instance;
        assert( topologyId < numTopologies );
        return instance.containers_[ topologyId ];
      }

    private:
      template< int topologyId >
      struct Initialize
      {
        typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
        static void apply ( Container (&containers)[ numTopologies ]  )
        {
          containers[ topologyId ].template initialize< Topology >();
        }
      };

      ReferenceMappings ()
      {
        ForLoop< Initialize, 0, numTopologies-1 >::apply( containers_ );
      }

      Container containers_[ numTopologies ];
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_REFERENCEMAPPINGS_HH
