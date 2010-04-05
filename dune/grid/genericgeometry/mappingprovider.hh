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

      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords )
      {
        assert( type.dim() == Mapping :: dimension );
        return new Mapping( coords );
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

    private:

      template< GeometryType :: BasicType type, class CoordVector >
      static Mapping *virtualMapping ( const CoordVector &coords )
      {
        typedef typename Convert< type, dim > :: type Topology;
        return new VirtualMapping< Topology, GeometryTraits >( coords );
      }

    public:
      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords )
      {
        assert( type.dim() == Mapping :: dimension );

        if (Mapping :: dimension >= 3) {  // static if: can be evaluated at compile time

          // If our dimension is at least 3 there are four possible types of elements
          switch( type.basicType() )
          {
          case GeometryType :: simplex :
            return virtualMapping< GeometryType :: simplex, CoordVector >( coords );

          case GeometryType :: cube :
            return virtualMapping< GeometryType :: cube, CoordVector >( coords );

          case GeometryType :: prism :
            return virtualMapping< GeometryType :: prism, CoordVector >( coords );

          case GeometryType :: pyramid :
            return virtualMapping< GeometryType :: pyramid, CoordVector >( coords );

          default :
            DUNE_THROW( RangeError, "Unknown basic geometry type: " << type );
          }

        }

        // Our dimension is less than 3.  There can be only triangles or quadrilaterals.
        /** \todo The case dim=1 is actually missing */
        switch( type.basicType() )
        {
        case GeometryType :: simplex :
          return virtualMapping< GeometryType :: simplex, CoordVector >( coords );

        case GeometryType :: cube :
          return virtualMapping< GeometryType :: cube, CoordVector >( coords );

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

    private:
      typedef VirtualMappingFactory< mydimension, GeometryTraits > Factory;

    public:
      typedef typename Factory :: Mapping Mapping;

      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords )
      {
        return Factory :: mapping( type, coords );
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
      mapping ( const GeometryType &type, const CoordVector &coords )
      {
        return Factory :: mapping( type, coords );
      }
    };

  }

}

#endif
