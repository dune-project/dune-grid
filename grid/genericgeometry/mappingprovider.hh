// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH
#define DUNE_GENERICGEOMETRY_MAPPINGPROVIDER_HH

#include <dune/grid/genericgeometry/submapping.hh>
#include <dune/grid/genericgeometry/cachedmapping.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // MappingProvider
    // ---------------

    template< class ElementMapping, unsigned int codim >
    struct MappingProvider
    {
      typedef GenericGeometry :: SubMappingTraits< ElementMapping, codim > SubMappingTraits;
      typedef typename SubMappingTraits :: SubMapping Mapping;

      typedef typename Mapping :: CachingType CachingType;

      template< bool > struct NonVirtual;
      template< bool > struct Virtual;

      template< GeometryType :: BasicType type, class CoordVector >
      static Mapping *virtualMapping ( const CoordVector &coords,
                                       const CachingType &cache )
      {
        typedef typename SubMappingTraits :: template VirtualMapping< type > :: type
        VirtualMapping;
        return new VirtualMapping( coords, cache );
      }

      template< class CoordVector >
      static Mapping *mapping ( const GeometryType &type,
                                const CoordVector &coords,
                                const CachingType &cache )
      {
        typedef ProtectedIf< SubMappingTraits :: isVirtual, Virtual, NonVirtual > Switch;
        return Switch :: mapping( type, coords, cache );
      }
    };


    template< class ElementMapping, unsigned int codim >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: NonVirtual
    {
      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, const CachingType &cache )
      {
        assert( type.dim() == Mapping :: dimG );
        return new Mapping( coords, cache );
      }
    };


    template< class ElementMapping, unsigned int codim >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: Virtual
    {
      template< bool > struct AllTypes;
      template< bool > struct OnlySimplexCube;

      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords, const CachingType &cache )
      {
        assert( type.dim() == Mapping :: dimG );
        typedef ProtectedIf< (Mapping :: dimG >= 3), AllTypes, OnlySimplexCube > Switch;
        return Switch :: mapping( type.basicType(), coords, cache );
      }
    };


    template< class ElementMapping, unsigned int codim >
    template< bool b >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: Virtual< b > :: AllTypes
    {
      template< class CoordVector >
      static Mapping *
      mapping ( GeometryType :: BasicType type, const CoordVector &coords, const CachingType &cache )
      {
        switch( type )
        {
        case GeometryType :: simplex :
          return virtualMapping< GeometryType :: simplex, CoordVector >( coords, cache );

        case GeometryType :: cube :
          return virtualMapping< GeometryType :: cube, CoordVector >( coords, cache );

        case GeometryType :: prism :
          return virtualMapping< GeometryType :: prism, CoordVector >( coords, cache );

        case GeometryType :: pyramid :
          return virtualMapping< GeometryType :: pyramid, CoordVector >( coords, cache );

        default :
          DUNE_THROW( RangeError, "Unknown basic geometry type: " << type );
        }
      }
    };


    template< class ElementMapping, unsigned int codim >
    template< bool b >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: Virtual< b > :: OnlySimplexCube
    {
      template< class CoordVector >
      static Mapping *
      mapping ( GeometryType :: BasicType type, const CoordVector &coords, const CachingType &cache )
      {
        switch( type )
        {
        case GeometryType :: simplex :
          return virtualMapping< GeometryType :: simplex, CoordVector >( coords, cache );

        case GeometryType :: cube :
          return virtualMapping< GeometryType :: cube, CoordVector >( coords, cache );

        default :
          DUNE_THROW( RangeError, "Unknown basic geometry type: " << type );
        }
      }
    };

  }

}

#endif
