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

      template< bool > struct NonVirtual;
      template< bool > struct Virtual;

      template< GeometryType :: BasicType type, class CoordVector >
      static Mapping *virtualMapping ( const CoordVector &coords )
      {
        typedef typename SubMappingTraits :: template VirtualMapping< type > :: type
        VirtualMapping;
        return new VirtualMapping( coords );
      }

      template< class CoordVector >
      static Mapping *mapping ( const GeometryType &type, const CoordVector &coords )
      {
        typedef ProtectedIf< SubMappingTraits :: isVirtual, Virtual, NonVirtual > Switch;
        return Switch :: mapping( type, coords );
      }
    };


    template< class ElementMapping, unsigned int codim >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: NonVirtual
    {
      template< class CoordVector >
      static Mapping *
      mapping ( const GeometryType &type, const CoordVector &coords )
      {
        assert( type.dim() == Mapping :: dimG );
        return new Mapping( coords );
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
      mapping ( const GeometryType &type, const CoordVector &coords )
      {
        assert( type.dim() == Mapping :: dimG );
        typedef ProtectedIf< (Mapping :: dimG >= 3), AllTypes, OnlySimplexCube > Switch;
        return Switch :: mapping( type.basicType(), coords );
      }
    };


    template< class ElementMapping, unsigned int codim >
    template< bool b >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: Virtual< b > :: AllTypes
    {
      template< class CoordVector >
      static Mapping *
      mapping ( GeometryType :: BasicType type, const CoordVector &coords )
      {
        switch( type )
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
    };


    template< class ElementMapping, unsigned int codim >
    template< bool b >
    template< bool >
    struct MappingProvider< ElementMapping, codim > :: Virtual< b > :: OnlySimplexCube
    {
      template< class CoordVector >
      static Mapping *
      mapping ( GeometryType :: BasicType type, const CoordVector &coords )
      {
        switch( type )
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

  }

}

#endif
