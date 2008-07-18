// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBMAPPING_HH
#define DUNE_GENERICGEOMETRY_SUBMAPPING_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // External Forward Declarations
    // -----------------------------

    template< class Topology, class CoordTraits, template< class > class Caching >
    class CachedMapping;

    template< unsigned int DimG, class CoordTraits, template< class > class Caching >
    class HybridMapping;

    template< class Topology, class CoordTraits, template< class > class Caching >
    class VirtualMapping;



    // MappingTraits
    // -------------

    template< unsigned int DimG, class CoordTraits >
    struct MappingTraits
    {
      enum { dimG = DimG };
      enum { dimW = CoordTraits :: dimCoord };
      enum { affine = CoordTraits :: affine };

      typedef typename CoordTraits :: FieldType FieldType;
      typedef typename CoordTraits :: template Vector< dimG > :: Type LocalCoordType;
      typedef typename CoordTraits :: template Vector< dimW > :: Type GlobalCoordType;

      typedef typename CoordTraits :: template Matrix< dimW, dimG > :: Type
      JacobianType;
      typedef typename CoordTraits :: template Matrix< dimG, dimW > :: Type
      JacobianTransposedType;
    };



    // CachedMappingTraits
    // -------------------

    template< unsigned int DimG, class CoordTraits, template< class > class Caching >
    struct CachedMappingTraits
      : public GenericGeometry :: MappingTraits< DimG, CoordTraits >
    {
      typedef MappingTraits< DimG, CoordTraits > BaseType;

      typedef Caching< BaseType > CachingType;
    };



    // SubMappingCoords
    // ----------------

    template< class Mapping, unsigned int codim >
    class SubMappingCoords
    {
      typedef typename Mapping :: GlobalCoordType GlobalCoordType;
      typedef typename Mapping :: ReferenceElement ReferenceElement;

      enum { dimension = ReferenceElement :: dimension };

      const Mapping &mapping_;
      const unsigned int i_;

    public:
      SubMappingCoords ( const Mapping &mapping, unsigned int i )
        : mapping_( mapping ), i_( i )
      {}

      const GlobalCoordType &operator[] ( unsigned int j ) const
      {
        const unsigned int k
          = ReferenceElement :: template subNumbering< codim, dimension - codim >( i_, j );
        return mapping_[ k ];
      }
    };



    // SubMappingTraits
    // ----------------

    template< class Mapping, unsigned int codim >
    struct SubMappingTraits;

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct SubMappingTraits< VirtualMapping< Topology, CoordTraits, Caching >, codim >
    {
      enum { dimension = Topology :: dimension };
      enum { isVirtual = true };

      typedef CachedMappingTraits< dimension - codim, CoordTraits, Caching > MappingTraits;
      typedef typename MappingTraits :: CachingType CachingType;

      typedef HybridMapping< dimension - codim, CoordTraits, Caching > HybridSubMapping;

      template< unsigned int i >
      struct VirtualSubMapping
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
        SubTopology;
        typedef VirtualMapping< SubTopology, CoordTraits, Caching > type;
      };

      typedef HybridSubMapping SubMapping;
    };

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct SubMappingTraits< CachedMapping< Topology, CoordTraits, Caching >, codim >
    {
      enum { dimension = Topology :: dimension };
      enum { isVirtual = IsCodimHybrid< Topology, codim > :: value };

      typedef CachedMappingTraits< dimension - codim, CoordTraits, Caching > MappingTraits;
      typedef typename MappingTraits :: CachingType CachingType;

      typedef HybridMapping< dimension - codim, CoordTraits, Caching > HybridSubMapping;

      template< unsigned int i >
      struct VirtualSubMapping
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
        SubTopology;
        typedef VirtualMapping< SubTopology, CoordTraits, Caching > type;
      };

      template< bool >
      struct Hybrid
      {
        typedef HybridSubMapping SubMapping;
      };

      template< bool >
      struct NonHybrid
      {
        typedef typename VirtualSubMapping< 0 > :: SubTopology SubTopology;
        typedef CachedMapping< SubTopology, CoordTraits, Caching > SubMapping;
      };

      typedef typename ProtectedIf< isVirtual, Hybrid, NonHybrid > :: SubMapping
      SubMapping;
    };

  }

}

#endif
