// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBMAPPING_HH
#define DUNE_GENERICGEOMETRY_SUBMAPPING_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/matrix.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // External Forward Declarations
    // -----------------------------

    template< class Mapping, unsigned int codim >
    class SubMappingCoords;

    template< class Topology, class GeometricMappingTraits >
    class CachedMapping;

    template< unsigned int DimG, class GeometricMappingTraits >
    class HybridMapping;

    template< class Topology, class GeometricMappingTraits >
    class VirtualMapping;



    // MappingTraits
    // -------------

    template< unsigned int DimG, class CoordTraits, bool alwaysAffine = false >
    struct MappingTraits
    {
      static const unsigned int dimG = DimG;
      static const unsigned int dimW = CoordTraits :: dimWorld;

      static const bool affine = alwaysAffine;

      typedef typename CoordTraits :: ctype FieldType;
      typedef typename CoordTraits :: template Vector< dimG > :: type LocalCoordType;
      typedef typename CoordTraits :: template Vector< dimW > :: type GlobalCoordType;

      typedef typename CoordTraits :: template Matrix< dimW, dimG > :: type
      JacobianType;
      typedef typename CoordTraits :: template Matrix< dimG, dimW > :: type
      JacobianTransposedType;

      typedef GenericGeometry :: MatrixHelper< CoordTraits > MatrixHelper;
    };



    // SubMappingTraits
    // ----------------

    template< class Mapping, unsigned int codim >
    struct SubMappingTraits;

    template< unsigned int DimG, class GeometricMappingTraits, unsigned int codim >
    struct SubMappingTraits< HybridMapping< DimG, GeometricMappingTraits >, codim >
    {
      enum { dimension = DimG };
      enum { isVirtual = true };

      typedef typename GeometricMappingTraits :: template Traits< dimension - codim >
      MappingTraits;
      typedef typename MappingTraits :: CachingType CachingType;

      typedef HybridMapping< dimension - codim, GeometricMappingTraits >
      HybridSubMapping;

      template< GeometryType :: BasicType btype >
      struct VirtualMapping
      {
        typedef typename Convert< btype, dimension - codim > :: type SubTopology;
        typedef GenericGeometry :: VirtualMapping< SubTopology, GeometricMappingTraits > type;
      };

      typedef HybridSubMapping SubMapping;
    };

    template< class Topology, class GeometricMappingTraits, unsigned int codim >
    struct SubMappingTraits< VirtualMapping< Topology, GeometricMappingTraits >, codim >
    {
      enum { dimension = Topology :: dimension };
      enum { isVirtual = true };

      typedef typename GeometricMappingTraits :: template Traits< dimension - codim >
      MappingTraits;
      typedef typename MappingTraits :: CachingType CachingType;

      typedef HybridMapping< dimension - codim, GeometricMappingTraits > HybridSubMapping;

      template< unsigned int i >
      struct VirtualSubMapping
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
        SubTopology;
        typedef VirtualMapping< SubTopology, GeometricMappingTraits > type;
      };

      template< GeometryType :: BasicType btype >
      struct VirtualMapping
      {
        typedef typename Convert< btype, dimension - codim > :: type SubTopology;
        typedef GenericGeometry :: VirtualMapping< SubTopology, GeometricMappingTraits > type;
      };

      typedef HybridSubMapping SubMapping;
    };

    template< class Topology, class GeometricMappingTraits, unsigned int codim >
    struct SubMappingTraits< CachedMapping< Topology, GeometricMappingTraits >, codim >
    {
      enum { dimension = Topology :: dimension };
      enum { isVirtual = IsCodimHybrid< Topology, codim > :: value };

      typedef typename GeometricMappingTraits :: CoordinateTraits CoordTraits;

      typedef typename GeometricMappingTraits :: template Traits< dimension - codim >
      MappingTraits;
      typedef typename MappingTraits :: CachingType CachingType;

      typedef HybridMapping< dimension - codim, GeometricMappingTraits > HybridSubMapping;

      template< unsigned int i >
      struct VirtualSubMapping
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
        SubTopology;
        typedef VirtualMapping< SubTopology, GeometricMappingTraits > type;
      };

      template< GeometryType :: BasicType btype >
      struct VirtualMapping
      {
        typedef typename Convert< btype, dimension - codim > :: type SubTopology;
        typedef GenericGeometry :: VirtualMapping< SubTopology, GeometricMappingTraits > type;
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
        typedef CachedMapping< SubTopology, GeometricMappingTraits > SubMapping;
      };

      typedef typename ProtectedIf< isVirtual, Hybrid, NonHybrid > :: SubMapping
      SubMapping;
    };



    // HybridSubMappingProvider
    // ------------------------

    template< class Mapping, unsigned int codim >
    class HybridSubMappingProvider
    {
      typedef typename Mapping :: ReferenceElement ReferenceElement;
      enum { numSubMappings = ReferenceElement :: template Codim< codim > :: size };

      struct CreatorInterface;
      template< int i > struct CreatorImplementation;

      typedef const CreatorInterface *CreatorPtr;
      CreatorPtr creator_[ numSubMappings ];

      typedef GenericGeometry :: SubMappingCoords< Mapping, codim > SubMappingCoords;

    public:
      typedef typename SubMappingTraits< Mapping, codim > :: HybridSubMapping SubMapping;

      static SubMapping *
      subMapping ( const Mapping &mapping, unsigned int i,
                   const typename SubMapping :: CachingType &cache )
      {
        assert( i < numSubMappings );
        SubMappingCoords coords( mapping, i );
        return instance().creator_[ i ]->create( coords, cache );
      }

    private:
      HybridSubMappingProvider ()
      {
        ForLoop< CreatorImplementation, 0, numSubMappings-1 > :: apply( creator_ );
      }

      static const HybridSubMappingProvider &instance ()
      {
        static HybridSubMappingProvider inst;
        return inst;
      }
    };

    template< class Mapping, unsigned int codim >
    struct HybridSubMappingProvider< Mapping, codim > :: CreatorInterface
    {
      virtual SubMapping *
      create ( const SubMappingCoords &coords,
               const typename SubMapping :: CachingType &cache ) const = 0;
      virtual ~CreatorInterface() {}
    };

    template< class Mapping, unsigned int codim >
    template< int i >
    struct HybridSubMappingProvider< Mapping, codim > :: CreatorImplementation
      : public CreatorInterface
    {
      typedef typename SubMappingTraits< Mapping, codim >
      :: template VirtualSubMapping< (unsigned int) i > :: type
      VirtualSubMapping;

      virtual SubMapping *
      create ( const SubMappingCoords &coords,
               const typename SubMapping :: CachingType &cache ) const
      {
        return new VirtualSubMapping( coords, cache );
      }

      static void apply ( CreatorPtr (&creator)[ numSubMappings ] )
      {
        creator[ i ] = new CreatorImplementation;
      }
    };



    // SubMappingProvider
    // ------------------

    template< class Mapping, unsigned int codim >
    class SubMappingProvider
    {
      typedef typename Mapping :: ReferenceElement ReferenceElement;
      enum { numSubMappings = ReferenceElement :: template Codim< codim > :: size };

      static const bool isVirtual = SubMappingTraits< Mapping, codim > :: isVirtual;

      typedef GenericGeometry :: HybridSubMappingProvider< Mapping, codim >
      HybridSubMappingProvider;

    public:
      typedef typename SubMappingTraits< Mapping, codim > :: SubMapping SubMapping;

    private:
      template< bool >
      struct Virtual
      {
        static SubMapping *
        subMapping ( const Mapping &mapping, unsigned int i,
                     const typename SubMapping :: CachingType &cache )
        {
          return HybridSubMappingProvider :: subMapping( mapping, i, cache );
        }
      };

      template< bool >
      struct NonVirtual
      {
        static SubMapping *
        subMapping ( const Mapping &mapping, unsigned int i,
                     const typename SubMapping :: CachingType &cache )
        {
          assert( i < numSubMappings );
          SubMappingCoords< Mapping, codim > coords( mapping, i );
          return new SubMapping( coords, cache );
        }
      };

    public:
      static SubMapping *
      subMapping ( const Mapping &mapping, unsigned int i,
                   const typename SubMapping :: CachingType &cache )
      {
        return ProtectedIf< isVirtual, Virtual, NonVirtual > :: subMapping( mapping, i, cache );
      }
    };

  }

}

#endif
