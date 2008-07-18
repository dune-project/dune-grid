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



    // If not affine only volume is cached (based on intElCompute)
    // otherwise all quantities can be cached using:
    //   geoCompute:    assign if method called using barycenter
    //   geoPreCompute: assign in constructor using barycenter
    //   geoIsComputed: assign in constructor using barycenter using callback
    enum { geoCompute, geoPreCompute, geoIsComputed };

    template <class Traits>
    struct ComputeAll {
      enum {jTCompute = geoCompute,
            jTInvCompute = geoCompute,
            intElCompute = geoCompute,
            normalCompute = geoCompute};
      void jacobianT(typename Traits::JacobianTransposedType& d) const {}
      void integrationElement(typename Traits::FieldType& intEl) const {}
      void jacobianInverseTransposed(typename Traits::JacobianType& dInv) const {}
      void normal(int face, typename Traits::GlobalCoordType& n) const {}
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

    template< unsigned int DimG, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct SubMappingTraits< HybridMapping< DimG, CoordTraits, Caching >, codim >
    {
      enum { dimension = DimG };
      enum { isVirtual = true };

      typedef CachedMappingTraits< dimension - codim, CoordTraits, Caching > MappingTraits;
      typedef typename MappingTraits :: CachingType CachingType;

      typedef HybridMapping< dimension - codim, CoordTraits, Caching > HybridSubMapping;

      typedef HybridSubMapping SubMapping;
    };

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



    // HYbridSubMappingProvider
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

      enum { isVirtual = SubMappingTraits< Mapping, codim > :: isVirtual };

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
