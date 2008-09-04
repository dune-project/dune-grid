// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH
#define DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/genericgeometry/submapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // External Forward Declarations
    // -----------------------------

    template< class Topology, class CoordTraits, template< class > class Caching >
    class CachedMapping;


    // HybridMapping
    // -------------

    template< unsigned int DimG, class CoordTraits, template< class > class Caching >
    class HybridMapping
      : public SmallObject
    {
      typedef HybridMapping< DimG, CoordTraits, Caching > ThisType;

    protected:
      typedef CachedMappingTraits< DimG, CoordTraits, Caching > Traits;

    public:
      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      typedef typename Traits :: CachingType CachingType;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename SubMappingTraits< ThisType, codim > :: SubMapping SubMapping;
        typedef typename SubMappingTraits< ThisType, codim > :: CachingType CachingType;
      };

      unsigned int referenceCount;

      virtual ~HybridMapping ()
      {}

      virtual unsigned int topologyId () const = 0;

      virtual int corners () const = 0;

      virtual const GlobalCoordType &operator[] ( int i ) const = 0;

      virtual GlobalCoordType global ( const LocalCoordType &local ) const = 0;

      virtual LocalCoordType local ( const GlobalCoordType &global ) const = 0;

      virtual bool checkInside ( const LocalCoordType &local ) const = 0;

      virtual bool affine () const = 0;

      virtual FieldType integrationElement ( const LocalCoordType &local ) const = 0;

      virtual FieldType volume () const = 0;

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &local ) const = 0;

      virtual GlobalCoordType
      normal ( int face, const LocalCoordType &local ) const = 0;

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *
      subMapping ( unsigned int i,
                   const typename Codim< codim > :: CachingType &cache ) const
      {
        typedef typename Codim< codim > :: SubMapping SubMapping;
        return reinterpret_cast< SubMapping * >( subMapping( codim, i, &cache ) );
      }

    protected:
      virtual void *
      subMapping ( unsigned int codim, unsigned int i, const void *cache ) const = 0;
    };



    template< class Topology, class CoordTraits, template< class > class Caching >
    class VirtualMapping
      : public HybridMapping< Topology :: dimension, CoordTraits, Caching >
    {
      typedef HybridMapping< Topology :: dimension, CoordTraits, Caching > BaseType;
      typedef VirtualMapping< Topology, CoordTraits, Caching > ThisType;

      typedef typename BaseType :: Traits Traits;

      typedef CachedMapping< Topology, CoordTraits, Caching > Mapping;

    public:
      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      typedef typename Traits :: CachingType CachingType;

      typedef typename Mapping :: ReferenceElement ReferenceElement;

      template< unsigned int codim >
      struct Codim
      {
        typedef typename SubMappingTraits< ThisType, codim > :: SubMapping SubMapping;
        typedef typename SubMappingTraits< ThisType, codim > :: CachingType CachingType;
      };

    private:
      Mapping mapping_;

    public:
      template< class CoordVector >
      explicit VirtualMapping ( const CoordVector &coordVector,
                                const CachingType &cache )
        : mapping_( coordVector, cache )
      {}

      virtual unsigned int topologyId () const
      {
        return mapping_.topologyId();
      }

      virtual int corners () const
      {
        return mapping_.corners();
      }

      virtual const GlobalCoordType &operator[] ( int i ) const
      {
        return mapping_[ i ];
      }

      virtual GlobalCoordType global ( const LocalCoordType &local ) const
      {
        return mapping_.global( local );
      }

      virtual LocalCoordType local ( const GlobalCoordType &global ) const
      {
        return mapping_.local( global );
      }

      virtual bool checkInside ( const LocalCoordType &local ) const
      {
        return mapping_.checkInside( local );
      }

      virtual bool affine () const
      {
        return mapping_.affine();
      }

      virtual FieldType integrationElement ( const LocalCoordType &local ) const
      {
        return mapping_.integrationElement( local );
      }

      virtual FieldType volume () const
      {
        return mapping_.volume();
      }

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &local ) const
      {
        return mapping_.jacobianInverseTransposed( local );
      }

      virtual GlobalCoordType
      normal ( int face, const LocalCoordType &local ) const
      {
        return mapping_.normal( face , local );
      }

      template< unsigned int codim >
      typename Codim< codim > :: SubMapping *
      subMapping ( unsigned int i,
                   const typename Codim< codim > :: CachingType &cache ) const
      {
        return SubMappingProvider< ThisType, codim > :: subMapping( *this, i, cache );
      }

    private:
      class CodimCaller;

    protected:
      void *subMapping ( unsigned int codim, unsigned int i, const void *cache ) const
      {
        return CodimCaller :: subMapping ( *this, codim, i, cache );
      }
    };


    template< class Topology, class CoordTraits, template< class > class Caching >
    class VirtualMapping< Topology, CoordTraits, Caching > :: CodimCaller
    {
      typedef VirtualMapping< Topology, CoordTraits, Caching > Mapping;

      struct CallerInterface;
      template< int codim > struct CallerImplementation;

      CallerInterface *caller_[ dimG+1 ];

      CodimCaller ()
      {
        ForLoop< CallerImplementation, 0, dimG > :: apply( caller_ );
      }

      static const CodimCaller &instance ()
      {
        static CodimCaller inst;
        return inst;
      }

    public:
      static void *
      subMapping ( const Mapping &mapping, unsigned int codim, unsigned int i, const void *cache )
      {
        assert( codim <= dimG );
        return instance().caller_[ codim ]->subMapping( mapping, i, cache );
      }
    };


    template< class Topology, class CoordTraits, template< class > class Caching >
    struct VirtualMapping< Topology, CoordTraits, Caching > :: CodimCaller :: CallerInterface
    {
      typedef VirtualMapping< Topology, CoordTraits, Caching > Mapping;

      virtual ~CallerInterface ()
      {}

      virtual void *
      subMapping ( const Mapping &mapping, unsigned int i, const void *cache ) const = 0;
    };

    template< class Topology, class CoordTraits, template< class > class Caching >
    template< int codim >
    struct VirtualMapping< Topology, CoordTraits, Caching > :: CodimCaller :: CallerImplementation
      : public CallerInterface
    {
      virtual void *
      subMapping ( const Mapping &mapping, unsigned int i, const void *cache ) const
      {
        typedef typename Codim< (unsigned int) codim > :: CachingType CachingType;
        const CachingType &caching = *reinterpret_cast< const CachingType * >( cache );
        return mapping.template subMapping< codim >( i, caching );
      }

      static void apply ( CallerInterface *(&caller)[ dimG+1 ] )
      {
        caller[ codim ] = new CallerImplementation< codim >;
      }
    };

  }

}

#endif
