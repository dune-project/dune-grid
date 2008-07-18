// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH
#define DUNE_GENERICGEOMTRY_HYBRIDMAPPING_HH

#include <dune/grid/genericgeometry/mappings.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< unsigned int DimG, class CoordTraits, template< class > class Caching >
    class HybridMapping
    {
      typedef CachedMappingTraits< DimG, CoordTraits, Caching > Traits;

    public:
      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      typedef typename Traits :: CachingType CachingType;

      virtual ~HybridMapping ()
      {}

      virtual GeometryType type () const = 0;

      virtual int corners () const = 0;

      virtual const GlobalCoordType &operator[] ( int i ) const = 0;

      virtual GlobalCoordType global ( const LocalCoordType &local ) const = 0;

      virtual LocalCoordType local ( const GlobalCoordType &global ) const = 0;

      virtual bool checkInside ( const LocalCoordType &local ) const = 0;

      virtual bool affine () const = 0;

      virtual Field integrationElement ( const LocalCoordType &local ) const = 0;

      virtual Field volume () const = 0;

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &local ) const = 0;
    };



    template< class Topology, class CoordTraits, template< class > class Caching >
    class VirtualMapping
      : public HybridMapping< Topology :: dimension, CoordTraits, Caching >
    {
      typedef HybridMapping< Topology :: dimension, CoordTraits, Caching > BaseType;

      typedef typename BaseType :: Traits Traits;

      typedef typename CachedMapping< Topology, Traits, Caching > Mapping;

    public:
      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;

      typedef typename Traits :: CachingType CachingType;

      typedef typename Mapping :: ReferenceElement ReferenceElement;

    private:
      Mapping mapping_;

    public:
      template< class CoordVector >
      explicit VirtualMapping ( const CoordVector &coordVector,
                                const CachingType &cache )
        : mapping_( coordVector, cache )
      {}

      virtual GeometryType type () const
      {
        return mapping_.type();
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

      virtual Field volume () const
      {
        return mapping_.volume();
      }

      virtual const JacobianType &
      jacobianInverseTransposed ( const LocalCoordType &local ) const
      {
        return mapping_.jacobianInverseTransposed( local );
      }
    };

  }

}

#endif
