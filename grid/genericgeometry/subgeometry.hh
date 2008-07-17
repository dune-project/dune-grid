// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBGEOMETRY_HH
#define DUNE_GENERICGEOMETRY_SUBGEOMETRY_HH

#include <dune/grid/genericgeometry/hybridgeometries.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, class Traits, template< class > class Caching,
        unsigned int codim >
    class HybridGeometryProvider
    {
      struct CreatorInterface;

      template< int i >
      struct CreatorImplementation;

      typedef const CreatorInterface *CreatorPtr;

      enum { numSubGeometries = Size< Topology, codim > :: value };

      CreatorPtr creator_[ numSubGeometries ];

      enum { dimension = Topology :: dimension - codim };

    public:
      typedef HybridGeometry< dimension, Traits > Geometry;

      static Geometry *
      geometry ( unsigned int i,
                 const typename Geometry :: CoordVector &coords,
                 const typename Geometry :: CachingType &cache )
      {
        assert( i < numSubGeometries );
        return instance().creator_[ i ]->create( coords, cache );
      }

    private:
      HybridGeometryProvider ()
      {
        ForLoop< CreatorImplementation, 0, numSubGeometries-1 > :: apply( creator_ );
      }

      static const HybridGeometryProvider &instance ()
      {
        static HybridGeometryProvider inst;
        return inst;
      }
    };


    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct HybridGeometryProvider< Topology, CoordTraits, Caching, codim >
    :: CreatorInterface
    {
      virtual Geometry *
      create ( const typename Geometry :: CoordVector &coords,
               const typename Geometry :: CachingType &cache ) const = 0;
    };

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    template< int i >
    struct HybridSubGeometryProvider< Topology, CoordTraits, Caching, codim >
    :: CreatorImplementation
      : public CreatorInterface
    {
      typedef typename GenericGeometry
      :: SubTopology< Topoogy, codim, (unsigned int) i > :: Type
      SubTopology;
      typedef GenericGeometry :: VirtualGeometry< SubTopology, CoordTraits, Caching >
      VirtualGeometry;

      virtual Geometry *
      create ( const typename Geometry :: CoordVector &coords
               const typename Geometry :: CachingType &cache ) const
      {
        return new VirtualGeometry( coords, cache )
      }

      static void apply ( CreatorPtr (&creator)[ numSubGeometries ] )
      {
        creator[ i ] = new CreatorImplementation;
      }
    };



    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim, bool isHybrid >
    struct NonHybridGeometryProvider;

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct NonHybridGeometryProvider< Topology, CoordTraits, Caching, codim, true >
      : public HybridGeometryProvider< Topology, CoordTraits, Caching, codim >
    {};

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct NonHybridGeometryProvider< Topology, CoordTraits, Caching, codim, false >
    {
      typedef GenericGeometry :: Geometry
      < SubTopology< Topology, codim, 0 >, Traits, Caching >
      Geometry;

      static Geometry *
      geometry ( unsigned int i,
                 const typename Geometry :: CoordVector &coords,
                 const typename Geometry :: CachingType &cache )
      {
        assert( i < numSubGeometries );
        return new Geometry( coords, cache );
      }
    };



    template< class ElementGeometry, unsigned int codim >
    struct GeometryProvider;

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct GeometryProvider< Geometry< Topology, CoordTraits, Caching >, codim >
      : public NonHybridGeometryProvider
        < Topology, CoordTraits, Caching, codim, IsCodimHybrid< Topology, codim > :: value >
    {};

    template< class Topology, class CoordTraits, template< class > class Caching,
        unsigned int codim >
    struct GeometryProvider< HybridGeometry< Topology, CoordTraits, Caching >, codim >
      : public HybridGeometryProvider< Topology, CoordTraits, Caching, codim >
    {};

  }

}

#endif
