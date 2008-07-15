// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBGEOMETRY_HH
#define DUNE_GENERICGEOMETRY_SUBGEOMETRY_HH

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, class Traits, unsigned int codim >
    class HybridSubGeometryProvider
    {
      struct CreatorInterface;

      template< int i >
      struct CreatorImplementation;

      typedef const CreatorInterface *CreatorPtr;

      enum { numSubGeometries = Size< Topology, codim > :: value };

      CreatorPtr creator_[ numSubGeometries ];

      enum { dimension = Topology :: dimension - codim };

    public:
      typedef HybridGeometry< dimension, Traits > SubGeometry;

      static SubGeometry *
      subGeometry ( unsigned int i, const typename SubGeometry :: CoordVector &coords )
      {
        assert( i < numSubGeometries );
        return instance().creator_[ i ]->create( coords );
      }

    private:
      HybridSubGeometryProvider ()
      {
        ForLoop< CreatorImplementation, 0, numSubGeometries-1 > :: apply( creator_ );
      }

      static const HybridSubGeometryProvider &instance ()
      {
        static HybridSubGeometryProvider inst;
        return inst;
      }
    };


    template< class Topology, class Traits, unsigned int codim >
    struct HybridSubGeometryProvider< Topology, Traits, codim >
    :: CreatorInterface
    {
      virtual SubGeometry *
      create ( const typename SubGeometry :: CoordVector &coords ) const = 0;
    };

    template< class Topology, class Traits, unsigned int codim >
    template< int i >
    struct HybridSubGeometryProvider< Topology, Traits, codim >
    :: CreatorImplementation
      : public CreatorInterface
    {
      typedef typename GenericGeometry
      :: SubTopology< Topoogy, codim, (unsigned int) i > :: Type
      SubTopology;

      virtual SubGeometry *
      create ( const typename SubGeometry :: CoordVector &coords ) const
      {
        return new VirtualGeometry< SubTopology, Traits >( coords );
      }

      static void apply ( CreatorPtr (&creator)[ numSubGeometries ] )
      {
        creator[ i ] = new CreatorImplementation;
      }
    };



    template< class Topology, class Traits, unsigned int codim, bool isHybrid >
    struct SubGeometryProviderBase;

    template< class Topology, class Traits, unsigned int codim >
    struct SubGeometryProviderBase< Topology, Traits, codim, true >
      : public HybridSubGeometryProvider< Topology, Traits, codim >
    {};

    template< class Topology, class Traits, unsigned int codim >
    struct SubGeometryProviderBase< Topology, Traits, codim, false >
    {
      typedef Geometry< SubTopology< Topology, codim, 0 >, Traits > SubGeometry;

      static SubGeometry *
      subGeometry ( unsigned int i, const typename SubGeometry :: CoordVector &coords )
      {
        assert( i < numSubGeometries );
        return new SubGeometry( coords );
      }
    };

    template< class Topology, class Traits, unsigned int codim >
    struct SubGeometryProvider
      : public SubGeometryProviderBase
        < Topology, Traits, codim, ((codim != 0) && IsHybrid< Topology > :: value) >
    {};

  }

}

#endif
