// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CACHEDCOORDFUNCTION_HH
#define DUNE_GEOGRID_CACHEDCOORDFUNCTION_HH

#include <cassert>
#include <map>
#include <vector>

#include <dune/common/interfaces.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/geogrid/coordfunction.hh>

//#define ALLCODIM_SUBINDEX

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class CachedCoordFunction;

  namespace GeoGrid
  {

    template< class HostGrid, class Coordindate, bool hasHierarchicIndexSet >
    class CoordCache;

  }



  // GeoGrid :: CoordCache
  // ---------------------

  namespace GeoGrid
  {

    template< class HostGrid, class Coordinate >
    class CoordCache< HostGrid, Coordinate, true >
    {
      typedef CoordCache< HostGrid, Coordinate, true > This;

      static const unsigned int dimension = HostGrid :: dimension;

      typedef typename HostGrid :: template Codim< dimension > :: Entity Vertex;

      typedef typename HostGrid :: HierarchicIndexSet HierarchicIndexSet;

      typedef std :: vector< Coordinate > DataCache;

      const HierarchicIndexSet &indexSet_;
      DataCache data_;

    public:
      explicit CoordCache ( const HostGrid &hostGrid )
        : indexSet_( hostGrid.hierarchicIndexSet() ),
          data_( indexSet_.size( dimension ) )
      {}

    private:
      CoordCache ( const This & );
      This &operator= ( const This & );

    public:
      template< class Entity >
      const Coordinate &operator() ( const Entity &entity, unsigned int corner ) const
      {
#ifdef ALLCODIM_SUBINDEX
        static const int codim = Entity::codimension;
        static const int subcodim = dimension - codim;
        return data_[ indexSet_.template subIndex< codim, subcodim >( entity, corner ) ];
#else
        return data_[ indexSet_.template subIndex< dimension >( entity, corner ) ];
#endif
      }

#ifndef ALLCODIM_SUBINDEX
      const Coordinate &operator() ( const Vertex &vertex, unsigned int corner ) const
      {
        assert( corner == 0 );
        return data_[ indexSet_.index( vertex ) ];
      }
#endif

      template< class Entity >
      Coordinate &operator() ( const Entity &entity, unsigned int corner )
      {
#ifdef ALLCODIM_SUBINDEX
        static const int codim = Entity :: codimension;
        static const int subcodim = dimension - codim;
        return data_[ indexSet_.template subIndex< codim, subcodim >( entity, corner ) ];
#else
        return data_[ indexSet_.template subIndex< dimension >( entity, corner ) ];
#endif
      }

#ifndef ALLCODIM_SUBINDEX
      Coordinate &operator() ( const Vertex &vertex, unsigned int corner )
      {
        assert( corner == 0 );
        return data_[ indexSet_.index( vertex ) ];
      }
#endif

      void adapt ()
      {
        data_.resize( indexSet_.size( dimension ) );
      }
    };


    template< class HostGrid, class Coordinate >
    class CoordCache< HostGrid, Coordinate, false >
    {
      typedef CoordCache< HostGrid, Coordinate, false > This;

      static const unsigned int dimension = HostGrid :: dimension;

      typedef typename HostGrid :: template Codim< dimension > :: Entity Vertex;

      typedef typename HostGrid :: Traits :: LocalIdSet LocalIdSet;
      typedef typename LocalIdSet :: IdType Id;

      typedef std :: map< Id, Coordinate > DataCache;

      const LocalIdSet &idSet_;
      mutable DataCache data_;

    public:
      explicit CoordCache ( const HostGrid &hostGrid )
        : idSet_( hostGrid.localIdSet() )
      {}

    private:
      CoordCache ( const This & );
      This &operator= ( const This & );

    public:
      template< class Entity >
      const Coordinate &operator() ( const Entity &entity, unsigned int corner ) const
      {
        const Id id = idSet_.template subId< dimension >( entity, corner );
        return data_[ id ];
      }

      const Coordinate &operator() ( const Vertex &vertex, unsigned int corner ) const
      {
        assert( corner == 0 );
        const Id id = idSet_.id( vertex );
        return data_[ id ];
      }

      template< class Entity >
      Coordinate &operator() ( const Entity &entity, unsigned int corner )
      {
        const Id id = idSet_.template subId< dimension >( entity, corner );
        return data_[ id ];
      }

      Coordinate &operator() ( const Vertex &vertex, unsigned int corner )
      {
        assert( corner == 0 );
        const Id id = idSet_.id( vertex );
        return data_[ id ];
      }

      void adapt ()
      {}
    };

  }



  // CachedCoordFunction
  // -------------------

  template< class HostGrid, class CoordFunction >
  class CachedCoordFunction
    : public DiscreteCoordFunction
      < typename CoordFunction::ctype, CoordFunction::dimRange,
          CachedCoordFunction< HostGrid, CoordFunction > >
  {
    typedef CachedCoordFunction< HostGrid, CoordFunction > This;
    typedef DiscreteCoordFunction
    < typename CoordFunction :: ctype, CoordFunction :: dimRange, This >
    Base;

  public:
    typedef typename Base :: RangeVector RangeVector;

  private:
    static const bool hasHierarchicIndexSet
      = Conversion< HostGrid, HasHierarchicIndexSet > :: exists;
    typedef GeoGrid::CoordCache< HostGrid, RangeVector, hasHierarchicIndexSet > Cache;

    const HostGrid &hostGrid_;
    const CoordFunction &coordFunction_;
    Cache cache_;

  public:
    explicit
    CachedCoordFunction ( const HostGrid &hostGrid,
                          const CoordFunction &coordFunction = CoordFunction() )
      : hostGrid_( hostGrid ),
        coordFunction_( coordFunction ),
        cache_( hostGrid )
    {
      buildCache();
    }

    void adapt ()
    {
      cache_.adapt();
      buildCache();
    }

    inline void buildCache ();

    template< class HostEntity >
    inline void insertEntity ( const HostEntity &hostEntity );

    template< class HostEntity >
    void evaluate ( const HostEntity &hostEntity, unsigned int corner,
                    RangeVector &y ) const
    {
      y = cache_( hostEntity, corner );
#ifndef NDEBUG
      RangeVector z;
      calculate( hostEntity.geometry(), corner, z );
      assert( ((y - z).two_norm() < 1e-6) );
#endif
    }

    template< class HostGeometry >
    void calculate ( const HostGeometry &hostGeometry, unsigned int corner,
                     RangeVector &y ) const
    {
      const int dimension = HostGeometry::dimension;
      typedef GenericGeometry::MapNumberingProvider< dimension > Map;

      const int tid = GenericGeometry::topologyId( hostGeometry.type() );
      const int i = Map::template dune2generic< dimension >( tid, corner );
      coordFunction_.evaluate( hostGeometry.corner( i ), y );
    }
  };



  template< class HostGrid, class CoordFunction >
  inline void CachedCoordFunction< HostGrid, CoordFunction > :: buildCache ()
  {
    typedef typename HostGrid :: template Codim< 0 > :: Entity Element;
    typedef typename HostGrid :: template Codim< 0 >
    :: template Partition< All_Partition > :: LevelIterator
    LevelIterator;
    typedef typename HostGrid :: Traits :: HierarchicIterator HierarchicIterator;

    const int maxLevel = hostGrid_.maxLevel();
    const LevelIterator macroEnd = hostGrid_.template lend< 0, All_Partition >( 0 );
    for( LevelIterator macroIt = hostGrid_.template lbegin< 0, All_Partition >( 0 );
         macroIt != macroEnd; ++macroIt )
    {
      const Element &macroElement = *macroIt;
      insertEntity( macroElement );

      const HierarchicIterator hEnd = macroElement.hend( maxLevel );
      for( HierarchicIterator hIt = macroElement.hbegin( maxLevel ); hIt != hEnd; ++hIt )
        insertEntity( *hIt );
    }
  }


  template< class HostGrid, class CoordFunction >
  template< class HostEntity >
  inline void CachedCoordFunction< HostGrid, CoordFunction >
  :: insertEntity ( const HostEntity &hostEntity )
  {
    typedef typename HostEntity::Geometry HostGeometry;

    const HostGeometry &hostGeo = hostEntity.geometry();
    const unsigned int numCorners = hostGeo.corners();
    for( unsigned int i = 0; i < numCorners; ++i )
    {
      RangeVector &y = cache_( hostEntity, i );
      calculate( hostGeo, i, y );
    }
  }

}

#endif
