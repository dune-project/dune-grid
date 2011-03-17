// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SIZECACHE_HH
#define DUNE_SIZECACHE_HH

#include <cassert>
#include <vector>

#include <dune/common/geometrytype.hh>

#include <dune/grid/common/gridenums.hh>

/** @file
   @author Robert Kloefkorn
   @brief Provides size cache classes to
   implement the grids size method efficiently.
 */

namespace Dune {

  //! organizes the caching of sizes for one grid and one GeometryType
  template <class GridImp>
  class SizeCache
  {
    typedef SizeCache<GridImp> ThisType;
    //! our dimension
    enum { dim    = GridImp::dimension   };

    //! number of codims
    enum { nCodim = GridImp::dimension+1 };

    typedef GridImp GridType;

    // stores all sizes of the levels
    mutable std::vector< int > levelSizes_[nCodim];

    // stores all sizes of the levels
    mutable std::vector< std::vector< int > > levelTypeSizes_[nCodim];

    // stores all sizes of leafs
    mutable int leafSizes_[nCodim];

    // stores all sizes of leafs
    mutable std::vector< int > leafTypeSizes_[nCodim];

    // the grid
    const GridType & grid_;

    // count elements of set by iterating the grid
    template <class SzCacheType ,PartitionIteratorType pitype, int codim >
    struct CountLevelEntities
    {
      static inline void count (const SzCacheType & sc, int level, int cd)
      {
        if( cd == codim )
        {
          sc.template countLevelEntities<pitype,codim> (level);
        }
        else
          CountLevelEntities < SzCacheType, pitype, codim-1> :: count (sc,level,cd);
      }
    };

    // count elements of set by iterating the grid
    template <class SzCacheType,PartitionIteratorType pitype>
    struct CountLevelEntities<SzCacheType,pitype,0>
    {
      static inline void count (const SzCacheType & sc, int level, int cd)
      {
        enum { codim = 0 };
        assert( cd == codim );
        sc.template countLevelEntities<pitype,codim> (level);
      }
    };


    // count elements of set by iterating the grid
    template <class SzCacheType , PartitionIteratorType pitype, int codim >
    struct CountLeafEntities
    {
      static inline void count (const SzCacheType & sc, int cd)
      {
        if( cd == codim )
        {
          sc.template countLeafEntities<pitype,codim> ();
        }
        else
          CountLeafEntities < SzCacheType, pitype, codim-1> :: count (sc,cd);
      }
    };

    // count elements of set by iterating the grid
    template <class SzCacheType,PartitionIteratorType pitype>
    struct CountLeafEntities<SzCacheType,pitype,0>
    {
      static inline void count (const SzCacheType & sc, int cd)
      {
        enum { codim = 0 };
        assert( cd == codim );
        sc.template countLeafEntities<pitype,codim> ();
      }
    };

    const int gtIndex( const GeometryType& type ) const
    {
      return type.id() >> 1 ;
    }

    const int sizeCodim( const int codim ) const
    {
      const int mydim = GridType :: dimension - codim;
      return ((1 << mydim) + 1) / 2;
    }

    // private copy constructor
    SizeCache (const SizeCache & );
  public:
    SizeCache (const GridType & grid) : grid_( grid )
    {
      for(int i=0; i<nCodim; i++)
      {
        leafSizes_[ i ] = -1;
        leafTypeSizes_[ i ].resize( sizeCodim( i ), -1 );
      }

      int numMxl = grid_.maxLevel()+1;
      for(int i=0; i<nCodim; i++)
      {
        std::vector<int> & vec = levelSizes_[i];
        vec.resize(numMxl);
        levelTypeSizes_[i].resize( numMxl );
        for(int level = 0; level<numMxl; level++)
        {
          vec[level] = -1;
          levelTypeSizes_[i][level].resize( sizeCodim( i ), -1 );
        }
      }
    }

    //********************************************************************
    // level sizes
    //********************************************************************
    /** \brief Number of grid entities per level and codim
     * because lbegin and lend are none const, and we need this methods
     * counting the entities on each level, you know.
     */
    int size (int level, int codim) const
    {
      assert( codim >= 0 );
      assert( codim < nCodim );
      assert( level >= 0 );
      if( level >= (int) levelSizes_[codim].size() ) return 0;

      if( levelSizes_[codim][level] < 0)
        CountLevelEntities<ThisType,All_Partition,dim>::count(*this,level,codim);
      assert( levelSizes_[codim][level] >= 0 );
      return levelSizes_[codim][level];
    }

    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const
    {
      const int codim = GridType ::dimension - type.dim();
      if( levelSizes_[codim][level] < 0)
        CountLevelEntities<ThisType,All_Partition,dim>::count(*this,level,codim);

      assert( levelTypeSizes_[codim][level][gtIndex( type )] >= 0 );
      return levelTypeSizes_[codim][level][gtIndex( type )];
    }

    //********************************************************************
    // leaf sizes
    //********************************************************************
    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      assert( codim >= 0 );
      assert( codim < nCodim );
      if( leafSizes_[codim] < 0 )
        CountLeafEntities<ThisType,All_Partition,dim>::count(*this,codim);
      assert( leafSizes_[codim] >= 0 );
      return leafSizes_[codim];
    };

    //! number of leaf entities per codim and geometry type in this process
    int size ( const GeometryType type ) const
    {
      const int codim = GridType :: dimension - type.dim();
      if( leafSizes_[codim] < 0 )
        CountLeafEntities<ThisType,All_Partition,dim>::count(*this,codim);

      assert( leafTypeSizes_[codim][ gtIndex( type )] >= 0 );
      return leafTypeSizes_[codim][ gtIndex( type )];
    }

  private:
    template <PartitionIteratorType pitype, int codim>
    void countLevelEntities(int level) const
    {
      typedef typename GridType::template Codim<codim> :: template Partition<pitype> :: LevelIterator LevelIterator;
      LevelIterator it  = grid_.template lbegin<codim,pitype> (level);
      LevelIterator end = grid_.template lend<codim,pitype>   (level);
      levelSizes_[codim][level] = countElements(it,end, levelTypeSizes_[codim][level]);
    }

    template <PartitionIteratorType pitype, int codim>
    void countLeafEntities() const
    {
      // count All_Partition entities
      typedef typename GridType::template Codim<codim> :: template Partition<pitype> :: LeafIterator LeafIterator;
      LeafIterator it  = grid_.template leafbegin<codim,pitype> ();
      LeafIterator end = grid_.template leafend<codim,pitype>   ();
      leafSizes_[codim] = countElements(it,end, leafTypeSizes_[codim] );
    }

    // counts entities with given type for given iterator
    template <class IteratorType>
    int countElements(IteratorType & it, const IteratorType & end, std::vector<int>& typeSizes) const
    {
      int overall = 0;
      const size_t types = typeSizes.size();
      for(size_t i=0; i<types; ++i) typeSizes[i] = 0;
      for( ; it != end; ++it )
      {
        const GeometryType type = it->type();
        ++typeSizes[ gtIndex( type ) ];
        ++overall;
      }

      int sumtypes = 0;
      for(size_t i=0; i<types; ++i) sumtypes += typeSizes[i];

      assert( overall == sumtypes );
      return overall;
    }
  };

  //! organizes the caching of sizes for one grid and one GeometryType
  template <class GridImp>
  struct SingleTypeSizeCache : public SizeCache< GridImp >
  {
    DUNE_DEPRECATED SingleTypeSizeCache( const GridImp& grid, bool isSimplex , bool isCube, bool notWorry = false ) : SizeCache< GridImp >( grid ) {}
  };

} // end namespace Dune
#endif
