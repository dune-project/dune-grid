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
  class SingleTypeSizeCache
  {
    typedef SingleTypeSizeCache<GridImp> ThisType;
    //! our dimension
    enum { dim    = GridImp::dimension   };

    //! number of codims
    enum { nCodim = GridImp::dimension+1 };

    typedef GridImp GridType;

    // stores all sizes of the levels
    mutable std::vector<int> levelSizes_[nCodim];

    // stores all sizes of leafs
    mutable int leafSizes_[nCodim];

    // stores all sizes of the levels
    mutable std::vector<int> ghostLevelSizes_[nCodim];

    // stores all sizes of leafs
    mutable int ghostLeafSizes_[nCodim];

    // the grid
    const GridType & grid_;

    // true if this class counts simplices
    const bool isSimplex_;
    // true if this class counts cubes
    const bool isCube_;

    bool notWorry_;

    // count elements of set by iterating the grid
    template <class SzCacheType ,PartitionIteratorType pitype, int codim >
    struct CountLevelEntities
    {
      static inline int count (const SzCacheType & sc, int level, int cd)
      {
        if( cd == codim )
        {
          return sc.template countLevelEntities<pitype,codim> (level);
        }
        else
          return CountLevelEntities < SzCacheType, pitype, codim-1> :: count (sc,level,cd);
      }
    };

    // count elements of set by iterating the grid
    template <class SzCacheType,PartitionIteratorType pitype>
    struct CountLevelEntities<SzCacheType,pitype,0>
    {
      static inline int count (const SzCacheType & sc, int level, int cd)
      {
        enum { codim = 0 };
        assert( cd == codim );
        return sc.template countLevelEntities<pitype,codim> (level);
      }
    };


    // count elements of set by iterating the grid
    template <class SzCacheType , PartitionIteratorType pitype, int codim >
    struct CountLeafEntities
    {
      static inline int count (const SzCacheType & sc, int cd)
      {
        if( cd == codim )
        {
          return sc.template countLeafEntities<pitype,codim> ();
        }
        else
          return CountLeafEntities < SzCacheType, pitype, codim-1> :: count (sc,cd);
      }
    };

    // count elements of set by iterating the grid
    template <class SzCacheType,PartitionIteratorType pitype>
    struct CountLeafEntities<SzCacheType,pitype,0>
    {
      static inline int count (const SzCacheType & sc, int cd)
      {
        enum { codim = 0 };
        assert( cd == codim );
        return sc.template countLeafEntities<pitype,codim> ();
      }
    };


    // private copy constructor
    SingleTypeSizeCache (const SingleTypeSizeCache & );
  public:
    SingleTypeSizeCache (const GridType & grid,
                         const bool isSimplex , const bool isCube, bool notWorry = false )
      : grid_(grid) , isSimplex_(isSimplex) , isCube_(isCube), notWorry_ ( notWorry )
    {
      assert( isSimplex_ != isCube_ );
      for(int i=0; i<nCodim; i++)
      {
        leafSizes_[i] = -1;
        ghostLeafSizes_[i] = -1;
      }

      int numMxl = grid_.maxLevel()+1;
      for(int i=0; i<nCodim; i++)
      {
        std::vector<int> & vec = levelSizes_[i];
        vec.resize(numMxl);
        for(int level = 0; level<numMxl; level++) vec[level] = -1;

        std::vector<int> & ghVec = ghostLevelSizes_[i];
        ghVec.resize(numMxl);
        for(int level = 0; level<numMxl; level++) ghVec[level] = -1;
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
        levelSizes_[codim][level] = CountLevelEntities<ThisType,All_Partition,dim>::count(*this,level,codim);
      return levelSizes_[codim][level];
    }

    //! number of entities per level, codim and geometry type in this process
    int size (int level, int codim, GeometryType type) const
    {
      // if isSimplex true, then this is a simplex counting one
      if( (isSimplex_) && (isSimplex_ != type.isSimplex()) ) return 0;
      // if isCube true, then this is a cube counting one
      if( (isCube_)    && (isCube_    != type.isCube()    ) ) return 0;
      return size(level,codim);
    }

    //********************************************************************
    // leaf sizes
    //********************************************************************
    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      assert( codim >= 0 );
      assert( codim < nCodim );
      if( leafSizes_[codim] < 0)
        leafSizes_[codim] = CountLeafEntities<ThisType,All_Partition,dim>::count(*this,codim);
      return leafSizes_[codim];
    };

    //! number of leaf entities per codim and geometry type in this process
    int size (int codim, GeometryType type) const
    {
      // if isSimplex true, then this is a simplex counting one
      if( (isSimplex_) && (isSimplex_ != type.isSimplex()) ) return 0;
      // if isCube true, then this is a cube counting one
      if( (isCube_)    && (isCube_    != type.isCube()   ) ) return 0;
      return size(codim);
    }

    //! return number of ghost entities for codimension codim
    int ghostSize (int codim) const {
      if( ghostLeafSizes_[codim] < 0)
        ghostLeafSizes_[codim] = CountLeafEntities<ThisType,Ghost_Partition,dim>::count(*this,codim);
      return ghostLeafSizes_[codim];
    }

    /** \brief ghostSize is zero for this grid  */
    int ghostSize (int level, int codim) const
    {
      // if level is larger than grid max level , return 0
      if( level >= (int) ghostLevelSizes_[codim].size() ) return 0;

      if( ghostLevelSizes_[codim][level] < 0)
        ghostLevelSizes_[codim][level] = CountLevelEntities<ThisType,Ghost_Partition,dim>::count(*this,level,codim);
      return ghostLevelSizes_[codim][level];
    }

  private:
    template <PartitionIteratorType pitype, int codim>
    int countLevelEntities(int level) const
    {
      typedef typename GridType::template Codim<codim> :: template Partition<pitype> :: LevelIterator LevelIterator;
      LevelIterator it  = grid_.template lbegin<codim,pitype> (level);
      LevelIterator end = grid_.template lend<codim,pitype>   (level);

      GeometryType type (((isSimplex_) ? GeometryType::simplex : GeometryType::cube ),dim-codim);
      assert( ((dim-codim) > 1) ? (type.isCube() == isCube_) : 1);
      if( notWorry_ ) return countElements(it,end,type);
      return countElements(it,end);
    }

    template <PartitionIteratorType pitype, int codim>
    int countLeafEntities() const
    {
      // count All_Partition entities
      typedef typename GridType::template Codim<codim> :: template Partition<pitype> :: LeafIterator LeafIterator;
      LeafIterator it  = grid_.template leafbegin<codim,pitype> ();
      LeafIterator end = grid_.template leafend<codim,pitype>   ();
      GeometryType type (((isSimplex_) ? GeometryType::simplex : GeometryType::cube ),dim-codim);
      assert( ((dim-codim) > 1) ? (type.isCube() == isCube_) : 1);
      if( notWorry_ ) return countElements(it,end,type);
      return countElements(it,end);
    }

    // counts entities with given type for given iterator
    template <class IteratorType>
    int countElements(IteratorType & it, const IteratorType & end ,
                      const GeometryType & type ) const
    {
      int count = 0;
      if((type.isSimplex()) || (type.isCube()))
      {
        for( ; it != end; ++it )
        {
          if(it->type() == type)
            ++ count ;
        }
      }
      return count;
    }

    // counts entities for given iterator
    template <class IteratorType>
    int countElements(IteratorType & it, const IteratorType & end) const
    {
      int count = 0;
      for( ; it != end; ++it )
      {
        ++ count ;
      }
      return count;
    }
  };

  //! oranizes the caching of sizes for one grid
  template <class GridImp>
  class SizeCache
  {
    typedef SizeCache<GridImp> ThisType;
    typedef GridImp GridType;

    SingleTypeSizeCache<GridType> simplexSize_;
    SingleTypeSizeCache<GridType>    cubeSize_;

  public:
    SizeCache (const GridType & grid) : simplexSize_(grid,true,false), cubeSize_(grid,false,true)
    {
      // check that used grid only has simplex and/or cube as geomTypes
      // to be revised
      const std::vector<GeometryType> & geomTypes = grid.geomTypes(0);
      int found  = 0;
      int others = 0;
      for(unsigned int i=0; i<geomTypes.size(); i++)
      {
        if( (geomTypes[i].isSimplex()) ||
            (geomTypes[i].isCube()   )    )
          found++;
        else
          others++;
      }
      // only use with simplex or cube types
      // if others are found assert
      assert( !others );

      // assert that at least one of our tpyes is found
      assert( found > 0 );
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
      return (simplexSize_.size(level,codim) + cubeSize_(level,codim));
    }

    //! number of entities per level, codim and geometry type in this process
    int size (int level, int codim, GeometryType type) const
    {
      if( type.isSimplex()) return simplexSize_.size(level,codim);
      if( type.isCube()   ) return cubeSize_(level,codim);
      return 0;
    }

    //********************************************************************
    // leaf sizes
    //********************************************************************
    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return (simplexSize_.size(codim) + cubeSize_(codim));
    };

    //! number of leaf entities per codim and geometry type in this process
    int size (int codim, GeometryType type) const
    {
      if( type.isSimplex() ) return simplexSize_.size(codim);
      if( type.isCube()    ) return cubeSize_(codim);
      return 0;
    }

    //! number of leaf entities per codim and geometry type in this process
    int ghostSize (int codim) const
    {
      return simplexSize_.ghostSize(codim) + cubeSize_.ghostSize(codim);
    }

    //! number of leaf entities per codim and geometry type in this process
    int ghostSize (int level, int codim) const
    {
      return simplexSize_.ghostSize(level,codim) + cubeSize_.ghostSize(level,codim);
    }
  };


} // end namespace Dune
#endif
