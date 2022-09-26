// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_SIZECACHE_HH
#define DUNE_GRID_COMMON_SIZECACHE_HH

#include <cassert>
#include <vector>
#include <set>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/capabilities.hh>

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
    constexpr static int dim = GridImp::dimension;

    //! number of codims
    constexpr static int nCodim = GridImp::dimension + 1;

    // type of grid
    typedef GridImp GridType;

    // coordinate type
    typedef typename GridType :: ctype ctype ;

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
    template < int codim, bool gridHasCodim >
    struct CountLevelEntitiesBase
    {
      template < class SzCacheType >
      static void apply(const SzCacheType & sc, int level, int cd)
      {
        if( cd == codim )
        {
          sc.template countLevelEntities<All_Partition,codim> (level);
        }
      }
    };

    template < int codim >
    struct CountLevelEntitiesBase< codim, false >
    {
      template < class SzCacheType >
      static void apply(const SzCacheType & sc, int level, int cd)
      {
        if( cd == codim )
        {
          sc.template countLevelEntitiesNoCodim<All_Partition,codim> (level);
        }
      }
    };

    template < int codim >
    struct CountLevelEntities
      : public CountLevelEntitiesBase< codim, Capabilities :: hasEntity< GridType, codim > :: v >
    {};

    // count elements of set by iterating the grid
    template < int codim, bool gridHasCodim >
    struct CountLeafEntitiesBase
    {
      template <class SzCacheType>
      static void apply(const SzCacheType & sc, int cd)
      {
        if( cd == codim )
        {
          sc.template countLeafEntities<All_Partition,codim> ();
        }
      }
    };

    // count elements of set by iterating the grid
    template < int codim >
    struct CountLeafEntitiesBase< codim, false >
    {
      template <class SzCacheType>
      static void apply(const SzCacheType & sc, int cd)
      {
        if( cd == codim )
        {
          sc.template countLeafEntitiesNoCodim<All_Partition,codim> ();
        }
      }
    };

    template < int codim >
    struct CountLeafEntities
      : public CountLeafEntitiesBase< codim, Capabilities :: hasEntity< GridType, codim > :: v >
    {};

    int gtIndex( const GeometryType& type ) const
    {
      return type.id() >> 1 ;
    }

    int sizeCodim( const int codim ) const
    {
      const int mydim = GridType :: dimension - codim;
      return ((1 << mydim) + 1) / 2;
    }

    // private copy constructor
    SizeCache (const SizeCache & );
  public:
    /** \brief constructor taking grid reference */
    SizeCache (const GridType & grid) : grid_( grid )
    {
      reset();
    }

    /** \brief reset all cached sizes */
    void reset()
    {
      for(int codim=0; codim<nCodim; ++codim)
      {
        leafSizes_[ codim ] = -1;
        leafTypeSizes_[ codim ].resize( sizeCodim( codim ), -1 );
      }

      const int numMxl = grid_.maxLevel()+1;
      for(int codim=0; codim<nCodim; ++codim)
      {
        std::vector<int> & vec = levelSizes_[codim];
        vec.resize(numMxl);
        levelTypeSizes_[codim].resize( numMxl );
        for(int level = 0; level<numMxl; ++level)
        {
          vec[level] = -1;
          levelTypeSizes_[codim][level].resize( sizeCodim( codim ), -1 );
        }
      }
    }

    //********************************************************************
    // level sizes
    //********************************************************************
    /** \copydoc Dune::Grid::size(int level,int codim) const */
    int size (int level, int codim) const
    {
      assert( codim >= 0 );
      assert( codim < nCodim );
      assert( level >= 0 );
      if( level >= (int) levelSizes_[codim].size() ) return 0;

      if( levelSizes_[codim][level] < 0)
        Hybrid::forEach( std::make_index_sequence< dim+1 >{}, [ & ]( auto i ){ CountLevelEntities< i >::apply( *this, level, codim ); } );

      //  CountLevelEntities<ThisType,All_Partition,dim>::count(*this,level,codim);

      assert( levelSizes_[codim][level] >= 0 );
      return levelSizes_[codim][level];
    }

    /** \copydoc Dune::Grid::size(int level,GeometryType type) const */
    int size (int level, GeometryType type) const
    {
      int codim = GridType ::dimension - type.dim();
      if( levelSizes_[codim][level] < 0)
        Hybrid::forEach( std::make_index_sequence< dim+1 >{}, [ & ]( auto i ){ CountLevelEntities< i >::apply( *this, level, codim ); } );

      assert( levelTypeSizes_[codim][level][gtIndex( type )] >= 0 );
      return levelTypeSizes_[codim][level][gtIndex( type )];
    }

    //********************************************************************
    // leaf sizes
    //********************************************************************
    /** \copydoc Dune::Grid::size(int codim) const */
    int size (int codim) const
    {
      assert( codim >= 0 );
      assert( codim < nCodim );
      if( leafSizes_[codim] < 0 )
        Hybrid::forEach( std::make_index_sequence< dim+1 >{}, [ & ]( auto i ){ CountLeafEntities< i >::apply( *this, codim ); } );

      assert( leafSizes_[codim] >= 0 );
      return leafSizes_[codim];
    };

    /** \copydoc Dune::Grid::size(GeometryType type) const */
    int size ( const GeometryType type ) const
    {
      int codim = GridType :: dimension - type.dim();
      if( leafSizes_[codim] < 0 )
        Hybrid::forEach( std::make_index_sequence< dim+1 >{}, [ & ]( auto i ){ CountLeafEntities< i >::apply( *this, codim ); } );

      assert( leafTypeSizes_[codim][ gtIndex( type )] >= 0 );
      return leafTypeSizes_[codim][ gtIndex( type )];
    }

  private:
    template <PartitionIteratorType pitype, int codim>
    void countLevelEntities(int level) const
    {
      typedef typename GridType :: LevelGridView GridView ;
      typedef typename GridView :: template Codim< codim > :: template Partition<pitype>  :: Iterator Iterator ;
      GridView gridView = grid_.levelGridView( level );
      Iterator it  = gridView.template begin<codim,pitype> ();
      Iterator end = gridView.template end<codim,pitype>   ();
      levelSizes_[codim][level] = countElements(it,end, levelTypeSizes_[codim][level]);
    }

    template <PartitionIteratorType pitype, int codim>
    void countLeafEntities() const
    {
      // count All_Partition entities
      typedef typename GridType :: LeafGridView GridView ;
      typedef typename GridView :: template Codim< codim > :: template Partition<pitype>  :: Iterator Iterator ;
      GridView gridView = grid_.leafGridView();
      Iterator it  = gridView.template begin<codim,pitype> ();
      Iterator end = gridView.template end<codim,pitype>   ();
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

#ifndef NDEBUG
      int sumtypes = 0;
      for(size_t i=0; i<types; ++i)
        sumtypes += typeSizes[i];

      assert( overall == sumtypes );
#endif

      return overall;
    }

    template <PartitionIteratorType pitype, int codim>
    void countLevelEntitiesNoCodim(int level) const
    {
      typedef typename GridType :: LevelGridView GridView ;
      typedef typename GridView :: template Codim< 0 > :: template Partition<pitype>  :: Iterator Iterator ;
      GridView gridView = grid_.levelGridView( level );
      Iterator it  = gridView.template begin< 0, pitype> ();
      Iterator end = gridView.template end< 0, pitype>   ();
      levelSizes_[codim][level] = countElementsNoCodim< codim >(it,end, levelTypeSizes_[codim][level]);
    }

    template <PartitionIteratorType pitype, int codim>
    void countLeafEntitiesNoCodim() const
    {
      // count All_Partition entities
      typedef typename GridType :: LeafGridView GridView ;
      typedef typename GridView :: template Codim< 0 > :: template Partition<pitype>  :: Iterator Iterator ;
      GridView gridView = grid_.leafGridView();
      Iterator it  = gridView.template begin< 0, pitype > ();
      Iterator end = gridView.template end< 0, pitype >   ();
      leafSizes_[codim] = countElementsNoCodim< codim >(it,end, leafTypeSizes_[codim] );
    }

    // counts entities with given type for given iterator
    template < int codim, class IteratorType >
    int countElementsNoCodim(IteratorType & it, const IteratorType & end, std::vector<int>& typeSizes) const
    {
      typedef typename GridType :: LocalIdSet LocalIdSet ;
      typedef typename LocalIdSet :: IdType IdType ;

      typedef ReferenceElements< ctype, dim > ReferenceElementContainerType;
      typedef typename ReferenceElementContainerType::ReferenceElement ReferenceElementType;

      typedef std::set< IdType > CodimIdSetType ;

      typedef typename IteratorType :: Entity ElementType ;

      // get id set
      const LocalIdSet& idSet = grid_.localIdSet();

      const size_t types = typeSizes.size();
      for(size_t i=0; i<types; ++i) typeSizes[ i ] = 0;

      std::vector< CodimIdSetType > typeCount( types );

      // count all elements of codimension codim
      for( ; it != end; ++it )
      {
        // get entity
        const ElementType& element = *it ;
        // get reference element
        ReferenceElementType refElem =
          ReferenceElementContainerType :: general( element.type() );

        // count all sub entities of codimension codim
        const int count = element.subEntities( codim );
        for( int i=0; i< count; ++ i )
        {
          // get geometry type
          const GeometryType geomType = refElem.type( i, codim );
          // get id of sub entity
          const IdType id = idSet.subId( element, i, codim );
          // insert id into set
          typeCount[ gtIndex( geomType ) ].insert( id );
        }
      }

      // accumulate numbers
      int overall = 0;
      for(size_t i=0; i<types; ++i)
      {
        typeSizes[ i ] = typeCount[ i ].size();
        overall += typeSizes[ i ];
      }

      return overall;
    }
  };

} // end namespace Dune
#endif
