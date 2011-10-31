// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DEFAULTINDEXSETS_HH
#define DUNE_DEFAULTINDEXSETS_HH

//- system includes
#include <vector>
#include <rpc/rpc.h>

//- Dune includes
#include <dune/common/misc.hh>
// #include <dune/common/interfaces.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh> // for compatibility only
#include <dune/grid/utility/persistentcontainer.hh>

/** @file
   @author Robert Kloefkorn
   @brief Provides default index set implementations for Level- and
   LeafIndexsets used by ALUGrid.
 */

namespace Dune {

  //! LevelIterator tpyes for all codims and partition types
  template <class GridImp>
  struct DefaultLevelIteratorTypes
  {
    //! The types
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator Iterator;
      };
    };
  };

  //! LeafIterator tpyes for all codims and partition types
  template <class GridImp>
  struct DefaultLeafIteratorTypes
  {
    //! The types of the iterator
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::Traits::template Codim<cd>::
        template Partition<pitype>::LeafIterator Iterator;
      };
    };
  };


  namespace DefaultIndexSetHelper
  {

    template< class Grid, class Index >
    struct ContainsIndex
    {
      static bool
      contains ( const PersistentContainer< Grid, Index > &container,
                 const size_t index )
      {
        DUNE_THROW( NotImplemented, "DefaultIndexSetHelper::ContainsIndex not implemented for this grid." );
        return true;
      }
    };

  } // namespace DefaultIndexSetHelper


  /*! \brief
     DefaultIndexSet creates an index set by using the grids persistent container
     an a given pair of iterators
   */
  template < class GridImp, class IteratorImp >
  class DefaultIndexSet :
    public IndexSet< GridImp, DefaultIndexSet <GridImp, IteratorImp>, unsigned int >
  {
    typedef GridImp GridType;
    enum { dim = GridType :: dimension };

  public:
    enum { ncodim = GridType::dimension + 1 };

    //! type of index
    typedef unsigned int IndexType;

    //! type of iterator to generate index set
    typedef IteratorImp IteratorType ;

    struct Index
    {
      int index_;
      Index() : index_( -1 ) {}
      int index() const { return index_; }
      void set( const int index ) { index_ = index; }
    };

    typedef PersistentContainer< GridType, Index > PersistentContainerType ;
    typedef std::vector< PersistentContainerType* > PersistentContainerVectorType;

  private:
    typedef DefaultIndexSet<GridType, IteratorType > ThisType;

    template <class EntityType, int codim>
    struct InsertEntity
    {
      template <class SizeVector>
      static void insert(const EntityType & en,
                         PersistentContainerVectorType &indexContainer,
                         SizeVector& sizes)
      {
        PersistentContainerType& codimContainer = *(indexContainer[ codim ]);
        for( int i = 0; i < en.template count< codim >(); ++i )
        {
          Index& idx = codimContainer( en , i );
          if( idx.index() < 0 )
          {
            idx.set( sizes[codim] );
            ++ sizes[ codim ];
          }
        }
        InsertEntity<EntityType,codim-1>::insert(en, indexContainer, sizes);
      }
    };

    template <class EntityType>
    struct InsertEntity<EntityType,0>
    {
      template <class SizeVector>
      static void insert(const EntityType & en,
                         PersistentContainerVectorType &indexContainer,
                         SizeVector& sizes)
      {
        enum { codim = 0 };
        PersistentContainerType& codimContainer = *(indexContainer[ codim ]);
        Index& idx = codimContainer[ en ];
        if( idx.index() < 0 )
        {
          idx.set( sizes[codim] );
          ++ sizes[ codim ];
        }
      }
    };

  public:
    //! import default implementation of subIndex<cc>
    //! \todo remove after next release
    using IndexSet<GridType, DefaultIndexSet>::subIndex;

    //! create index set by using the given begin and end iterator
    //! for the given level (level == -1 means leaf level)
    DefaultIndexSet( const GridType & grid ,
                     const IteratorType& begin,
                     const IteratorType& end,
                     const int level = -1 )
      : grid_(grid),
        indexContainers_( ncodim, (PersistentContainerType *) 0),
        size_( ncodim, -1 ),
        level_(level)
    {
      for( int codim=0; codim < ncodim; ++codim )
        indexContainers_[ codim ] = new PersistentContainerType( grid, codim );

      calcNewIndex (begin, end);
    }

    //! desctructor deleting persistent containers
    ~DefaultIndexSet ()
    {
      for( int codim=0; codim < ncodim; ++codim )
        delete indexContainers_[ codim ];
    }

    const PersistentContainerType& indexContainer( const size_t codim ) const
    {
      assert( codim < indexContainers_.size() );
      assert( indexContainers_[ codim ] );
      return *( indexContainers_[ codim ] );
    }

    PersistentContainerType& indexContainer( const size_t codim )
    {
      assert( codim < indexContainers_.size() );
      assert( indexContainers_[ codim ] );
      return *( indexContainers_[ codim ] );
    }

    //! return LevelIndex of given entity
    template<class EntityType>
    IndexType index (const EntityType & en) const
    {
      enum { cd = EntityType :: codimension };
      // this must not be true for vertices
      // therefore only check other codims
#ifndef NDEBUG
      const int codim = cd;
      assert( (codim == dim) ? (1) : ( level_ < 0 ) || (level_ == en.level() ));
      assert( indexContainer( codim )[ en ].index() >= 0 );
#endif
      return indexContainer( cd )[ en ].index();
    }

    //! return LevelIndex of given entity
    template<int cd>
    IndexType index (const typename GridImp::template Codim<cd>::Entity& en) const
    {
      // this must not be true for vertices
      // therefore only check other codims
#ifndef NDEBUG
      const int codim = cd;
      //const bool isLeaf = (codim == 0) ? en.isLeaf() : true ;
      assert( (codim == dim) ? (1) : ( level_ < 0 ) || (level_ == en.level() ));
      assert( indexContainer( cd )[ en ].index() >= 0 );
#endif
      return indexContainer( cd )[ en ].index();
    }

    //! return subIndex (LevelIndex) for a given Entity of codim = 0 and a
    //! given SubEntity codim and number of SubEntity
    template< int cc >
    IndexType subIndex ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      assert( (codim != 0) || (level_ < 0) || ( level_ == e.level() ) );
      assert( indexContainer( codim ) ( e, i ).index() >= 0 );
      return indexContainer( codim ) ( e, i ).index();
    }

    //! returns true if this set provides an index for given entity
    template<class EntityType>
    bool contains (const EntityType& en) const
    {
      enum { cd = EntityType :: codimension };
      return (indexContainer( cd )[ en ].index() >= 0 );
    }

    //! return size of IndexSet for a given level and codim
    IndexType size ( int codim ) const
    {
      assert( codim >= 0 && codim <= GridType::dimension );
      return size_[ codim ];
    }

    //! return size of IndexSet for a given level and codim
    //! this method is to be revised
    IndexType size ( GeometryType type ) const
    {
      if( typeNotValid(type) ) return 0;
      return size_[GridType::dimension-type.dim()];
    }

    //! do calculation of the index set, has to be called when grid was
    //! changed or if index set is created
    void calcNewIndex (const IteratorType& begin, const IteratorType& end )
    {
      // resize arrays to new size
      // and set size to zero
      for(int cd=0; cd<ncodim; ++cd)
      {
        indexContainer( cd ).clear();
        size_[cd] = 0;
      }

      // grid walk to setup index set
      for(IteratorType it = begin; it != end; ++it)
      {
        assert( ( level_ < 0 ) ? it->isLeaf() : (it->level() == level_) );
        insertEntity( *it, size_ );
      }

      // remember the number of entity on level and cd = 0
      for(int cd=0; cd<ncodim; ++cd)
      {
        //if( size_[cd] != grid_.size( level_, cd) )
        //{
        //  std::cout << size_[cd] << " s | g " << grid_.size( level_, cd) << std::endl;
        //}
#ifndef NDEBUG
        const int gridSize = ( level_ < 0 ) ? grid_.size( cd ) : grid_.size( level_, cd);
        const int mySize = size_[cd];
        assert( mySize <= gridSize );
#endif
      }
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_.geomTypes( codim );
    }

    //! returns true if this set provides an index for given entity
    bool containsIndex ( const int cd , const int idx) const
    {
      typedef DefaultIndexSetHelper::ContainsIndex< GridType, Index > ContainsIndex;
      return ContainsIndex::contains( indexContainer( cd ), idx );
    }

  private:
    // return whether set has this type stored or not
    bool typeNotValid (const GeometryType & type) const
    {
      int codim = GridType :: dimension - type.dim();
      const std::vector<GeometryType> & geomT = geomTypes(codim);
      for(size_t i=0; i<geomT.size(); ++i) if(geomT[i] == type) return false;
      return true;
    }

    // calculate index for the codim
    template <class EntityType, class SizeVector>
    void insertEntity(EntityType & en, SizeVector& sizes)
    {
      InsertEntity<EntityType,dim>::insert( en, indexContainers_, sizes);
    }

    // grid this index set belongs to
    const GridType& grid_;

    //! vector with PersistentContainer for each codim
    PersistentContainerVectorType indexContainers_;

    // number of entitys of each level an codim
    std::vector< int > size_;

    // the level for which this index set is created
    const int level_;

  };


} // end namespace Dune
#endif
