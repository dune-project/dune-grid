// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <cassert>
#include <map>
#include <vector>

#include <dune/common/misc.hh>
#include <dune/common/forloop.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune
{
  /** \brief A class for storing data during an adaptation cycle.
   *
   * This container allows to store data which is to remain persistent
   * even during adaptation cycles. There is a default implementation based
   * on std::map but any grid implementation can provide a specialized implementation.
   *
   * The container expects that the class Data has a default constructor. This default
   * value is returned if data for an entity is requested which has not yet been added
   * to the container. Therefore no find as for std::map is required and provided.
   * The class furthermore provides a size method and iterators. The
   * guarantee being that the iteration is at least over all entries in the container
   * which are not set to default but the iteration could also include entries with the default
   * value. The size method returns the number of elements visited by an iteration so
   * that size() is larger or equal to the entries which are not set to the default value.
   * The method clear can be used to set all entries in the container to the default value.
   *
   * After the grid has changed update() or reserve() should to be called.
   * In contrast to the method reserve, the method update can compress the data to
   * reduce memory consumption, whereas reserve() can be implemented not to reduce any
   * memory even if the size of the grid has changed.
   *
   * Specialized implementations of this class are provided for some grid
   * implementations. Especially grids with an id type with a hash function could use
   * the std::unordered_map to store the data (the implementation Dune::PersistentContainerMap
   * can be used). For grids providing an id set suitable for storage in vector like
   * structures (id type is int and size method provided) Dune::PersistentContainerVector can be
   * used.
   */
  template < class Grid, class Data, class Allocator=std::allocator<Data> >
  class PersistentContainer;

  /** \brief An implementation for the PersistentContainer based on a container
   * satisfying the std::vector interface and using a class providing an IndexSet
   * for storing the Data*/
  template< class Grid, class IndexSet, class Vector >
  class PersistentContainerVector
  {
  public:
    typedef typename Vector::value_type Data;
    typedef Grid GridType;

    //! \brief entity of codimension 0
    typedef typename GridType :: template Codim< 0 > :: Entity ElementType;

    //! \brief iterator type
    typedef typename Vector :: iterator Iterator ;
    //! \brief const iterator type
    typedef typename Vector :: const_iterator ConstIterator ;

    //! \brief constructor creating vector container filled with default value
    //         store data on entities of given codim using indexSet to store data in vector.
    //         The overEstimate parameter can be used to allocate more memory than
    //         required to store the data.
    PersistentContainerVector ( const GridType &grid, const int codim,
                                const IndexSet &indexSet,
                                const double overEstimate,
                                const typename Vector::allocator_type &allocator )
      : codim_( codim ),
        indexSet_( indexSet ),
        overEstimate_( overEstimate ),   // this is not yet the right approach - will be revised
        data_( indexSet_.size( codim ), Data(), allocator )
    {}

    //! pass on index set used for Container
    const IndexSet &indexSet ()
    {
      return indexSet_;
    }

    //! \brief random access to entity data with correct codimension
    template <class Entity>
    Data& operator [] (const Entity& entity )
    {
      assert( Entity :: codimension == codim_ );
      assert( (typename IndexSet::IndexType)indexSet_.index( entity ) < (typename IndexSet::IndexType)data_.size() );
      return data_[ indexSet_.index( entity ) ];
    }

    //! \brief random access to entity data with correct codimension
    template <class Entity>
    const Data& operator [] (const Entity& entity ) const
    {
      assert( Entity :: codimension == codim_ );
      assert( (typename IndexSet::IndexType)indexSet_.index( entity ) < (typename IndexSet::IndexType)data_.size() );
      return data_[ indexSet_.index( entity ) ];
    }

    //! \brief access for sub entity data
    template< class Entity >
    Data &operator () ( const Entity &entity, const int subEntity )
    {
      const typename IndexSet::IndexType idx = indexSet_.subIndex( entity, subEntity, codim_ );
      assert( idx < (typename IndexSet::IndexType)data_.size() );
      return data_[ idx ];
    }

    //! \brief access for sub entity data
    template< class Entity >
    const Data &operator () ( const Entity &entity, const int subEntity ) const
    {
      const typename IndexSet::IndexType idx = indexSet_.subIndex( entity, subEntity, codim_ );
      assert( idx < (typename IndexSet::IndexType)data_.size() );
      return data_[ idx ];
    }

    //! \brief const iterator begin
    Iterator begin()
    {
      return data_.begin();
    }

    //! \brief const iterator begin
    ConstIterator begin() const
    {
      return data_.begin();
    }

    //! \brief iterator end
    Iterator end()
    {
      return data_.end();
    }

    //! \brief const iterator end
    ConstIterator end() const
    {
      return data_.end();
    }

    //! \brief return size of allocated data
    size_t size() const { return data_.size(); }

    //! \brief enlarge container, compress is not necessary but could be done
    void reserve( ) // need a better name
    {
      if( (typename IndexSet::IndexType)indexSet_.size( codim_ ) > (typename IndexSet::IndexType)data_.size() )
        update( );
    }

    //! \brief adjust container to correct size and set all values to default
    void clear( )
    {
      // clear all entries
      data_.clear();
      // resize with new default value
      const size_t newSize = indexSet_.size( codim_ );
      data_.resize( newSize, Data() );
    }

    //! \brief adjust container to correct size including compress
    void update( )
    { // this could be more sophisticated (although std::vector is not stupid and
      // overestimated on its own...
      const size_t newSize = indexSet_.size( codim_ );
      if (newSize < data_.capacity())
        data_.resize(newSize, Data() );
      else
      {
        data_.reserve(newSize*overEstimate_);
        data_.resize(newSize, Data() );
      }
    }

  protected:
    const int codim_;
    const IndexSet &indexSet_;
    const double overEstimate_;
    Vector data_;
  };

  /** \brief An implementation for the PersistentContainer based on a container
   * satisfying the std::map interface and using a class providing an IdSet
   * for storing the Data*/
  template< class Grid, class IdSet, class Map >
  class PersistentContainerMap
  {
    typedef PersistentContainerMap< Grid, IdSet, Map > ThisType;

  protected:
    typedef typename Map :: mapped_type Data;
    typedef typename IdSet :: IdType IdType;
    typedef Grid GridType;
    const GridType& grid_;
    const int codim_;
    const IdSet& idSet_;
    mutable Map data_;

    typedef typename Map :: iterator iterator ;
    typedef typename Map :: const_iterator const_iterator ;

    template <class D, class IteratorType >
    struct DataExtractor ;

    // Data type for iterator
    template <class D>
    struct DataExtractor< D, iterator >
    {
      typedef D Type ;
    };

    // Data type for const iterator
    template <class D>
    struct DataExtractor< D, const_iterator >
    {
      typedef const D Type ;
    };

    template <class IteratorType>
    class MyIterator
    {
      IteratorType it_;
    public:
      // get correct data type (const or non-const)
      typedef typename DataExtractor<Data, IteratorType> :: Type value_type ;

      MyIterator(const IteratorType& it) : it_( it ) {}
      MyIterator(const MyIterator& other) : it_( other.it_ ) {}

      bool operator == (const MyIterator& other) const { return it_ == other.it_; }
      bool operator != (const MyIterator& other) const { return it_ != other.it_; }

      MyIterator& operator ++ ()
      {
        ++it_;
        return *this;
      }
      value_type& operator * () { return (*it_).second; }
      value_type* operator -> () { return &((*it_).second); }
      MyIterator& operator = (const MyIterator& other)
      {
        it_ = other.it_;
        return *this;
      }
    };

    template< int codim , bool gridHasCodim >
    struct AdaptCodimBase
    {
      static void apply ( ThisType &container, const Data& value , const int myCodim)
      {
        if( codim == myCodim )
          container.template adaptCodim< codim > ( value );
      }
    };

    template< int codim >
    struct AdaptCodimBase< codim, false >
    {
      static void apply ( ThisType &container, const Data& value , const int myCodim)
      {}
    };

    template< int codim >
    struct AdaptCodim
      : public AdaptCodimBase< codim, Capabilities :: hasEntity < GridType, codim > :: v >
    {};

  public:
    typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
    typedef MyIterator< iterator > Iterator;
    typedef MyIterator< const_iterator > ConstIterator;

    typedef typename Map::key_compare Compare;
    typedef typename Map::allocator_type Allocator;

    //! \brief constructor creating container filled with default values.
    //
    //         Container is to be used to store data on entities of given codim using idSet to store data in map.
    PersistentContainerMap ( const GridType &grid, const int codim, const IdSet &idSet,
                             const Compare &comp = Compare(), const Allocator &allocator = Allocator() )
      : grid_( grid ),
        codim_( codim ),
        idSet_( idSet ),
        data_( comp, allocator )
    {}

    //! \brief random access entity with correct codimension
    template <class Entity>
    Data& operator [] (const Entity& entity )
    {
      assert( Entity :: codimension == codim_ );
      return data_[ idSet_.id( entity ) ];
    }

    //! \brief random access entity with correct codimension
    template <class Entity>
    const Data& operator [] (const Entity& entity ) const
    {
      assert( Entity :: codimension == codim_ );
      return data_[ idSet_.id( entity ) ];
    }

    //! \brief access for sub entity data
    template< class Entity >
    Data &operator () ( const Entity &entity, const int subEntity )
    {
      return data_[ idSet_.subId( entity, subEntity, codim_ ) ];
    }

    //! \brief access for sub entity data
    template< class Entity >
    const Data &operator () ( const Entity &entity, const int subEntity ) const
    {
      return data_[ idSet_.subId( entity, subEntity, codim_ ) ];
    }

    //! \brief iterator begin for iterating over data actually stored in container
    Iterator begin()
    {
      return Iterator( data_.begin() );
    }

    //! \brief const iterator begin
    ConstIterator begin() const
    {
      return ConstIterator( data_.begin() );
    }

    //! \brief iterator end
    Iterator end()
    {
      return Iterator( data_.end() );
    }

    //! \brief const iterator end
    ConstIterator end() const
    {
      return ConstIterator( data_.end() );
    }

    //! \brief return size of allocated data
    size_t size() const { return data_.size(); }

    //! \brief enlarge container, compress is not necessary but could be done
    void reserve()
    {}

    //! \brief adjust container to correct size and set all values to default
    void clear( )
    {
      data_.clear();
    }

    //! \brief adjust container to correct size including compress
    void update( )
    { // this version could be implemented differently by only compressing
      update( Data() );
    }
  protected:
    //! \brief adjust container to correct size including compress
    void update( const Data& value )
    {
      // loop over all codimensions (needed to make codim_ static)
      ForLoop< AdaptCodim, 0, GridType :: dimension > :: apply( *this, value, codim_ );
    }

    template <int codim>
    void adaptCodim( const Data& value )
    {
      assert( codim_ == codim );
      // create empty map and swap it with current map (no need to copy twice)
      Map oldData;
      std::swap( oldData, data_ );

      const iterator olddataend = oldData.end();
      typedef typename GridType :: template Codim< codim > :: LevelIterator LevelIterator ;
      typedef typename LevelIterator :: Entity Entity;
      for(int l = 0; l <= grid_.maxLevel(); ++ l)
      {
        const LevelIterator endit = grid_.template lend< codim > ( l );
        for( LevelIterator it = grid_.template lbegin< codim > ( l ); it != endit; ++ it )
        {
          const Entity& entity = * it ;
          const IdType id = idSet_.id( entity );
          Data& data = data_[ id ];
          iterator entry = oldData.find( id );
          if( entry != olddataend )
            data = (*entry).second;
        }
      }
    }
  };

  // PersistentContainer (default is to use PersistentContainerMap)
  // -------------------
  template < class Grid, class Data, class Allocator>
  class PersistentContainer
    : public PersistentContainerMap< Grid, typename Grid::Traits::LocalIdSet,
          std::map<const typename Grid::Traits::LocalIdSet::IdType, Data,
              std::less<const typename Grid::Traits::LocalIdSet::IdType>,
              typename Allocator::template rebind<typename Grid::Traits::LocalIdSet::IdType>::other > >
  {
  public:
    typedef Grid GridType;
  protected:
    typedef typename Grid::Traits::LocalIdSet IdSet;
    typedef typename IdSet::IdType IdType;
    typedef typename Allocator::template rebind<IdType>::other IdAllocator;
    typedef std::map<const IdType, Data, std::less<const IdType>,
        IdAllocator> Map;
    typedef PersistentContainerMap< Grid, IdSet, Map > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator=Allocator() )
      : BaseType( grid, codim, grid.localIdSet(), std::less<const IdType>(), allocator )
    {}
  };

#if 0 // the following implementation can be used for a grid providing a hash for the id type
#include <unordered_map>
  template < class MyGrid, class Data, class Allocator >
  class PersistentContainer
    : public PersistentContainerMap< MyGrid, typename MyGrid::Traits::LocalIdSet,
          std::unordered_map<const typename MyGrid::Traits::LocalIdSet::IdType, Data,
              std::hash<typename MyGrid::Traits::LocalIdSet::IdType>,
              std::equal_to<const typename MyGrid::Traits::LocalIdSet::IdType>, Allocator> >
  {
    typedef MyGrid GridType;
    typedef typename GridType::Traits::LocalIdSet IdSet;
    typedef typename IdSet::IdType IdType;
    typedef std::unordered_map<const IdType, Data, std::hash<IdType>, std::equal_to<const IdType>, Allocator> Map;
    typedef PersistentContainerMap< GridType, IdSet, Map > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    //! Note that we can not use a non-default allocator for this implementation due to
    //! the missing constructor (PersistentContainerMap has to be reimplemented for the
    //! hash-map.
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator=Allocator() )
      : BaseType( grid, codim, grid.localIdSet() )
    {}
  };
#endif
} // end namespace Dune

#endif // end DUNE_PERSISTENTCONTAINER_HH
