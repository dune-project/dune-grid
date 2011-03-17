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

  /*!
     The DefaultGridIndexSet is a wrapper for the grid index which can be
     index of entity or globalIndex. The DofMapper uses an IndexSet for
     mapping the dofs, so we can hide the real grid index behind the index
     set. Furthermore if a grid doesn't provide the consecutive index set
     then this can be calculated in the IndexSet. These two following index
     sets are just the identiy to the grid indices.

     The DefaultGridIndexSetBase defines some methods that are needed for
     index sets that cope with adaptation, but aren't needed for the following
     index set, so most of this methods do notin'.
   */
  class DefaultEmptyIndexSet
  {
    // dummy value
    enum { myType = -1 };

  protected:
    const bool adaptive_;

  public:
    //! default constructor
    DefaultEmptyIndexSet (bool adaptive) : adaptive_(adaptive) {}

    //! return false mean the no memory has to be allocated
    //! and no compress of date has to be done
    bool compress () {
      return false;
    }

    //! \brief returns true if index set can be used for adaptive
    //! calculations
    bool adaptive () const { return adaptive_; }

    //! returns true if index set gernally needs compress after adaptation
    bool needsCompress () const { return false; }

    //! do nothing here, because fathers index should already exist
    template <class EntityType>
    void insertNewIndex( const EntityType & en ) {
      assert(adaptive_) ;
    }

    //! do nothing here, because fathers index should already exist
    template <class EntityType>
    void removeOldIndex( const EntityType & en ) {
      assert(adaptive_) ;
    }

    //! nothing to do here
    void resize () {}

    //! no extra memory for restriction is needed
    int additionalSizeEstimate () const { return 0; }

    static int type() { return myType; }

    //! we have no old size
    int numberOfHoles ( int codim ) const { return 0; }

    //! return old index, for dof manager only
    int oldIndex (int hole, int codim ) const { return 0; }

    //! return new index, for dof manager only
    int newIndex (int hole, int codim ) const { return 0; }

    //! write index set to xdr file
    bool write_xdr(const std::basic_string<char> filename , int timestep)
    {
      FILE  *file;
      XDR xdrs;
      const char *path = "";

      std::basic_string<char> fnstr  = genFilename(path,filename, timestep);
      const char * fn = fnstr.c_str();
      file = fopen(fn, "wb");
      if (!file)
      {
        std::cerr << "\aERROR in DefaultGridIndexSet::write_xdr(..): could not open <"
                  << filename << ">!" << std::endl;
        return false;
      }

      xdrstdio_create(&xdrs, file, XDR_ENCODE);
      this->processXdr(&xdrs);

      xdr_destroy(&xdrs);
      fclose(file);

      return true;
    }

    //! read index set to xdr file
    bool read_xdr(const std::basic_string<char> filename , int timestep)
    {
      FILE   *file;
      XDR xdrs;
      const char *path = "";

      std::basic_string<char> fnstr = genFilename(path,filename, timestep);
      const char * fn  = fnstr.c_str();
      std::cout << "Reading <" << fn << "> \n";
      file = fopen(fn, "rb");
      if(!file)
      {
        std::cerr << "\aERROR in DefaultGridIndexSet::read_xdr(..): could not open <"
                  << filename << ">!" << std::endl;
        return(false);
      }

      // read xdr
      xdrstdio_create(&xdrs, file, XDR_DECODE);
      this->processXdr(&xdrs);

      xdr_destroy(&xdrs);
      fclose(file);
      return true;
    }

  protected:
    // read/write from/to xdr stream
    bool processXdr(XDR *xdrs)
    {
      int type = myType;
      xdr_int ( xdrs, &type);
      if(type != myType)
      {
        std::cerr << "\nERROR: DefaultGridIndex: wrong type choosen! \n\n";
        assert(type == myType);
      }
      return true;
    }
  };
  //! default base class for index set implementations for FR numerics
  template <class GridType>
  class DefaultGridIndexSetBase : public DefaultEmptyIndexSet
  {
    // dummy value
    enum { myType = -1 };
  public:
    enum { ncodim = GridType::dimension + 1 };

    //! Conschdrugdor
    DefaultGridIndexSetBase (const GridType & grid )
    // here false, because methods have to be overloaded
      : DefaultEmptyIndexSet(false)
        , grid_ (grid) {}
  protected:
    // the corresponding grid
    const GridType & grid_;
  };
  //*********************************************************************
  /*! \brief
   * DefaultLevelIndexSet generates a level index set for a grid out of a
   * grids hierarchic index set by storing for each entity in the grid
   * a number in an array.
   */
  //*********************************************************************

  template <class DefaultLevelIndexSetType, int codim>
  struct CheckLevelForCodim
  {
    static void recursive(DefaultLevelIndexSetType & d)
    {
      d.template checkLevelIndexForCodim<codim> ();
      CheckLevelForCodim<DefaultLevelIndexSetType,codim-1>::recursive(d);
    }
  };

  template <class DefaultLevelIndexSetType>
  struct CheckLevelForCodim<DefaultLevelIndexSetType,0>
  {
    static void recursive(DefaultLevelIndexSetType & d)
    {
      d.template checkLevelIndexForCodim<0> ();
    }
  };


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

  /*! \brief
     DefaultLevelIndexSet creates a LevelIndexSet for a Grid by using its
     HierarchicIndexSet
   */
  template <class GridImp>
  class DefaultLevelIndexSet :
    public IndexSet< GridImp, DefaultLevelIndexSet <GridImp>, unsigned int >

  {
    typedef GridImp GridType;
    enum { dim = GridType :: dimension };

  public:
    enum { ncodim = GridType::dimension + 1 };

    //! type of index
    typedef unsigned int IndexType;

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
    typedef DefaultLevelIndexSet<GridType> ThisType;

    template <int dummy, class Grid>
    struct ContainsIndex
    {
      static bool contains( const PersistentContainerType& container,
                            const size_t index )
      {
        DUNE_THROW(NotImplemented,"DefaultLevelIndexSet::containsIndex is only implemented for ALU3dGrid" );
        return true ;
      }
    };

    template <int dummy, ALU3dGridElementType elType, class Comm>
    struct ContainsIndex< dummy, ALU3dGrid< elType, Comm > >
    {
      static bool contains( const PersistentContainerType& container,
                            const size_t index )
      {
        return container.getData( index ).index() >= 0;
      }
    };

    template <class EntityType, int codim>
    struct InsertEntity
    {
      static void insert(const EntityType & en,
                         PersistentContainerVectorType &indexContainer,
                         int (&num)[ncodim])
      {
        PersistentContainerType& codimContainer = *(indexContainer[ codim ]);
        for( int i = 0; i < en.template count< codim >(); ++i )
        {
          Index& idx = codimContainer( en , i );
          if( idx.index() < 0 )
          {
            idx.set( num[codim] );
            ++ num[ codim ];
          }
        }
        InsertEntity<EntityType,codim-1>::insert(en, indexContainer, num);
      }
    };

    template <class EntityType>
    struct InsertEntity<EntityType,0>
    {
      static void insert(const EntityType & en,
                         PersistentContainerVectorType &indexContainer,
                         int (&num)[ncodim])
      {
        enum { codim = 0 };
        PersistentContainerType& codimContainer = *(indexContainer[ codim ]);
        Index& idx = codimContainer[ en ];
        if( idx.index() < 0 )
        {
          idx.set( num[codim] );
          ++ num[ codim ];
        }
      }
    };

  public:
    //! import default implementation of subIndex<cc>
    //! \todo remove after next release
    using IndexSet<GridType, DefaultLevelIndexSet>::subIndex;

    //! create LevelIndex by using the HierarchicIndexSet of a grid
    //! for the given level
    DefaultLevelIndexSet( const GridType & grid ,
                          const int level )
      : grid_(grid),
        indexContainers_( ncodim, (PersistentContainerType *) 0),
        level_(level)
    {
      for( int codim=0; codim < ncodim; ++codim )
        indexContainers_[ codim ] = new PersistentContainerType( grid, codim );

      calcNewIndex ();
    }

    //! desctructor deleting persistent containers
    ~DefaultLevelIndexSet ()
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
      assert( (codim == dim) ? (1) : (level_ == en.level() ));
      assert( indexContainer( codim )[ en ] >= 0 );
#endif
      return indexContainer( codim )[ en ];
    }

    //! return LevelIndex of given entity
    template<int cd>
    IndexType index (const typename GridImp::template Codim<cd>::Entity& en) const
    {
      // this must not be true for vertices
      // therefore only check other codims
#ifndef NDEBUG
      const int codim = cd;
      assert( (codim == dim) ? (1) : (level_ == en.level() ));
      assert( indexContainer( codim )[ en ].index() >= 0 );
#endif
      return indexContainer( codim )[ en ].index();
    }

    //! return subIndex (LevelIndex) for a given Entity of codim = 0 and a
    //! given SubEntity codim and number of SubEntity
    template< int cc >
    IndexType subIndex ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      assert( (codim != 0) || (level_ == e.level()) );
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
    void calcNewIndex ()
    {
      // resize arrays to new size
      for(int cd=0; cd<ncodim; ++cd)
      {
        indexContainer( cd ).clear();
      }

      // walk grid and store index
      typedef typename DefaultLevelIteratorTypes<GridImp>:: template Codim<0>::
      template Partition<All_Partition> :: Iterator IteratorType;

      // we start with zero for all codims
      for(int cd=0; cd<ncodim; ++cd) size_[cd] = 0;

      IteratorType endit  = this->template end  < 0, All_Partition > ();
      for(IteratorType it = this->template begin< 0, All_Partition > ();
          it != endit; ++it)
      {
        assert( it->level() == level_ );
        insertEntity( *it, size_ );
      }

      // remember the number of entity on level and cd = 0
      for(int cd=0; cd<ncodim; ++cd)
      {
        //if( size_[cd] != grid_.size( level_, cd) )
        //{
        //  std::cout << size_[cd] << " s | g " << grid_.size( level_, cd) << std::endl;
        //}
        assert( size_[cd] == grid_.size( level_, cd) );
      }

#ifndef NDEBUG
      CheckLevelForCodim<ThisType,dim>::recursive(*this);
#endif
    }

    // calculate index for the codim
    template <int cd>
    void checkLevelIndexForCodim ()
    {
      PersistentContainerType& codimContainer = indexContainer( cd );

      // resize memory if necessary
      // walk grid and store index
      typedef typename DefaultLevelIteratorTypes<GridImp>:: template Codim<cd>::
      template Partition<All_Partition> :: Iterator LevelIterator;

      LevelIterator endit  = this->template end  < cd , All_Partition > ();
      for(LevelIterator it = this->template begin< cd , All_Partition > (); it != endit; ++it)
      {
        assert( codimContainer[ *it ].index() >= 0 );
      }
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_.geomTypes( codim );
    }

    /** @brief Iterator to first entity of given codimension and partition type.
     */
    template<int cd, PartitionIteratorType pitype>
    typename DefaultLevelIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
    {
      return this->grid_.template lbegin<cd,pitype> (level_);
    }

    /** @brief Iterator to one past the last entity of given codim for partition type
     */
    template<int cd, PartitionIteratorType pitype>
    typename DefaultLevelIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
    {
      return this->grid_.template lend<cd,pitype> (level_);
    }

    //! returns true if this set provides an index for given entity
    bool containsIndex ( const int cd , const int idx) const
    {
      return ContainsIndex<0, GridType> :: contains( indexContainer( cd ), idx );
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
    template <class EntityType>
    void insertEntity(EntityType & en, int (&num)[ncodim])
    {
      InsertEntity<EntityType,dim>::insert( en, indexContainers_, num);
    }

    // grid this level set belongs to
    const GridType & grid_;

    //! vector with PersistentContainer for each codim
    PersistentContainerVectorType indexContainers_;

    // number of entitys of each level an codim
    int size_[ ncodim ];

    // the level for which this index set is created
    const int level_;

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


  //! Default LeafIndexSet
  template <class GridImp>
  class DefaultLeafIndexSet :
    public IndexSet< GridImp, DefaultLeafIndexSet <GridImp>, unsigned int >

  {
    typedef GridImp GridType;
    enum { dim = GridType :: dimension };

  public:
    enum { ncodim = dim + 1 };
    typedef typename GridType :: HierarchicIndexSet HierarchicIndexSetType;

    //! type of index
    typedef unsigned int IndexType;
  private:
    typedef std::vector<int> IndexArrayType;

    typedef DefaultLeafIndexSet<GridType> ThisType;

    template <class EntityType, int codim>
    struct InsertEntity
    {
      template <class HierarchicIndexSet>
      static void insert(const EntityType & en,
                         const HierarchicIndexSet & hset,
                         IndexArrayType (&index)[ncodim],
                         int (&num)[ncodim])
      {
        IndexArrayType & idx = index[codim];
        for( int i = 0; i < en.template count< codim >(); ++i )
        {
          const int id = hset.subIndex( en, i, codim );
          if( idx[ id ] < 0 )
          {
            idx[id] = num[codim];
            ++num[codim];
          }
        }
        InsertEntity<EntityType,codim-1>::insert(en,hset,index,num);
      }
    };

    template <class EntityType>
    struct InsertEntity<EntityType,0>
    {
      template <class HierarchicIndexSet>
      static void insert(const EntityType & en,
                         const HierarchicIndexSet & hset,
                         IndexArrayType (&index)[ncodim],
                         int (&num)[ncodim])
      {
        enum { codim = 0 };
        IndexArrayType & idx = index[codim];
        const int id = hset.index(en);
        if( idx[id] < 0 )
        {
          idx[id] = num[codim];
          ++num[codim];
        }
      }
    };

  public:
    //! import default implementation of subIndex<cc>
    //! \todo remove after next release
    using IndexSet<GridType, DefaultLeafIndexSet>::subIndex;

    //! create LevelIndex by using the HierarchicIndexSet of a grid
    //! for the given level
    DefaultLeafIndexSet(const GridType & grid)
      : grid_(grid)
        , hIndexSet_ ( grid.hierarchicIndexSet() )
        , size_ ( ncodim )
    {
      calcNewIndex ();
    }

    //! return LeafIndex of given entity
    template<class EntityType>
    IndexType index (const EntityType & en) const
    {
      enum { cd = EntityType :: codimension };
      // this must not be true for vertices
      // therefore only check other codims
      assert( index_[cd][ hIndexSet_.index(en) ] >= 0 );
      return index_[cd][ hIndexSet_.index(en) ];
    }

    template<int cd>
    IndexType index (const typename GridImp::template Codim<cd>::Entity& en) const
    {
      // this must not be true for vertices
      // therefore only check other codims
      assert( index_[cd][ hIndexSet_.index(en) ] >= 0 );
      return index_[cd][ hIndexSet_.index(en) ];
    }

    //! return subIndex (LevelIndex) for a given Entity of codim = 0 and a
    //! given SubEntity codim and number of SubEntity
    template< int cc >
    IndexType subIndex ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      const int hIndex = hIndexSet_.subIndex( e, i, codim );
      assert( index_[ codim ][ hIndex ] >= 0 );
      return index_[ codim ][ hIndex ];
    }

    //! returns true if this set provides an index for given entity
    template<class EntityType>
    bool contains (const EntityType& en) const
    {
      enum { cd = EntityType :: codimension };
      return (index_[cd][ hIndexSet_.index(en) ] >= 0 );
    }

    //! return size of IndexSet for a given codim
    IndexType size ( int codim ) const
    {
      assert( codim >= 0 && codim <= GridType::dimension );
      return size_[codim];
    }

    //! return size of IndexSet for a codim
    //! this method is to be revised
    IndexType size ( GeometryType type ) const
    {
      if( typeNotValid(type) ) return 0;
      return size_[GridType::dimension-type.dim()];
    }

    //! do calculation of the index set, has to be called when grid was
    //! changed or if index set is created
    void calcNewIndex ()
    {
      // resize arrays to new size
      for(int cd=0; cd<ncodim; ++cd)
      {
        resizeVectors(index_[cd], hIndexSet_.size(cd));
      }

      // walk grid and store index
      typedef typename DefaultLeafIteratorTypes<GridImp>:: template Codim<0>::
      template Partition<All_Partition> :: Iterator IteratorType;

      // we start with zero for all codims
      int num[ncodim];
      for(int cd=0; cd<ncodim; ++cd) num[cd] = 0;

      int elems=0;

      IteratorType endit  = this->template end  < 0, All_Partition > ();
      for(IteratorType it = this->template begin< 0, All_Partition > ();
          it != endit; ++it)
      {
        ++elems;
        insertEntity(*it,num);
      }

      // remember the number of entity on level and cd = 0
      for(int cd=0; cd<ncodim; ++cd)
      {
        size_[cd] = num[cd];
        //if( size_[cd] != grid_.size(cd) )
        //  std::cout << size_[cd] << " calc | grid " << grid_.size(cd)  << "\n";
        //assert( size_[cd] == grid_.size(cd) );
      }
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return hIndexSet_.geomTypes( codim );
    }


    /** @brief Iterator to first entity of given codimension and partition type.
     */
    template<int cd, PartitionIteratorType pitype>
    typename DefaultLeafIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
    {
      return this->grid_.template leafbegin<cd,pitype> ();
    }

    /** @brief Iterator to one past the last entity of given codim for partition type
     */
    template<int cd, PartitionIteratorType pitype>
    typename DefaultLeafIteratorTypes<GridImp>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
    {
      return this->grid_.template leafend<cd,pitype> ();
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
    template <class EntityType>
    void insertEntity(EntityType & en, int (&num)[ncodim])
    {
      InsertEntity<EntityType,dim>::insert(en,hIndexSet_,index_,num);
    }

    // resize vectors of index set
    void resizeVectors(IndexArrayType &a, int newNumberOfEntries)
    {
      if( newNumberOfEntries > 0 )
      {
        a.resize(newNumberOfEntries);
      }
      for(size_t i=0; i<a.size(); ++i) a[i] = -1;
    }

    // method prints indices of given codim, for debugging
    void print (int codim) const
    {
      for(size_t i=0; i<index_[codim].size(); i++)
      {
        std::cout << "levelind[" << i << "] = " << index_[codim][i] << "\n";
      }
    }

    // grid this level set belongs to
    const GridType & grid_;

    // the grids HierarchicIndexSet
    const HierarchicIndexSetType & hIndexSet_;

    // number of entitys of each level an codim
    IndexArrayType size_;

    //*********************************************************
    // Methods for mapping the hierarchic Index to index on Level
    IndexArrayType index_[ncodim];
  };

} // end namespace Dune

#endif
