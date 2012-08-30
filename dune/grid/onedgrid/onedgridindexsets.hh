// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_INDEXSETS_HH
#define DUNE_ONEDGRID_INDEXSETS_HH

/** \file
    \brief The index and id sets for the OneDGrid class
 */

#include <vector>

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/onedgrid/onedgridlist.hh>
#include <dune/grid/onedgrid/onedgridentity.hh>

namespace Dune
{

  template< class GridImp >
  class OneDGridLevelIndexSet
    : public IndexSet< GridImp, OneDGridLevelIndexSet< GridImp > >
  {
  public:
    /** \brief Constructor for a given level of a given grid
     */
    OneDGridLevelIndexSet ( const GridImp &grid, int level )
      : grid_( &grid ),
        level_( level )
    {}

    //! get index of an entity
    template<int cd>
    int index (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      return grid_->getRealImplementation(e).levelIndex();
    }

    //! get index of subentity of an entity
    template<int cc>
    unsigned int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e,
                           int i,
                           unsigned int codim) const
    {
      if( cc == 0 )
        return grid_->getRealImplementation(e).subLevelIndex(i,codim);
      else
        return grid_->getRealImplementation(e).levelIndex();
    }

    //! get number of entities of given type and on this level
    int size (GeometryType type) const
    {
      return grid_->size(level_,type);
    }

    //! get number of entities of given codim, type and on this level
    int size (int codim) const
    {
      return grid_->size(level_,codim);
    }

    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return myTypes_[codim];
    }

    /** \brief Return true if e is contained in the index set.

       Warning: this implementation takes O(n) time!  It also assumes that e belongs
       to the correct grid.
     */
    template <class EntityType>
    bool contains (const EntityType& e) const
    {
      // If the entity isn't on the same level as this index set it cannot be contained in the index set
      if (e.level() != level_)
        return false;

      enum { cd = EntityType::codimension };
      typedef typename GridImp::template Codim<cd>::template Partition<All_Partition>::LevelIterator IteratorType;
      IteratorType iend = grid_->template lend<cd,All_Partition>(level_);
      for (IteratorType it = grid_->template lbegin<cd,All_Partition>(level_); it != iend; ++it)
      {
        if (this->template index<cd>(*it) == this->template index<cd>(e))
          return true;
      }
      return false;
    }

    /** \brief Sets the corresponding internal fields.  Used by the GridFactory */
    void setSizesAndTypes(unsigned int numVertices, unsigned int numElements)
    {
      numVertices_ = numVertices;
      numElements_ = numElements;

      // ///////////////////////////////////////////////
      //   Update the list of geometry types present
      // ///////////////////////////////////////////////
      if (numElements_>0)
      {
        myTypes_[0].resize(1);
        myTypes_[0][0] = GeometryType(1);
      }
      else
        myTypes_[0].resize(0);

      if (numVertices_>0)
      {
        myTypes_[1].resize(1);
        myTypes_[1][0] = GeometryType(0);
      }
      else
        myTypes_[1].resize(0);
    }

    /** \todo Should be private */
    void update ()
    {
      // ///////////////////////////////
      //   Init the element indices
      // ///////////////////////////////
      numElements_ = 0;
      OneDGridList<OneDEntityImp<1> >::const_iterator eIt;
      for (eIt = grid_->elements(level_).begin(); eIt != grid_->elements(level_).end(); eIt = eIt->succ_)
      {
        /** \todo Remove this const cast */
        const_cast<OneDEntityImp<1>*>(eIt)->levelIndex_ = numElements_++;
      }

      // //////////////////////////////
      //   Init the vertex indices
      // //////////////////////////////

      numVertices_ = 0;
      OneDGridList<OneDEntityImp<0> >::const_iterator vIt;
      for (vIt = grid_->vertices(level_).begin(); vIt != grid_->vertices(level_).end(); vIt = vIt->succ_)
      {
        /** \todo Remove this const cast */
        const_cast<OneDEntityImp<0>*>(vIt)->levelIndex_ = numVertices_++;
      }

      // set the list of geometry types
      setSizesAndTypes(numVertices_, numElements_);
    }

  private:
    const GridImp* grid_;
    int level_;

    int numElements_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
    std::vector<GeometryType> myTypes_[2];
  };



  template< class GridImp >
  class OneDGridLeafIndexSet
    : public IndexSet< GridImp, OneDGridLeafIndexSet< GridImp > >
  {
  public:
    //! constructor stores reference to a grid and level
    OneDGridLeafIndexSet ( const GridImp &g )
      : grid_( g )
    {}

    //! get index of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    int index (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return grid_.getRealImplementation(e).leafIndex();
    }

    //! get index of subentity of an entity
    template<int cc>
    int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e,
                  int i,
                  unsigned int codim) const
    {
      if( cc == 0 )
        return grid_.getRealImplementation(e).subLeafIndex(i,codim);
      else
        return grid_.getRealImplementation(e).leafIndex();
    }

    //! get number of entities of given codim, type on the leaf level
    int size ( GeometryType type ) const
    {
      if( type.isVertex() )
        return numVertices_;
      else if( type.isLine() )
        return numElements_;
      else
        return 0;
    }

    //! get number of entities of given codim, type on the leaf level
    int size ( int codim ) const
    {
      switch( codim )
      {
      case 0 :
        return numElements_;

      case 1 :
        return numVertices_;

      default :
        return 0;
      }
    }

    /** deliver all geometry types used in this grid */
    const std::vector< GeometryType > &geomTypes ( int codim ) const
    {
      return myTypes_[ codim ];
    }

    /** \brief Return true if e is contained in the index set.

       Warning: this implementation takes O(n) time!  It also assumes that e belongs
       to the correct grid.
     */
    template <class EntityType>
    bool contains (const EntityType& e) const
    {
      enum { cd = EntityType::codimension };
      typedef typename GridImp::template Codim<cd>::template Partition<All_Partition>::LeafIterator IteratorType;
      IteratorType iend = grid_.template leafend<cd,All_Partition>();
      for( IteratorType it = grid_.template leafbegin<cd,All_Partition>(); it != iend; ++it )
      {
        if (it->level() == e.level() && this->template index<cd>(*it) == this->template index<cd>(e))
          return true;
      }
      return false;
    }

    /** \brief Sets the corresponding internal fields.  Used by the GridFactory */
    void setSizesAndTypes(unsigned int numVertices, unsigned int numElements)
    {
      numVertices_ = numVertices;
      numElements_ = numElements;

      // ///////////////////////////////////////////////
      //   Update the list of geometry types present
      // ///////////////////////////////////////////////
      if (numElements_>0) {
        myTypes_[0].resize(1);
        myTypes_[0][0] = GeometryType(1);
      } else
        myTypes_[0].resize(0);

      if (numVertices_>0) {
        myTypes_[1].resize(1);
        myTypes_[1][0] = GeometryType(0);
      } else
        myTypes_[1].resize(0);
    }

    /** \todo Should be private */
    void update ()
    {
      // ///////////////////////////////
      //   Init the element indices
      // ///////////////////////////////
      numElements_ = 0;
      typename GridImp::Traits::template Codim<0>::LeafIterator eIt    = grid_.template leafbegin<0>();
      typename GridImp::Traits::template Codim<0>::LeafIterator eEndIt = grid_.template leafend<0>();

      for (; eIt!=eEndIt; ++eIt) {

        grid_.getRealImplementation(*eIt).target_->leafIndex_ = numElements_++;

      }

      // //////////////////////////////
      //   Init the vertex indices
      // //////////////////////////////

      numVertices_ = 0;

      for (int i=grid_.maxLevel(); i>=0; i--)
      {
        const OneDEntityImp<0>* vIt;
        for (vIt = grid_.vertices(i).begin(); vIt != grid_.vertices(i).end(); vIt = vIt->succ_) {

          /** \todo Remove the const casts */
          if (vIt->isLeaf())
            const_cast<OneDEntityImp<0>*>(vIt)->leafIndex_ = numVertices_++;
          else
            const_cast<OneDEntityImp<0>*>(vIt)->leafIndex_ = vIt->son_->leafIndex_;
        }
      }

      // set the list of geometry types
      setSizesAndTypes(numVertices_, numElements_);
    }

  private:
    const GridImp &grid_;

    int numElements_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
    std::vector<GeometryType> myTypes_[2];
  };



  template< class Grid >
  class OneDGridIdSet
    : public IdSet< Grid, OneDGridIdSet< Grid >, unsigned int >
  {
    typedef OneDGridIdSet< Grid > This;
    typedef IdSet< Grid, This, unsigned int > Base;

    /* We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    //! define the type used for persistent indices
    typedef typename Base::IdType IdType;

    //! get id of an entity
    template< int cd >
    IdType id ( const typename Traits::template Codim< cd >::Entity &e ) const
    {
      return Grid::getRealImplementation( e ).globalId();
    }

    //! get id of subentity
    template< int cd >
    IdType subId ( const typename Traits::template Codim< cd >::Entity &e, int i, unsigned int codim ) const
    {
      return Grid::getRealImplementation( e ).subId( i, codim );
    }
  };

}  // namespace Dune

#endif // #ifndef DUNE_ONEDGRID_INDEXSETS_HH
