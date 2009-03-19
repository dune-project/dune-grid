// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_INDEXSETS_HH
#define DUNE_ONEDGRID_INDEXSETS_HH

/** \file
    \brief The index and id sets for the OneDGrid class
 */

#include <vector>

namespace Dune {

  template <class GridImp>
  struct OneDGridLevelIndexSetTypes
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


  template<class GridImp>
  class OneDGridLevelIndexSet : public IndexSet<GridImp,OneDGridLevelIndexSet<GridImp> >
  {
    typedef IndexSet<GridImp,OneDGridLevelIndexSet<GridImp> > Base;
  public:

    /** \brief Constructor for a given level of a given grid
     */
    OneDGridLevelIndexSet (const GridImp& grid, int level)
      : grid_(&grid), level_(level)
    {}

    //! get index of an entity
    template<int cd>
    int index (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      return grid_->getRealImplementation(e).levelIndex();
    }

    //! get index of subentity of a codim 0 entity
    template<int cc>
    int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i) const
    {
      return grid_->getRealImplementation(e).template subLevelIndex<cc>(i);
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

    /** \todo Should be private */
    void update() {

      // ///////////////////////////////
      //   Init the element indices
      // ///////////////////////////////
      numElements_ = 0;
      OneDGridList<OneDEntityImp<1> >::const_iterator eIt;
      for (eIt = grid_->elements[level_].begin(); eIt!=grid_->elements[level_].end(); eIt = eIt->succ_)
        /** \todo Remove this const cast */
        const_cast<OneDEntityImp<1>*>(eIt)->levelIndex_ = numElements_++;

      // //////////////////////////////
      //   Init the vertex indices
      // //////////////////////////////

      numVertices_ = 0;
      OneDGridList<OneDEntityImp<0> >::const_iterator vIt;
      for (vIt = grid_->vertices[level_].begin(); vIt!=grid_->vertices[level_].end(); vIt = vIt->succ_)
        /** \todo Remove this const cast */
        const_cast<OneDEntityImp<0>*>(vIt)->levelIndex_ = numVertices_++;

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

  private:
    const GridImp* grid_;
    int level_;

    int numElements_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
    std::vector<GeometryType> myTypes_[2];
  };

  template <class GridImp>
  struct OneDGridLeafIndexSetTypes
  {
    //! The types
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator Iterator;
      };
    };
  };

  template<class GridImp>
  class OneDGridLeafIndexSet :
    public IndexSet<GridImp,OneDGridLeafIndexSet<GridImp> >
  {
    typedef IndexSet<GridImp,OneDGridLeafIndexSet<GridImp> > Base;
  public:
    //! constructor stores reference to a grid and level
    OneDGridLeafIndexSet (const GridImp& g) : grid_(g)
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

    //! get index of subentity of a codim 0 entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cc>
    int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
    {
      return grid_.getRealImplementation(e).template subLeafIndex<cc>(i);
    }

    //! get number of entities of given codim, type on the leaf level
    int size (GeometryType type) const
    {
      if (type.isVertex()) {

        return numVertices_;

      } else if (type.isLine()) {

        return numElements_;

      }

      return 0;
    }

    //! get number of entities of given codim, type on the leaf level
    int size(int codim) const
    {
      return Base::size(codim);
    }

    /** deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return myTypes_[codim];
    }

    /** \todo Should be private */
    void update() {

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

      for (int i=grid_.maxLevel(); i>=0; i--) {

        const OneDEntityImp<0>* vIt;
        for (vIt = grid_.vertices[i].begin(); vIt!=grid_.vertices[i].end(); vIt = vIt->succ_) {

          /** \todo Remove the const casts */
          if (vIt->isLeaf())
            const_cast<OneDEntityImp<0>*>(vIt)->leafIndex_ = numVertices_++;
          else
            const_cast<OneDEntityImp<0>*>(vIt)->leafIndex_ = vIt->son_->leafIndex_;

        }

      }

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

  private:

    const GridImp& grid_;

    int numElements_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
    std::vector<GeometryType> myTypes_[2];
  };


  template<class GridImp>
  class OneDGridIdSet : public IdSet<GridImp,OneDGridIdSet<GridImp>,unsigned int>
  {
  public:
    //! define the type used for persistent indices
    typedef unsigned int GlobalIdType;
    typedef unsigned int LocalIdType;

    //! constructor stores reference to a grid
    OneDGridIdSet (const GridImp& g) : grid_(g) {}

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    GlobalIdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return grid_.getRealImplementation(e).globalId();
    }

    //! get id of subentity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    GlobalIdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
    {
      return grid_.getRealImplementation(e).template subId<cd>(i);
    }

    /** \todo Should be private */
    void update() {}

  private:

    const GridImp& grid_;
  };

}  // namespace Dune


#endif
