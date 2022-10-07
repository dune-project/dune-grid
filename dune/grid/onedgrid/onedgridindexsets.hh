// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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

namespace Dune {

  template<class GridImp>
  class OneDGridLevelIndexSet : public IndexSet<GridImp,OneDGridLevelIndexSet<GridImp> >
  {
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
      return e.impl().levelIndex();
    }

    //! get index of subentity of an entity
    template<int cc>
    unsigned int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e,
                           int i,
                           unsigned int codim) const
    {
      return e.impl().subLevelIndex(i,codim);
    }

    //! get number of entities of given type and on this level
    std::size_t size (GeometryType type) const
    {
      return grid_->size(level_,type);
    }

    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const
    {
      return grid_->size(level_,codim);
    }

    std::vector< GeometryType > types ( int codim ) const { return myTypes_[ codim ]; }

    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return myTypes_[codim];
    }

    /** \brief Return true if e is contained in the index set.

       \warning This implementation takes O(n) time!  It also assumes that e belongs
       to the correct grid.
     */
    template <class EntityType>
    bool contains (const EntityType& e) const
    {
      // If the entity isn't on the same level as this index set it cannot be contained in the index set
      if (e.level() != level_)
        return false;

      constexpr static int cd = EntityType::codimension;
      typedef typename GridImp::template Codim<cd>::template Partition<All_Partition>::LevelIterator IteratorType;
      IteratorType iend = grid_->levelGridView(level_).template end<cd,All_Partition>();
      for (IteratorType it = grid_->levelGridView(level_).template begin<cd,All_Partition>();
           it != iend; ++it)
      {
        if (this->template index<cd>(*it) == this->template index<cd>(e))
          return true;
      }
      return false;
    }

    /** \brief Sets the corresponding internal fields.  Used by the GridFactory */
    void setSizesAndTypes(unsigned int numVertices, unsigned int numElements) {

      numVertices_ = numVertices;
      numElements_ = numElements;

      // ///////////////////////////////////////////////
      //   Update the list of geometry types present
      // ///////////////////////////////////////////////
      if (numElements_>0) {
        myTypes_[0].resize(1);
        myTypes_[0][0] = GeometryTypes::line;
      } else
        myTypes_[0].resize(0);

      if (numVertices_>0) {
        myTypes_[1].resize(1);
        myTypes_[1][0] = GeometryTypes::vertex;
      } else
        myTypes_[1].resize(0);

    }

    /** \todo Should be private */
    void update() {

      // ///////////////////////////////
      //   Init the element indices
      // ///////////////////////////////
      numElements_ = 0;
      OneDGridList<OneDEntityImp<1> >::const_iterator eIt;
      for (eIt = grid_->elements(level_).begin(); eIt != grid_->elements(level_).end(); eIt = eIt->succ_)
        /** \todo Remove this const cast */
        const_cast<OneDEntityImp<1>*>(eIt)->levelIndex_ = numElements_++;

      // //////////////////////////////
      //   Init the vertex indices
      // //////////////////////////////

      numVertices_ = 0;
      OneDGridList<OneDEntityImp<0> >::const_iterator vIt;
      for (vIt = grid_->vertices(level_).begin(); vIt != grid_->vertices(level_).end(); vIt = vIt->succ_)
        /** \todo Remove this const cast */
        const_cast<OneDEntityImp<0>*>(vIt)->levelIndex_ = numVertices_++;

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

  template<class GridImp>
  class OneDGridLeafIndexSet :
    public IndexSet<GridImp,OneDGridLeafIndexSet<GridImp> >
  {
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
    int index (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return e.impl().leafIndex();
    }

    //! get index of subentity of an entity
    template<int cc>
    int subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e,
                  int i,
                  unsigned int codim) const
    {
      return e.impl().subLeafIndex(i,codim);
    }

    //! get number of entities of given codim, type on the leaf level
    std::size_t size (GeometryType type) const
    {
      if (type.isVertex()) {

        return numVertices_;

      } else if (type.isLine()) {

        return numElements_;

      }

      return 0;
    }

    //! get number of entities of given codim, type on the leaf level
    std::size_t size(int codim) const
    {
      if (codim==1) {

        return numVertices_;

      } else if (codim==0) {

        return numElements_;

      }

      return 0;
    }

    std::vector< GeometryType > types ( int codim ) const { return myTypes_[ codim ]; }

    /** deliver all geometry types used in this grid */
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
      constexpr static int cd = EntityType::codimension;
      for (const auto& entity : entities(grid_.leafGridView(), Codim<cd>()))
      {
        if (entity.level() == e.level() && this->template index<cd>(entity) == this->template index<cd>(e))
          return true;
      }
      return false;
    }

    /** \brief Sets the corresponding internal fields.  Used by the GridFactory */
    void setSizesAndTypes(unsigned int numVertices, unsigned int numElements) {

      numVertices_ = numVertices;
      numElements_ = numElements;

      // ///////////////////////////////////////////////
      //   Update the list of geometry types present
      // ///////////////////////////////////////////////
      if (numElements_>0) {
        myTypes_[0].resize(1);
        myTypes_[0][0] = GeometryTypes::line;
      } else
        myTypes_[0].resize(0);

      if (numVertices_>0) {
        myTypes_[1].resize(1);
        myTypes_[1][0] = GeometryTypes::vertex;
      } else
        myTypes_[1].resize(0);

    }

    /** \todo Should be private */
    void update() {

      // ///////////////////////////////
      //   Init the element indices
      // ///////////////////////////////
      numElements_ = 0;
      for (const auto& element : elements(grid_.leafGridView()))
        element.impl().target_->leafIndex_ = numElements_++;

      // //////////////////////////////
      //   Init the vertex indices
      // //////////////////////////////

      numVertices_ = 0;

      for (int i=grid_.maxLevel(); i>=0; i--) {

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
    typedef unsigned int IdType;

    //! constructor stores reference to a grid
    OneDGridIdSet (const GridImp& g) : grid_(g) {}

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return e.impl().globalId();
    }

    //! get id of subentity
    IdType subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e,
                  int i,
                  unsigned int codim) const
    {
      return e.impl().subId(i,codim);
    }

  private:

    const GridImp& grid_;
  };

}  // namespace Dune


#endif
