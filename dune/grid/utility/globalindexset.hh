// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
 *
 *  \brief Provides a globally unique index for all entities of a distributed Dune grid
 *
 * Such functionality is relevant for a number of applications:
 *  - Map a degree of freedom associated with an entity to its
 *         location in a global matrix or global vector
 *  - Such indices for elements are needed as input to external mesh partitioners
 *  - Using matrix and vector routines from the PETSc or trilinos parallel linear algebra
 *    packages for distributed memory parallel computers.
 *
 *  Method: (1) The UniqueEntityPartition class assigns an owner process to each entity
 *
 *          (2) Compute the number of entities that are owned by each process
 *
 *          (3) we communicate the index of entities that are owned by the process to processes
 *              that also contain these entities but do not own them, so that on a non-owner process
 *              we have information on the index of the entity that it got from the owner-process;
 *
 *  Copyright by Benedikt Oswald and Patrick Leidenberger, 2002-2010. All rights reserved.
 *
 *  \author    Benedikt Oswald, Patrick Leidenberger, Oliver Sander
 *
 *  \attention globally unique indices are ONLY provided for entities of the
 *             InteriorBorder_Partition type, NOT for the Ghost_Partition type !!!
 *
 *  \note The interface in this file is experimental, and may change without prior notice.
 */

#ifndef DUNE_GRID_UTILITY_GLOBALINDEXSET_HH
#define DUNE_GRID_UTILITY_GLOBALINDEXSET_HH

/** \brief Include standard header files. */
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <algorithm>

/** include base class functionality for the communication interface */
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

// Include Dune header files
#include <dune/common/version.hh>

/** include parallel capability */
#if HAVE_MPI
  #include <dune/common/parallel/mpihelper.hh>
#endif

namespace Dune
{

  /** \brief Calculate globally unique index over all processes in a Dune grid
   */
  template<class GridView, int CODIM>
  class GlobalIndexSet
  {
  public:
    /** \brief The number type used for global indices  */
    typedef int Index;

  private:
    /** define data types */
    typedef typename GridView::Grid Grid;

    typedef typename GridView::Grid::GlobalIdSet GlobalIdSet;
    typedef typename GridView::Grid::GlobalIdSet::IdType IdType;
    typedef typename GridView::Traits::template Codim<0>::Iterator Iterator;
    typedef typename GridView::Traits::template Codim<CODIM>::Entity Entity;

    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    typedef std::map<IdType,int> MapId2Index;
    typedef std::map<int,int>    IndexMap;

    /*********************************************************************************************/
    /* calculate unique partitioning for all entities of a given codim in a given GridView,      */
    /* assuming they all have the same geometry, i.e. codim, type                                */
    /*********************************************************************************************/
    class UniqueEntityPartition
    {
    private:
      /* A DataHandle class to calculate the minimum of a std::vector which is accompanied by an index set */
      template<class IS, class V> // mapper type and vector type
      class MinimumExchange
      : public Dune::CommDataHandleIF<MinimumExchange<IS,V>,typename V::value_type>
      {
      public:
        //! export type of data for message buffer
        typedef typename V::value_type DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, int codim) const
        {
          return codim==CODIM ;
        }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedsize (int dim, int codim) const
        {
          return true ;
        }

        /*! how many objects of type DataType have to be sent for a given entity
         *
         *  Note: Only the sender side needs to know this size. */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
          return 1 ;
        }

        /*! pack data from user to message buffer */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
          buff.write(v_[indexset_.index(e)]);
        }

        /** \brief Unpack data from message buffer to user
         *
         * \param n The number of objects sent by the sender
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
          DataType x;
          buff.read(x);
          if (x>=0) // other is -1 means, he does not want it
            v_[indexset_.index(e)] = std::min(x,v_[indexset_.index(e)]);
        }

        //! constructor
        MinimumExchange (const IS& indexset, V& v)
        : indexset_(indexset),
          v_(v)
        {}

      private:
        const IS& indexset_;
        V& v_;
      };

    public:
      /*! \brief Constructor needs to know the grid function space
       */
      UniqueEntityPartition (const GridView& gridview)
      : assignment_(gridview.size(CODIM))
      {
        /** extract types from the GridView data type */
        typedef typename GridView::IndexSet IndexSet;

        // assign own rank to entities that I might have
        for (auto it = gridview.template begin<0>(); it!=gridview.template end<0>(); ++it)
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
          for (int i=0; i<it->subEntities(CODIM); i++)
#else
          for (int i=0; i<it->template count<CODIM>(); i++)
#endif
          {
            assignment_[gridview.indexSet().subIndex(*it,i,CODIM)]
              = ( (it->template subEntity<CODIM>(i)->partitionType()==Dune::InteriorEntity) || (it->template subEntity<CODIM>(i)->partitionType()==Dune::BorderEntity) )
              ? gridview.comm().rank()  // set to own rank
              : - 1;   // it is a ghost entity, I will not possibly own it.
          }

        /** exchange entity index through communication */
        MinimumExchange<IndexSet,std::vector<int> > dh(gridview.indexSet(),assignment_);

        gridview.communicate(dh,Dune::All_All_Interface,Dune::ForwardCommunication);

        /* convert vector of minimum ranks to assignment vector */
        for (size_t i=0; i<assignment_.size(); i++)
          assignment_[i] = (assignment_[i] == gridview.comm().rank()) ? 1 : 0;
      }

      /** answer question if entity belongs to me, to this process */
      bool owner(size_t i)
      {
        return assignment_[i];
      }

      size_t numOwners() const
      {
        return std::accumulate(assignment_.begin(), assignment_.end(), 0);
      }

    private:
      std::vector<int> assignment_;
    };

  private:
    /* A DataHandle class to communicate the global index from the
     * owning to the non-owning entity; the class is based on the MinimumExchange
     * class in the parallelsolver.hh header file.
     */
    class IndexExchange
    : public Dune::CommDataHandleIF<IndexExchange,int>
    {
    public:
      //! returns true if data for this codim should be communicated
      bool contains (int dim, int codim) const
      {
        return codim==CODIM;
      }

      //! returns true if size per entity of given dim and codim is a constant
      bool fixedsize (int dim, int codim) const
      {
        return true;
      }

      /** \brief How many objects of type DataType have to be sent for a given entity
       *
       * \note Only the sender side needs to know this size.
       */
      template<class EntityType>
      size_t size (EntityType& e) const
      {
        return 1;
      }

      /*! pack data from user to message buffer */
      template<class MessageBuffer, class EntityType>
      void gather (MessageBuffer& buff, const EntityType& e) const
      {
        IdType id=globalidset_.id(e);

        if (CODIM==0)
          buff.write(mapid2entity_[id]);
        else
          buff.write((*mapid2entity_.find(id)).second);
      }

      /** \brief Unpack data from message buffer to user
       *
       * \param n The number of objects sent by the sender
       */
      template<class MessageBuffer, class EntityType>
      void scatter (MessageBuffer& buff, const EntityType& entity, size_t n)
      {
        int x;
        buff.read(x);

        /** only if the incoming index is a valid one,
         *  i.e. if it is greater than zero, will it be
         *  inserted as the global index; it is made
         *  sure in the upper class, i.e. GlobalIndexSet,
         *  that non-owning processes use -1 to mark an entity
         *  that they do not own.
         */
        if(x >= 0) {
          const IdType id = globalidset_.id(entity);

          if (CODIM==0)
            mapid2entity_[id] = x;
          else
          {
            mapid2entity_.erase(id);
            mapid2entity_.insert(std::make_pair(id,x));

            const int lindex = indexSet_.index(entity);
            localGlobalMap_[lindex] = x;
            globalLocalMap_[x]      = lindex;
          }
        }
      }

      //! constructor
      IndexExchange (const GlobalIdSet& globalidset, MapId2Index& mapid2entity,
                     const typename GridView::IndexSet& localIndexSet, IndexMap& localGlobal, IndexMap& globalLocal)
      : globalidset_(globalidset),
      mapid2entity_(mapid2entity),
      indexSet_(localIndexSet),
      localGlobalMap_(localGlobal),
      globalLocalMap_(globalLocal)
      {}

    private:
      const GlobalIdSet& globalidset_;
      MapId2Index& mapid2entity_;

      const typename GridView::IndexSet& indexSet_;
      IndexMap& localGlobalMap_;
      IndexMap& globalLocalMap_;
    };

  public:
    /** \brief Constructor for a given GridView
     *
     * This constructor calculates the complete set of global unique indices so that we can then
     *  later query the global index, by directly passing the entity in question.
     */
    GlobalIndexSet(const GridView& gridview)
    : gridview_(gridview)
    {
      int rank = gridview.comm().rank();
      int size = gridview.comm().size();

      const typename GridView::IndexSet& indexSet = gridview.indexSet();

      std::unique_ptr<UniqueEntityPartition> uniqueEntityPartition;
      if (CODIM!=0)
        uniqueEntityPartition = std::unique_ptr<UniqueEntityPartition>(new UniqueEntityPartition(gridview));

      nLocalEntity_ = (CODIM==0)
                    ? std::distance(gridview.template begin<0, Dune::Interior_Partition>(), gridview.template end<0, Dune::Interior_Partition>())
                    : uniqueEntityPartition->numOwners();

      // Compute the global, non-redundant number of entities, i.e. the number of entities in the set
      // without double, aka. redundant entities, on the interprocessor boundary via global reduce. */
      nGlobalEntity_ = gridview.comm().template sum<int>(nLocalEntity_);

      /* communicate the number of locally owned entities to all other processes so that the respective offset
       * can be calculated on the respective processor; we use the Dune mpi collective communication facility
       * for this; first, we gather the number of locally owned entities on the root process and, second, we
       * broadcast the array to all processes where the respective offset can be calculated. */

      std::vector<int> offset(size);
      std::fill(offset.begin(), offset.end(), 0);

      if (CODIM==0)  // This case is simpler
      {
        /** Share number of locally owned entities */
        gridview_.comm().template allgather<int>(&nLocalEntity_, 1, offset.data());
      }
      else
      {
        // gather number of locally owned entities on root process
        gridview.comm().template gather<int>(&nLocalEntity_,offset.data(),1,0);

        // broadcast the array containing the number of locally owned entities to all processes
        gridview.comm().template broadcast<int>(offset.data(),size,0);
      }

      indexOffset_.clear();

      // Do we need both cases here?
      if (CODIM==0)
        indexOffset_.resize(size + 1, 0);
      else
        indexOffset_.resize(size, 0);

      for (unsigned int ii=0; ii<indexOffset_.size(); ++ii)
        for (unsigned int jj=0; jj<ii; ++jj)
          indexOffset_[ii] += offset[jj];

      /*  compute globally unique index over all processes; the idea of the algorithm is as follows: if
       *  an entity is owned by the process, it is assigned an index that is the addition of the offset
       *  specific for this process and a consecutively incremented counter; if the entity is not owned
       *  by the process, it is assigned -1, which signals that this specific entity will get its global
       *  unique index through communication afterwards;
       *
       *  thus, the calculation of the globally unique index is divided into 2 stages:
       *
       *  (1) we calculate the global index independently;
       *
       *  (2) we achieve parallel adjustment by communicating the index
       *      from the owning entity to the non-owning entity.
       *
       */

      // 1st stage of global index calculation: calculate global index for owned entities
      // initialize map that stores an entity's global index via it's globally unique id as key
      globalIndex_.clear();

      const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet();      /** retrieve globally unique Id set */

      const int myoffset = indexOffset_[rank];

      int globalcontrib = 0;      /** initialize contribution for the global index */

      if (CODIM==0)  // This case is simpler
      {
        for (Iterator iter = gridview_.template begin<0>(); iter!=gridview_.template end<0>(); ++iter)
        {
          const IdType id = globalIdSet.id(*iter);      /** retrieve the entity's id */

          /** if the entity is owned by the process, go ahead with computing the global index */
          if (iter->partitionType() == Dune::InteriorEntity)
          {
            const int gindex = myoffset + globalcontrib;    /** compute global index */

            globalIndex_[id] = gindex;                      /** insert pair (key, datum) into the map */
            globalcontrib++;                                /** increment contribution to global index */
          }

          /** if entity is not owned, insert -1 to signal not yet calculated global index */
          else
          {
            globalIndex_[id] = -1;     /** insert pair (key, datum) into the map */
          }
        }
      }
      else  // if (CODIM==0) else
      {
      std::vector<bool> firstTime(gridview_.size(CODIM));
      std::fill(firstTime.begin(), firstTime.end(), true);

      for(Iterator iter = gridview_.template begin<0>();iter!=gridview_.template end<0>(); ++iter)
      {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
        for (size_t i=0; i<iter->subEntities(CODIM); i++)
#else
        for (size_t i=0; i<iter->template count<CODIM>(); i++)
#endif
        {
          IdType id=globalIdSet.subId(*iter,i,CODIM);                 /** retrieve the entity's id */

          int idx = gridview_.indexSet().subIndex(*iter,i,CODIM);

          if (!firstTime[idx] )
            continue;

          firstTime[idx] = false;

          if (uniqueEntityPartition->owner(idx) == true)  /** if the entity is owned by the process, go ahead with computing the global index */
          {
            const int gindex = myoffset + globalcontrib;    /** compute global index */
            globalIndex_.insert(std::make_pair(id,gindex)); /** insert pair (key, value) into the map */

            const int lindex = idx;
            localGlobalMap_[lindex] = gindex;
            globalLocalMap_[gindex] = lindex;

            globalcontrib++;                                /** increment contribution to global index */
          }
          else /** if entity is not owned, insert -1 to signal not yet calculated global index */
          {
            globalIndex_.insert(std::make_pair(id,-1));
          }
        }

      }
      }

      // 2nd stage of global index calculation: communicate global index for non-owned entities

      // Create the data handle and communicate.
      IndexExchange dataHandle(globalIdSet,globalIndex_,indexSet,localGlobalMap_,globalLocalMap_);
      gridview_.communicate(dataHandle, Dune::All_All_Interface, Dune::ForwardCommunication);
    }


    /** \brief Given a local index, retrieve its globally unique index */
    int globalIndex(const int& localIndex) const {
      return localGlobalMap_.find(localIndex)->second;
    }

    int localIndex(const int& globalIndex) const {
      return globalLocalMap_.find(globalIndex)->second;
    }

    int globalIndex(const typename GridView::template Codim<CODIM>::Entity& entity) const
    {
      if (CODIM==0)
      {
        /** global unique index is only applicable for inter or border type entities */
        const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet(); /** retrieve globally unique Id set */
        const IdType id = globalIdSet.id(entity);                        /** obtain the entity's id */
        const int gindex = globalIndex_.find(id)->second;                /** retrieve the global index in the map with the id as key */

        return gindex;
      }
      else
        return localGlobalMap_.find(gridview_.indexSet().index(entity))->second;
    }

    inline unsigned int nGlobalEntity() const
    {
      return nGlobalEntity_;
    }

    inline unsigned int nOwnedLocalEntity() const
    {
      return nLocalEntity_;
    }

    const std::vector<int>& indexOffset()
    {
      return indexOffset_;
    }

  protected:
    const GridView gridview_;

    //! Number of entities that are owned by the local process
    int nLocalEntity_;

    //! Global number of entities, i.e. number of entities without rendundant entities on interprocessor boundaries
    int nGlobalEntity_;

    //! Offset of entity index on every process
    std::vector<int> indexOffset_;

    IndexMap localGlobalMap_;
    IndexMap globalLocalMap_;

    /** \brief Stores global index of entities with entity's globally unique id as key
     */
    MapId2Index globalIndex_;
  };

}  // namespace Dune

#endif /* DUNE_GRID_UTILITY_GLOBALINDEXSET_HH */
