// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UG_MESSAGE_BUFFER_HH
#define UG_MESSAGE_BUFFER_HH

#include <algorithm>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/gridenums.hh>

namespace Dune {

  /** converts the UG speak message buffers to DUNE speak and vice-versa */
  template <class DataHandle, int gridDim, int codim>
  class UGMessageBuffer
  {
    typedef UGMessageBuffer<DataHandle, gridDim, codim>  ThisType;
    typedef UGGrid<gridDim>                              GridType;
    typedef typename DataHandle::DataType DataType;

    using Entity = typename UGGrid<gridDim>::template Codim<codim>::Entity;
    using EntityImp = typename Entity::Implementation;

    UGMessageBuffer(void *ugData)
    {
      ugData_ = static_cast<char*>(ugData);
    };

  public:
    void write(const DataType &t)
    { this->writeRaw_<DataType>(t);  }

    void read(DataType &t)
    { this->readRaw_<DataType>(t);  }

  private:
    friend class Dune::UGGrid<gridDim>;

    template <class ValueType>
    void writeRaw_(const ValueType &v)
    {
      std::copy(
        reinterpret_cast<const char*>(&v),
        reinterpret_cast<const char*>(&v) + sizeof(ValueType),
        reinterpret_cast<char*>(ugData_)
        );
      ugData_ += sizeof(ValueType);
    }

    template <class ValueType>
    void readRaw_(ValueType &v)
    {
      std::copy(
        reinterpret_cast<const char*>(ugData_),
        reinterpret_cast<const char*>(ugData_) + sizeof(ValueType),
        reinterpret_cast<char*>(&v)
        );
      ugData_ += sizeof(ValueType);
    }

    // called by DDD_IFOneway to serialize the data structure to
    // be sent
    static int ugGather_(
#if DUNE_UGGRID_HAVE_DDDCONTEXT
      DDD::DDDContext&,
#endif
      typename UG_NS<gridDim>::DDD_OBJ obj, void* data)
    {
      // cast the DDD object to a UG entity pointer
      auto ugEP = reinterpret_cast<typename Dune::UG_NS<gridDim>::template Entity<codim>::T*>(obj);

      // construct a DUNE entity from the UG entity pointer
      Entity entity(EntityImp(ugEP, grid_));

      // safety check to only communicate what is needed
      if ((level == -1 && isToplevelLeaf(ugEP)) || entity.level() == level)
      {
        ThisType msgBuf(static_cast<DataType*>(data));
        if (!duneDataHandle_->fixedSize(gridDim, codim))
          msgBuf.template writeRaw_<unsigned>(duneDataHandle_->size(entity));
        duneDataHandle_->gather(msgBuf, entity);
      }

      return 0;
    }

    // called by DDD_IFOneway to deserialize the data structure
    // that has been received
    static int ugScatter_(
#if DUNE_UGGRID_HAVE_DDDCONTEXT
      DDD::DDDContext&,
#endif
      typename UG_NS<gridDim>::DDD_OBJ obj, void* data)
    {
      // cast the DDD object to a UG entity pointer
      auto ugEP = reinterpret_cast<typename Dune::UG_NS<gridDim>::template Entity<codim>::T*>(obj);

      // construct a DUNE entity from the UG entity pointer
      Entity entity(EntityImp(ugEP, grid_));

      // safety check to only communicate what is needed
      if ((level == -1 && isToplevelLeaf(ugEP)) || entity.level() == level)
      {
        ThisType msgBuf(static_cast<DataType*>(data));
        int size;
        if (!duneDataHandle_->fixedSize(gridDim, codim))
          msgBuf.readRaw_(size);
        else
          size = duneDataHandle_->size(entity);
        if (size > 0)
          duneDataHandle_->scatter(msgBuf, entity, size);

      }

      return 0;
    }

    /** \brief Test whether entity is leaf and has no copies on higher (i.e., finer) levels
     *
     * Entities of the leaf grid view are equivalence class of entities of the level view.
     * Therefore, the isLeaf method of the grid interface sometimes returns true for entities
     * that have copies on a finer level.  In communication, we really want to call gather
     * and scatter for these entities only once.  Therefore we use the following special
     * implementation that returns true only if an entity is leaf and it has now children.
     */
    static bool isToplevelLeaf(const typename Dune::UG_NS<gridDim>::template Entity<codim>::T* ugEP)
    {
      // safety check to only communicate what is needed
      bool isLeaf = UG_NS<gridDim>::isLeaf(ugEP);

      // Edges: Simply checking isLeaf here will return 'true' for each top-level
      // leaf and *and* its lower-level copies (if it has any).  However, we really
      // want only the highest-level copy, because otherwise gather will be called
      // more than once for that edge (interpreted as a leaf grid edge).  Therefore
      // we need a tigher criterion.
      if constexpr (gridDim-codim==1)
      {
        if (isLeaf)
        {
          // isLeaf is true, i.e. we either are a top-level leaf oder or a copy of one.
          // If we are are lower-level copy then both of our end nodes have children,
          // and these two children are connected.
          auto nodeA = ugEP->links[0].nbnode;
          auto nodeB = ugEP->links[1].nbnode;

          if (nodeA->son && nodeB->son && UG_NS<gridDim>::GetEdge(nodeA->son, nodeB->son))
            isLeaf = false;
        }
      }

      // Facets: No such handling is currently needed, because UGGridLeafIndexSet<GridImp>::update
      // does not assign leaf indices of leaf facets to their copy-fathers.
      // This may be a bug...

      return isLeaf;
    }

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize(const GridView &gv)
    {
      // If the data handle claims to have the same number of data items for each entity
      // of the given codimension, then just ask for that numbering giving the first entity.
      if (duneDataHandle_->fixedSize(gridDim, codim))
      {
        auto element = gv.template begin<0, InteriorBorder_Partition>();
        if (element == gv.template end<0, InteriorBorder_Partition>())
        {
          return 0;
        }
        return sizeof(DataType)
               * duneDataHandle_->size(element->template subEntity<codim>(0));
      }

      // iterate over all entities, find the maximum size for
      // the current rank
      std::size_t maxSize = 0;
      if constexpr (codim==gridDim || codim==0)
      {
        for (const auto& entity : entities(gv, Codim<codim>(), Dune::Partitions::all))
          maxSize = std::max(maxSize, duneDataHandle_->size(entity));
      }
      else
      {
        for (const auto& element : elements(gv, Dune::Partitions::all))
        {
          int numberOfSubentities = element.subEntities(codim);
          for (int k = 0; k < numberOfSubentities; k++)
          {
            const auto subEntity = element.template subEntity<codim>(k);

            maxSize = std::max(maxSize, duneDataHandle_->size(subEntity));
          }
        }
      }

      // find maximum size for all ranks
      maxSize = gv.comm().max(maxSize);
      if (maxSize==0)
        return 0;

      // Add the size of an unsigned integer to the actual buffer size,
      // for storing the actual number of objects for each entity.
      return sizeof(unsigned) + sizeof(DataType)*maxSize;
    }

    static const GridType* grid_;
    static DataHandle *duneDataHandle_;
    static int level;
    char *ugData_;
  };

}   // end namespace Dune

#endif  // UG_MESSAGE_BUFFER_HH
