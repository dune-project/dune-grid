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
    // be send
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
      if ((level == -1 && UG_NS<gridDim>::isLeaf(ugEP)) || entity.level() == level)
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
      if ((level == -1 && UG_NS<gridDim>::isLeaf(ugEP)) || entity.level() == level)
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

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize(const GridView &gv)
    {
      // If the data handle claims to have the same number of data items for each entity
      // of the given codimension, then just ask for that numbering giving the first entity.
      if (duneDataHandle_->fixedSize(gridDim, codim))
      {
        auto element = gv.template begin<0, InteriorBorder_Partition>();
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
