// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UG_MESSAGE_BUFFER_HH
#define UG_MESSAGE_BUFFER_HH

#include <algorithm>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/gridenums.hh>

namespace Dune {

  /** converts the UG speak message buffers to DUNE speak and vice-versa */
  template <class DataHandle, int GridDim, int codim>
  class UGMessageBufferBase {
  protected:
    typedef UGMessageBufferBase<DataHandle, GridDim, codim>  ThisType;
    typedef UGGrid<GridDim>                              GridType;
    typedef typename DataHandle::DataType DataType;

    enum {
      dim = GridDim
    };

    using Entity = typename UGGrid<dim>::template Codim<codim>::Entity;
    using EntityImp = typename Entity::Implementation;

    UGMessageBufferBase(void *ugData)
    {
      ugData_ = static_cast<char*>(ugData);
    };

  public:
    void write(const DataType &t)
    { this->writeRaw_<DataType>(t);  }

    void read(DataType &t)
    { this->readRaw_<DataType>(t);  }

  protected:
    friend class Dune::UGGrid<dim>;

    template <class ValueType>
    void writeRaw_(const ValueType &v)
    {
      *reinterpret_cast<ValueType*>(ugData_) = v;
      ugData_ += sizeof(ValueType);
    }

    template <class ValueType>
    void readRaw_(ValueType &v)
    {
      v = *reinterpret_cast<ValueType*>(ugData_);
      ugData_ += sizeof(ValueType);
    }

    // called by DDD_IFOneway to serialize the data structure to
    // be send
    static int ugGather_(
#if DUNE_UGGRID_HAVE_DDDCONTEXT
      DDD::DDDContext&,
#endif
      typename UG_NS<dim>::DDD_OBJ obj, void* data)
    {
      // cast the DDD object to a UG entity pointer
      auto ugEP = reinterpret_cast<typename Dune::UG_NS<dim>::template Entity<codim>::T*>(obj);

      // construct a DUNE entity from the UG entity pointer
      /** \bug The nullptr argument should actually be the UGGrid object.  But that is hard to obtain here,
       * and the argument is (currently) only used for the boundarySegmentIndex method, which we don't call. */
      Entity entity(EntityImp(ugEP, nullptr));

      // safety check to only communicate what is needed
      if ((level == -1 && UG_NS<dim>::isLeaf(ugEP)) || entity.level() == level)
      {
        ThisType msgBuf(static_cast<DataType*>(data));
        if (!duneDataHandle_->fixedSize(dim, codim))
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
      typename UG_NS<dim>::DDD_OBJ obj, void* data)
    {
      // cast the DDD object to a UG entity pointer
      auto ugEP = reinterpret_cast<typename Dune::UG_NS<dim>::template Entity<codim>::T*>(obj);

      // construct a DUNE entity from the UG entity pointer
      /** \bug The nullptr argument should actually the UGGrid object.  But that is hard to obtain here,
       * and the argument is (currently) only used for the boundarySegmentIndex method, which we don't call. */
      Entity entity(EntityImp(ugEP, nullptr));

      // safety check to only communicate what is needed
      if ((level == -1 && UG_NS<dim>::isLeaf(ugEP)) || entity.level() == level)
      {
        ThisType msgBuf(static_cast<DataType*>(data));
        int size;
        if (!duneDataHandle_->fixedSize(dim, codim))
          msgBuf.readRaw_(size);
        else
          size = duneDataHandle_->size(entity);
        if (size > 0)
          duneDataHandle_->scatter(msgBuf, entity, size);

      }

      return 0;
    }

    static DataHandle *duneDataHandle_;
    static int level;
    char *ugData_;
  };

  template <class DataHandle, int GridDim, int codim>
  class UGMessageBuffer
    : public UGMessageBufferBase<DataHandle, GridDim, codim>
  {
    typedef typename DataHandle::DataType DataType;
    typedef UGMessageBufferBase<DataHandle, GridDim, codim> Base;
    enum { dim = GridDim };

  protected:
    friend class Dune::UGGrid<dim>;

    UGMessageBuffer(void *ugData)
      : Base(ugData)
    {}

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize_(const GridView &gv)
    {
      if (Base::duneDataHandle_->fixedSize(dim, codim)) {
        return sizeof(DataType)
               * Base::duneDataHandle_->size(*gv.template begin<codim,InteriorBorder_Partition>());
      }

      // iterate over all entities, find the maximum size for
      // the current rank
      int maxSize = 0;
      for (const auto& entity : entities(gv, Codim<codim>(), Dune::Partitions::all))
      {
        maxSize = std::max((int) maxSize,
                           (int) Base::duneDataHandle_->size(entity));
      }

      // find maximum size for all ranks
      maxSize = gv.comm().max(maxSize);
      if (!maxSize)
        return 0;

      // add the size of an unsigned integer to the actual
      // buffer size. (we somewhere have to store the actual
      // number of objects for each entity.)
      return sizeof(unsigned) + sizeof(DataType)*maxSize;
    }
  };

  template <class DataHandle, int GridDim, int codim>
  class UGEdgeAndFaceMessageBuffer
    : public UGMessageBufferBase<DataHandle, GridDim, codim>
  {
    enum {dim = GridDim};
    typedef typename DataHandle::DataType DataType;
    typedef UGMessageBufferBase<DataHandle, GridDim, codim> Base;
  protected:
    friend class Dune::UGGrid<dim>;

    UGEdgeAndFaceMessageBuffer(void *ugData)
      : Base(ugData)
    {}

    // returns number of bytes required for the UG message buffer
    template <class GridView>
    static unsigned ugBufferSize_(const GridView &gv)
    {
      if (Base::duneDataHandle_->fixedSize(dim, codim)) {
        auto element = gv.template begin<0, InteriorBorder_Partition>();
        return sizeof(DataType)
               * Base::duneDataHandle_->size(element->template subEntity<codim>(0));
      }

      // iterate over all entities, find the maximum size for
      // the current rank
      int maxSize = 0;
      for (const auto& element : elements(gv, Dune::Partitions::all))
      {
        int numberOfSubentities = element.subEntities(codim);
        for (int k = 0; k < numberOfSubentities; k++)
        {
          const auto subEntity = element.template subEntity<codim>(k);

          maxSize = std::max((int) maxSize,
                             (int) Base::duneDataHandle_->size(subEntity));
        }
      }

      // find maximum size for all ranks
      maxSize = gv.comm().max(maxSize);
      if (!maxSize)
        return 0;

      // add the size of an unsigned integer to the actual
      // buffer size. (we somewhere have to store the actual
      // number of objects for each entity.)
      return sizeof(unsigned) + sizeof(DataType)*maxSize;
    }
  };

  template <class DataHandle>
  class UGMessageBuffer<DataHandle, 2, 1>
    : public UGEdgeAndFaceMessageBuffer<DataHandle, 2, 1>
  {};

  template <class DataHandle>
  class UGMessageBuffer<DataHandle, 3, 2>
    : public UGEdgeAndFaceMessageBuffer<DataHandle, 3, 2>
  {};

  template <class DataHandle>
  class UGMessageBuffer<DataHandle, 3, 1>
    : public UGEdgeAndFaceMessageBuffer<DataHandle, 3, 1>
  {};

}   // end namespace Dune

#endif  // UG_MESSAGE_BUFFER_HH
