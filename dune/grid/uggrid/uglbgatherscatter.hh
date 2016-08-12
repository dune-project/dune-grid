// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGLBGATHERSCATTER_HH
#define DUNE_UGLBGATHERSCATTER_HH

namespace Dune {

  /** \brief Gather/scatter methods for dynamic loadbalancing with UGGrid
   */
  class UGLBGatherScatter
  {
    class LBMessageBuffer
    {
    public:
      template <class DataType>
      void read(DataType& x)
      {
        count_--;
        memcpy(&x, data_ + count_*sizeof(DataType), sizeof(DataType));
        if (!count_)
          free(data_);
      }

      template <class DataType>
      void write(const DataType& x)
      {
        count_++;
        char* moreData_ = (char*)realloc(data_, count_*sizeof(DataType));
        if (moreData_)
          data_ = moreData_;
        memcpy(data_ + (count_-1)*sizeof(DataType), &x, sizeof(DataType));
      }

      LBMessageBuffer()
        : count_(0), data_(nullptr)
      {}

    private:
      int count_;
      char* data_;
    };

  public:

    /** \brief Gather data before load balancing
     */
    template <int codim, class GridView, class DataHandle>
    static void gather(const GridView& gridView, DataHandle& dataHandle)
    {
      const int dim = GridView::dimension;

      // do nothing if nothing has to be done
      if (!dataHandle.contains(dim, codim))
        return;

      // write the data into a global vector on process 0
      // write the macrogrid index of each entity into the corresponding UG vector
      for (const auto entity : entities(gridView, Codim<codim>()))
      {
        int numberOfParams = dataHandle.size(entity);
        if (!numberOfParams)
          continue;

        // obtain data from DUNE handle and write it into the UG message buffer
        LBMessageBuffer lbMessageBuffer;
        dataHandle.gather(lbMessageBuffer, entity);

        char*& buffer = gridView.grid().getRealImplementation(entity).getTarget()->message_buffer;
        assert(not buffer);

        typedef typename DataHandle::DataType DataType;
        buffer = (char*)malloc(sizeof(int) + numberOfParams*sizeof(DataType));
        *((int*)buffer) = numberOfParams*sizeof(DataType);       // Size of the actual payload

        for (int paramIdx = 0; paramIdx < numberOfParams; paramIdx++)
        {
          DataType *dataPointer = (DataType*)(buffer + sizeof(int) + paramIdx*sizeof(DataType));
          lbMessageBuffer.read(*dataPointer);
        }
      }
    }

    /** \brief Scatter data after load balancing
     */
    template <int codim, class GridView, class DataHandle>
    static void scatter(const GridView& gridView, DataHandle& dataHandle)
    {
      const int dim = GridView::dimension;

      // do nothing if nothing has to be done
      if (!dataHandle.contains(dim, codim))
        return;

      // obtain the data from the global vector with help of
      // the macro index and scatter it
      for (const auto entity : entities(gridView, Codim<codim>()))
      {
        int numberOfParams = dataHandle.size(entity);
        if (!numberOfParams)
          continue;

        // get data from UG message buffer and write to DUNE message buffer
        char*& buffer = gridView.grid().getRealImplementation(entity).getTarget()->message_buffer;
        assert(buffer);

        LBMessageBuffer lbMessageBuffer;

        for(int paramIdx = 0; paramIdx < numberOfParams; paramIdx++)
        {
          typedef typename DataHandle::DataType DataType;
          DataType *dataPointer = (DataType*)(buffer + sizeof(int) + paramIdx*sizeof(DataType));
          lbMessageBuffer.write(*dataPointer);
        }

        // call the data handle with the message buffer
        dataHandle.scatter(lbMessageBuffer, entity, numberOfParams);

        // free object's local message buffer
        free (buffer);
        buffer = nullptr;
      }
    }
  };

} // namespace Dune

#endif
