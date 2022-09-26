// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGLBGATHERSCATTER_HH
#define DUNE_UGLBGATHERSCATTER_HH

namespace Dune {

  /** \brief Gather/scatter methods for dynamic loadbalancing with UGGrid
   */
  class UGLBGatherScatter
  {
    // a last-in, first-out message buffer/queue
    template <class DataType>
    class LBMessageBuffer
    {
      friend class UGLBGatherScatter;
      static constexpr std::size_t size = sizeof(DataType);
    public:
      void read(DataType& x)
      {
        assert(bufferSize_ >= size);
        memcpy(&x, data_+bufferSize_-size, size);
        // decrease size without shrinking the allocated buffer
        bufferSize_ -= size;
        // if all item are read, deallocate the buffer
        if (bufferSize_ == 0)
          resize(0);
      }

      void write(const DataType& x)
      {
        resize(bufferSize_+size);
        memcpy(data_+bufferSize_-size, &x, size);
      }

      LBMessageBuffer()
        : bufferSize_(0), data_(nullptr)
      {}

    private:
      // copy num chars from this buffer to the ug buffer
      void copyDuneIntoUGBuffer(char* ugBuffer, std::size_t num)
      {
        // last-in, first-out
        assert(bufferSize_ >= num);
        for (std::size_t i = 0; i < num; i += size)
          std::copy(data_+bufferSize_-i-size, data_+bufferSize_-i, ugBuffer+i);
        resize(bufferSize_-num);
      }

      // copy num chars from the ug buffer to this buffer
      void copyUGIntoDuneBuffer(const char* ugBuffer, std::size_t num)
      {
        resize(bufferSize_+num);
        std::copy(ugBuffer, ugBuffer+num, data_+bufferSize_-num);
      }

      // resize internal buffer to newSize
      void resize(std::size_t newSize)
      {
        if (newSize == 0)
        {
          std::free(data_);
          data_ = nullptr;
        }
        else
        {
          if (void* mem = std::realloc(data_, newSize))
            data_ = static_cast<char*>(mem);
          else
            throw std::bad_alloc();
        }
        bufferSize_ = newSize;
      }

      std::size_t bufferSize_;
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
      for (const auto &entity : entities(gridView, Codim<codim>()))
      {
        int numberOfParams = dataHandle.size(entity);
        if (!numberOfParams)
          continue;

        // obtain data from DUNE handle and write it into the UG message buffer
        using DataType = typename DataHandle::DataType;
        LBMessageBuffer<DataType> lbMessageBuffer;
        dataHandle.gather(lbMessageBuffer, entity);

        auto ugEntity = entity.impl().getTarget();
        assert(not ugEntity->message_buffer());

        std::size_t buffer_size = numberOfParams * sizeof(DataType);
        char* buffer = static_cast<char*>(std::malloc(buffer_size));
        ugEntity->message_buffer(buffer, buffer_size);
        lbMessageBuffer.copyDuneIntoUGBuffer(buffer, buffer_size);
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
      for (const auto &entity : entities(gridView, Codim<codim>()))
      {
        // get data from UG message buffer and write to DUNE message buffer
        auto ugEntity = entity.impl().getTarget();

        typedef typename DataHandle::DataType DataType;
        auto numberOfParams = ugEntity->message_buffer_size() / sizeof(DataType);
        if (!numberOfParams)
          continue;

        auto buffer = ugEntity->message_buffer();
        assert(buffer);

        LBMessageBuffer<DataType> lbMessageBuffer;
        lbMessageBuffer.copyUGIntoDuneBuffer(buffer, numberOfParams*sizeof(DataType));

        // call the data handle with the message buffer
        dataHandle.scatter(lbMessageBuffer, entity, numberOfParams);

        // free object's local message buffer
        ugEntity->message_buffer_free();
      }
    }
  };

} // namespace Dune

#endif
