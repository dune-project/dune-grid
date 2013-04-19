// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGLBGATHERSCATTER_HH
#define DUNE_UGLBGATHERSCATTER_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

  /** \brief Gather/scatter methods for dynamic loadbalancing with UGGrid
   *
   * \tparam dim Grid dimension (2 or 3)
   */
  template<int dim>
  class UGLBGatherScatter
  {
    class LBMessageBuffer
    {
    public:
      template <class DataType>
      void read(DataType& x)
      {
        memcpy(&x, data_, sizeof(DataType));
        free(data_);
      }

      template <class DataType>
      void write(const DataType& x)
      {
        data_ = (char*)malloc(sizeof(DataType));
        memcpy(data_, &x, sizeof(DataType));
      }

    private:
      char* data_;
    };

    template <int commCodim>
    struct LayoutWrapper
    {
      template <int dimension>
      struct Layout
      {
        bool contains(Dune::GeometryType gt)
        { return gt.dim() == dimension - commCodim;  }
      };
    };

  public:

    /** \brief Gather data before load balancing
     * \param dataVector global vector to be filled
     */
    template <int codim, class GridView, class DataHandle>
    static void gather(const GridView& gridView,
                       DataHandle& dataHandle)
    {
      typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
          LayoutWrapper<codim>::template Layout
          > MapperType;
      MapperType mapper(gridView);

      typedef typename GridView::template Codim<codim>::Iterator Iterator;
      Iterator it = gridView.template begin<codim>();
      const Iterator& endIt = gridView.template end<codim>();
      // obtain size of the data handle
      int numberOfParams = 0;
      if (it != endIt && dataHandle.contains(dim, codim))
        numberOfParams = dataHandle.size(*it);

      // do nothing if nothing has to be done
      if (gridView.comm().max(numberOfParams) == 0)
        return;

      // write the data into a global vector on process 0
      // write the macrogrid index of each entity into the corresponding UG vector
      for (; it != endIt; ++it) {
        // obtain data from handle
        // and write it into the message buffer
        LBMessageBuffer lbMessageBuffer;
        dataHandle.gather(lbMessageBuffer, *it);

        char*& buffer = gridView.grid().getRealImplementation(*it).getTarget()->message_buffer;
        assert(not buffer);

        typedef typename DataHandle::DataType DataType;
        buffer = (char*)malloc(sizeof(int) + numberOfParams*sizeof(DataType));
        *((int*)buffer) = numberOfParams*sizeof(DataType);       // Size of the actual payload

        DataType dummy;
        lbMessageBuffer.read(dummy);
        // Copy the data into the message buffer
        memcpy(buffer+sizeof(int), &dummy, sizeof(DataType));
      }
    }

    /** \brief Scatter data after load balancing
     * \param dataVector global vector with data
     */
    template <int codim, class GridView, class DataHandle>
    static void scatter(const GridView& gridView,
                        DataHandle& dataHandle)
    {
      typedef typename GridView::template Codim<codim>::Iterator Iterator;
      Iterator it = gridView.template begin<codim>();
      const Iterator& endIt = gridView.template end<codim>();

      // obtain size of the data handle
      int numberOfParams = 0;
      if (it != endIt && dataHandle.contains(dim, codim))
        numberOfParams = dataHandle.size(*it);

      // do nothing if nothing has to be done
      if (gridView.comm().max(numberOfParams) == 0)
        return;

      // obtain the data from the global vector with help of
      // the macro index and scatter it
      for (; it != endIt; ++it) {
        char*& buffer = gridView.grid().getRealImplementation(*it).getTarget()->message_buffer;
        assert(buffer);

        // get data from global vector and write to message buffer
        typedef typename DataHandle::DataType DataType;
        DataType dummy;
        memcpy(&dummy, buffer+sizeof(int), sizeof(DataType));

        LBMessageBuffer lbMessageBuffer;
        lbMessageBuffer.write(dummy);

        // provide the data handle with the message buffer
        dataHandle.scatter(lbMessageBuffer, *it, numberOfParams);

        // free object's local message buffer
        free (buffer);
        buffer = nullptr;
      }
    }
  };

} // namespace Dune

#endif
