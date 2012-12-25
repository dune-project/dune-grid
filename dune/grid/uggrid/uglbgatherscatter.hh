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
        x = data_[0];
        data_.erase(data_.begin());
      }

      template <class DataType>
      void write(const DataType& x)
      {
        data_.push_back(x);
      }

    private:
      std::vector<double> data_;
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
    template <int codim, class GridView, class DataHandle, class DataVector>
    static void gather(const GridView& gridView,
                       DataHandle& dataHandle,
                       DataVector& dataVector)
    {
      typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
          LayoutWrapper<codim>::template Layout
          > MapperType;
      MapperType mapper(gridView);

      typedef typename GridView::template Codim<codim>::Iterator Iterator;
      Iterator it = gridView.template begin<codim>();
      // obtain size of the data handle
      int numberOfParams = dataHandle.size(*it);

      // do nothing if nothing has to be done
      if (gridView.comm().max(numberOfParams) == 0)
        return;

      // resize the data vector according to the data handle
      dataVector.resize(gridView.size(codim)*numberOfParams);

      // write the data into a global vector on process 0
      // write the macrogrid index of each entity into the corresponding UG vector
      const Iterator& endIt = gridView.template end<codim>();
      for (; it != endIt; ++it) {
        // obtain data from handle
        // and write it into the message buffer
        LBMessageBuffer lbMessageBuffer;
        dataHandle.gather(lbMessageBuffer, *it);

        // assign data to the global vector
        int macroIdx = mapper.map(*it);
        for (int i = 0; i < numberOfParams; i++)
          lbMessageBuffer.read(dataVector[macroIdx*numberOfParams + i]);

        // set value of UG entity vector to the macro grid index
        UG_NS<dim>::EntityVector(GridView::Grid::getRealImplementation(*it).getTarget())->value[0] = macroIdx;
      }

      // broadcast the global data vector
      int dataSize = dataVector.size();
      gridView.comm().broadcast(&dataSize, 1, 0);
      dataVector.resize(dataSize);
      gridView.comm().broadcast(&dataVector[0], dataVector.size(), 0);
    }

    /** \brief Scatter data after load balancing
     * \param dataVector global vector with data
     */
    template <int codim, class GridView, class DataHandle, class DataVector>
    static void scatter(const GridView& gridView,
                        DataHandle& dataHandle,
                        const DataVector& dataVector)
    {
      typedef typename GridView::template Codim<codim>::Iterator Iterator;
      Iterator it = gridView.template begin<codim>();
      // obtain size of the data handle
      int numberOfParams = dataHandle.size(*it);

      // do nothing if nothing has to be done
      if (gridView.comm().max(numberOfParams) == 0)
        return;

      // write the data into a global vector on process 0
      // write the macrogrid index of each entity into the corresponding UG vector
      const Iterator& endIt = gridView.template end<codim>();
      for (; it != endIt; ++it) {
        // get macrogrid index from element vector
        int macroIdx = UG_NS<dim>::EntityVector(GridView::Grid::getRealImplementation(*it).getTarget())->value[0];

        // get data from global vector and write to message buffer
        LBMessageBuffer lbMessageBuffer;
        for (int i = 0; i < numberOfParams; i++)
          lbMessageBuffer.write(dataVector[macroIdx*numberOfParams + i]);

        // provide the data handle with the message buffer
        dataHandle.scatter(lbMessageBuffer, *it, numberOfParams);
      }
    }
  };

} // namespace Dune

#endif
