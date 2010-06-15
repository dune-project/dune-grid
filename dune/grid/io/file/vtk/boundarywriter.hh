// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_BOUNDARYWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_BOUNDARYWRITER_HH

#include <string>

#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/vtk/basicwriter.hh>
#include <dune/grid/io/file/vtk/boundaryiterators.hh>
#include <dune/grid/io/file/vtk/skeletonfunction.hh>

namespace Dune {
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    template<typename GV>
    class NonConformingBoundaryWriter
      : public NonConformingBoundaryIteratorFactory<GV>,
        public BasicWriter<NonConformingBoundaryIteratorFactory<GV> >
    {
      typedef NonConformingBoundaryIteratorFactory<GV> Factory;
      typedef BasicWriter<Factory> Base;

      const GV& gv;

    public:
      NonConformingBoundaryWriter(const GV& gv_)
        : Factory(gv_), Base(static_cast<const Factory&>(*this)), gv(gv_)
      { }

      using Base::addCellData;

      template<typename Func>
      void addCellData(const shared_ptr<Func>& p, const std::string& name) {
        addCellData(shared_ptr<typename Base::FunctionWriter>
                      (new SkeletonFunctionWriter<Func>(p, name)));
      }

      template<typename Func>
      void addCellData(Func* p, const std::string& name) {
        addCellData(shared_ptr<Func>(p), name);
      }

      using Base::addPointData;

      template<typename Func>
      void addPointData(const shared_ptr<Func>& p, const std::string& name) {
        addPointData(shared_ptr<typename Base::FunctionWriter>
                       (new SkeletonFunctionWriter<Func>(p, name)));
      }

      template<typename Func>
      void addPointData(Func* p, const std::string& name) {
        addPointData(shared_ptr<Func>(p), name);
      }

    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_BOUNDARYWRITER_HH
