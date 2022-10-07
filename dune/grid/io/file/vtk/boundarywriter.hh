// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_BOUNDARYWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_BOUNDARYWRITER_HH

#include <memory>
#include <string>

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
      void addCellData(const std::shared_ptr<Func>& p, const std::string& name) {
        addCellData(std::shared_ptr<typename Base::FunctionWriter>
                      (new SkeletonFunctionWriter<Func>(p, name)));
      }

      template<typename Func>
      void addCellData(Func* p, const std::string& name) {
        addCellData(std::shared_ptr<Func>(p), name);
      }

      using Base::addPointData;

      template<typename Func>
      void addPointData(const std::shared_ptr<Func>& p, const std::string& name) {
        addPointData(std::shared_ptr<typename Base::FunctionWriter>
                       (new SkeletonFunctionWriter<Func>(p, name)));
      }

      template<typename Func>
      void addPointData(Func* p, const std::string& name) {
        addPointData(std::shared_ptr<Func>(p), name);
      }

    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_BOUNDARYWRITER_HH
