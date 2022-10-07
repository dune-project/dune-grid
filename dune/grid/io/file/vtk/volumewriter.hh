// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_VOLUMEWRITER_HH
#define DUNE_GRID_IO_FILE_VTK_VOLUMEWRITER_HH

#include <memory>

#include <dune/grid/io/file/vtk/basicwriter.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/functionwriter.hh>
#include <dune/grid/io/file/vtk/volumeiterators.hh>

namespace Dune {
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    template<typename GV>
    class ConformingVolumeWriter
      : public ConformingVolumeIteratorFactory<GV>,
        public BasicWriter<ConformingVolumeIteratorFactory<GV> >
    {
      typedef ConformingVolumeIteratorFactory<GV> Factory;
      typedef BasicWriter<Factory> Base;

      const GV& gv;

    public:
      typedef Dune::VTKFunction< GV > VTKFunction;
      typedef std::shared_ptr<VTKFunction> VTKFunctionPtr;

      ConformingVolumeWriter(const GV& gv_)
        : Factory(gv_), Base(static_cast<const Factory&>(*this)), gv(gv_)
      { }

      using Base::addPointData;

      void addCellData(const VTKFunctionPtr& p) {
        Base::addCellData(std::shared_ptr<typename Base::FunctionWriter>
                      (new VTKFunctionWriter<VTKFunction>(p)));
      }

      void addCellData(VTKFunction* p) {
        addCellData(VTKFunctionPtr(p));
      }

      template<typename V>
      void addCellData(const V &v, const std::string &name, int ncomps=1) {
        addCellData(new P0VTKFunction<GV, V>(gv, v, name, ncomps));
      }

      void addVertexData(const VTKFunctionPtr& p) {
        addPointData(std::shared_ptr<typename Base::FunctionWriter>
                       (new VTKFunctionWriter<VTKFunction>(p)));
      }

      void addVertexData(VTKFunction* p) {
        addVertexData(VTKFunctionPtr(p));
      }

      template<typename V>
      void addVertexData(const V &v, const std::string &name, int ncomps=1) {
        addVertexData(new P1VTKFunction<GV, V>(gv, v, name, ncomps));
      }

    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_VOLUMEWRITER_HH
