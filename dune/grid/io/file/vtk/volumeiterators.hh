// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_VOLUMEITERATORS_HH
#define DUNE_GRID_IO_FILE_VTK_VOLUMEITERATORS_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/io/file/vtk/corner.hh>
#include <dune/grid/io/file/vtk/corneriterator.hh>
#include <dune/grid/io/file/vtk/functionwriter.hh>
#include <dune/grid/io/file/vtk/pointiterator.hh>

namespace Dune {
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    template<typename GV>
    class ConformingVolumeIteratorFactory {
      const GV& gv;

      typedef typename GV::IndexSet IndexSet;
      const IndexSet& indexSet() const { return gv.indexSet(); }
      friend class ConformingConnectivityWriter<
          ConformingVolumeIteratorFactory<GV> >;

    public:
      static const unsigned dimCell = GV::dimension;

      typedef typename GV::template Codim<0>::Entity Cell;
      typedef typename GV::template Codim<0>::
      template Partition<InteriorBorder_Partition>::Iterator CellIterator;

      typedef VTK::Corner<Cell> Corner;
      typedef VTK::CornerIterator<CellIterator> CornerIterator;

      typedef VTK::Corner<Cell> Point;
      typedef VTK::PointIterator<CellIterator,
          typename GV::IndexSet> PointIterator;

      typedef ConformingConnectivityWriter<ConformingVolumeIteratorFactory<GV>
          > ConnectivityWriter;
      typedef typename GV::Communication Communication;

      explicit ConformingVolumeIteratorFactory(const GV& gv_)
        : gv(gv_)
      { }

      CellIterator beginCells() const {
        return gv.template begin<0, InteriorBorder_Partition>();
      }
      CellIterator endCells() const {
        return gv.template end<0, InteriorBorder_Partition>();
      }

      CornerIterator beginCorners() const {
        return CornerIterator(beginCells(), endCells());
      }
      CornerIterator endCorners() const {
        return CornerIterator(endCells());
      }

      PointIterator beginPoints() const {
        return PointIterator(beginCells(), endCells(), gv.indexSet());
      }
      PointIterator endPoints() const {
        return PointIterator(endCells());
      }

      ConnectivityWriter makeConnectivity() const {
        return ConnectivityWriter(*this);
      }
      const Communication& comm() const {
        return gv.comm();
      }
    };

    template<typename GV>
    class NonConformingVolumeIteratorFactory {
      const GV& gv;

    public:
      static const unsigned dimCell = GV::dimension;

      typedef typename GV::template Codim<0>::Entity Cell;
      typedef typename GV::template Codim<0>::
      template Partition<InteriorBorder_Partition>::Iterator CellIterator;

      typedef VTK::Corner<Cell> Corner;
      typedef VTK::CornerIterator<CellIterator> CornerIterator;

      typedef Corner Point;
      typedef CornerIterator PointIterator;

      typedef NonConformingConnectivityWriter<Cell> ConnectivityWriter;
      typedef typename GV::Communication Communication;

      explicit NonConformingVolumeIteratorFactory(const GV& gv_)
        : gv(gv_)
      { }

      CellIterator beginCells() const {
        return gv.template begin<0, InteriorBorder_Partition>();
      }
      CellIterator endCells() const {
        return gv.template end<0, InteriorBorder_Partition>();
      }

      CornerIterator beginCorners() const {
        return CornerIterator(beginCells(), endCells());
      }
      CornerIterator endCorners() const {
        return CornerIterator(endCells());
      }

      PointIterator beginPoints() const { return beginCorners(); }
      PointIterator endPoints() const { return endCorners(); }

      ConnectivityWriter makeConnectivity() const {
        return ConnectivityWriter();
      }
      const Communication& comm() const {
        return gv.comm();
      }
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_VOLUMEITERATORS_HH
