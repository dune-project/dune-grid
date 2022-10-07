// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/python/grid/commops.hh>
#include <dune/python/grid/enums.hh>

#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _grid, module )
{
  // enumeration types from dune-grid

  pybind11::enum_< Dune::PartitionType > partitionType( module, "PartitionType" );
  partitionType.value( "Interior", Dune::InteriorEntity );
  partitionType.value( "Border", Dune::BorderEntity );
  partitionType.value( "Overlap", Dune::OverlapEntity );
  partitionType.value( "Front", Dune::FrontEntity );
  partitionType.value( "Ghost", Dune::GhostEntity );

  pybind11::enum_< Dune::InterfaceType > interfaceType( module, "InterfaceType" );
  interfaceType.value( "InteriorBorder_InteriorBorder", Dune::InteriorBorder_InteriorBorder_Interface );
  interfaceType.value( "InteriorBorder_All", Dune::InteriorBorder_All_Interface );
  interfaceType.value( "Overlap_OverlapFront", Dune::Overlap_OverlapFront_Interface );
  interfaceType.value( "Overlap_All", Dune::Overlap_All_Interface );
  interfaceType.value( "All_All", Dune::All_All_Interface );

  pybind11::enum_< Dune::CommunicationDirection > communicationDirection( module, "CommunicationDirection" );
  communicationDirection.value( "Forward", Dune::ForwardCommunication );
  communicationDirection.value( "Backward", Dune::BackwardCommunication );

  pybind11::enum_< Dune::VTK::OutputType > vtkOutputType( module, "OutputType" );
  vtkOutputType.value( "ascii", Dune::VTK::OutputType::ascii );
  vtkOutputType.value( "base64", Dune::VTK::OutputType::base64 );
  vtkOutputType.value( "appendedraw", Dune::VTK::OutputType::appendedraw );
  vtkOutputType.value( "appendedbase64", Dune::VTK::OutputType::appendedbase64 );

  // enumeration types added by dune-python

  pybind11::enum_< Dune::Python::detail::CommOp > commOps( module, "CommOp" );
  commOps.value( "set", Dune::Python::detail::CommOp::set );
  commOps.value( "add", Dune::Python::detail::CommOp::add );

  pybind11::enum_< Dune::Python::Marker > marker( module, "Marker" );
  marker.value( "coarsen", Dune::Python::Marker::Coarsen );
  marker.value( "keep", Dune::Python::Marker::Keep );
  marker.value( "refine", Dune::Python::Marker::Refine );

  pybind11::enum_< Dune::Python::VTKDataType > vtkDataType( module, "DataType" );
  vtkDataType.value( "CellData", Dune::Python::VTKDataType::CellData );
  vtkDataType.value( "PointData", Dune::Python::VTKDataType::PointData );
  vtkDataType.value( "CellVector", Dune::Python::VTKDataType::CellVector );
  vtkDataType.value( "PointVector", Dune::Python::VTKDataType::PointVector );

  pybind11::enum_< Dune::Python::Reader > reader( module, "reader" );
  reader.value( "dgf", Dune::Python::Reader::dgf );
  reader.value( "dgfString", Dune::Python::Reader::dgfString );
  reader.value( "gmsh", Dune::Python::Reader::gmsh );
  reader.value( "structured", Dune::Python::Reader::structured );
}
