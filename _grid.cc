// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/python/grid/commops.hh>

#include <dune/python/pybind11/pybind11.h>


PYBIND11_MODULE( _grid, module )
{
  pybind11::enum_< Dune::Python::detail::CommOp > commOps( module, "CommOp" );
  commOps.value( "set", Dune::Python::detail::CommOp::set );
  commOps.value( "add", Dune::Python::detail::CommOp::add );
}
