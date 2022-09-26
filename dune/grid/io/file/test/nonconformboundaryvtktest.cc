// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/boundarywriter.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/skeletonfunction.hh>
#include <dune/grid/yaspgrid.hh>

template< class GridView >
class ScalarFunction
{
public:
  typedef Dune::VTK::SkeletonFunctionTraits<GridView, typename GridView::ctype>
  Traits;

  //! return number of components
  unsigned dimRange() const { return 1; }

  void evaluate(const typename Traits::Cell& c,
                const typename Traits::Domain& xl,
                typename Traits::Range& result) const
  {
    result.resize(1, c.geometry().global(xl).two_norm());
  }
};

template< class GridView >
class VectorFunction
{
public:
  typedef Dune::VTK::SkeletonFunctionTraits<GridView, typename GridView::ctype>
  Traits;

  //! return number of components
  unsigned dimRange() const { return GridView::dimensionworld; }

  void evaluate(const typename Traits::Cell& c,
                const typename Traits::Domain& xl,
                typename Traits::Range& result) const
  {
    const Dune::FieldVector<typename Traits::DomainField,
        GridView::dimensionworld>& normal = c.unitOuterNormal(xl);
    result.assign(normal.begin(), normal.end());
  }

};

template< class GridView >
void doWrite( const GridView &gridView )
{
  constexpr static int dim = GridView :: dimension;

  Dune::VTK::NonConformingBoundaryWriter< GridView > vtk( gridView );

  std::shared_ptr<ScalarFunction<GridView> > scalarFunc
    (new ScalarFunction<GridView>);
  vtk.addCellData(scalarFunc, "cellScalar");
  vtk.addPointData(scalarFunc, "pointScalar");

  std::shared_ptr<VectorFunction<GridView> > vectorFunc
    (new VectorFunction<GridView>);
  vtk.addCellData(vectorFunc, "cellVector");
  vtk.addPointData(vectorFunc, "pointVector");

  char name[256];
  snprintf(name,256,"nonconformboundaryvtktest-%iD-ascii", dim);
  vtk.write(name, Dune::VTK::ascii);

  snprintf(name,256,"nonconformboundaryvtktest-%iD-base64", dim);
  vtk.write(name, Dune::VTK::base64);

  snprintf(name,256,"nonconformboundaryvtktest-%iD-appendedraw", dim);
  vtk.write(name, Dune::VTK::appendedraw);

  snprintf(name,256,"nonconformboundaryvtktest-%iD-appendedbase64", dim);
  vtk.write(name, Dune::VTK::appendedbase64);
}

template<int dim>
void vtkCheck(const std::array<int,dim>& n,
              const Dune::FieldVector<double,dim>& h)
{
  std::cout << std::endl << "vtkCheck dim=" << dim << std::endl << std::endl;
  Dune::YaspGrid<dim> grid(h, n);
  grid.globalRefine(1);

  doWrite( grid.leafGridView() );
  doWrite( grid.levelGridView( 0 ) );
  doWrite( grid.levelGridView( grid.maxLevel() ) );
}

int main(int argc, char **argv)
{
  try
  {
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if(mpiHelper.rank() == 0)
      std::cout << "subsamplingvtktest: MPI_Comm_size == " << mpiHelper.size()
                << std::endl;

    {
      std::array<int,2> n = { { 5, 5 } };
      Dune::FieldVector<double,2> h = { 1.0, 2.0 };
      vtkCheck<2>(n,h);
    }
    {
      std::array<int,3> n = { { 5, 5, 5 } };
      Dune::FieldVector<double,3> h = { 1.0, 2.0, 3.0 };
      vtkCheck<3>(n,h);
    }

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
