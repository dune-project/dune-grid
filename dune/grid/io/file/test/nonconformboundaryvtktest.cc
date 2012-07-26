// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/boundarywriter.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/skeletonfunction.hh>
#include <dune/grid/sgrid.hh>

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
  enum { dim = GridView :: dimension };

  Dune::VTK::NonConformingBoundaryWriter< GridView > vtk( gridView );

  Dune::shared_ptr<ScalarFunction<GridView> > scalarFunc
    (new ScalarFunction<GridView>);
  vtk.addCellData(scalarFunc, "cellScalar");
  vtk.addPointData(scalarFunc, "pointScalar");

  Dune::shared_ptr<VectorFunction<GridView> > vectorFunc
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
void vtkCheck(int* n, double* h)
{
  const Dune :: PartitionIteratorType VTK_Partition = Dune :: InteriorBorder_Partition;
  std::cout << std::endl << "vtkCheck dim=" << dim << std::endl << std::endl;
  Dune::SGrid<dim,dim> g(n, h);
  g.globalRefine(1);

  doWrite( g.template leafView< VTK_Partition >() );
  doWrite( g.template levelView< VTK_Partition >( 0 ) );
  doWrite( g.template levelView< VTK_Partition >( g.maxLevel() ) );
}

int main(int argc, char **argv)
{
  try {

    int n[] = { 5, 5, 5, 5 };
    double h[] = { 1.0, 2.0, 3.0, 4.0 };

    vtkCheck<2>(n,h);
    vtkCheck<3>(n,h);

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (std::exception &e) {
    std::cerr << e.what << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
