// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id:$

#include "config.h" // autoconf defines, needed by the dune headers

// dune headers
#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <algorithm>
#include <vector>
#include <unistd.h>

template< class GridView >
class VTKVectorFunction
  : public Dune :: VTKWriter< GridView > :: VTKFunction
{
  // extract types
  enum { n = GridView :: dimension };
  enum { w = GridView :: dimensionworld };
  typedef typename GridView :: Grid :: ctype DT;
  typedef typename GridView :: template Codim< 0 > :: Entity Entity;
  const char *type;
public:
  /** @brief Make a new VTKVectorFunction
   *
   * @param type_ Type of the function for use in its name (hint: "cell" or
   *              "vertex")
   */
  VTKVectorFunction(const char *type_)
    : type(type_)
  { }

  //! return number of components
  virtual int ncomps () const { return n; }

  //! evaluate single component comp in the entity e at local coordinates xi
  /*! Evaluate the function in an entity at local coordinates.
     @param[in]  comp   number of component to be evaluated
     @param[in]  e      reference to grid entity of codimension 0
     @param[in]  xi     point in local coordinates of the reference element of e
     \return            value of the component
   */
  virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
  {
    Dune::FieldVector<DT,n> global = e.geometry().global( xi );
    return global.two_norm2();
  }

  // get name
  virtual std::string name () const
  {
    char _name[256];
    snprintf(_name, 256, "%s-vector-%iD", type, ncomps());
    return std::string(_name);
  }


};

template< class GridView >
void doWrite( const GridView &gridView, bool coerceToSimplex)
{
  enum { dim = GridView :: dimension };

  const typename GridView :: IndexSet &is = gridView.indexSet();
  std::vector<int> celldata(is.size(0));
  for(std::size_t i = 0; i < celldata.size(); ++i) celldata[i] = i;

  Dune :: SubsamplingVTKWriter< GridView > vtk( gridView, 1, coerceToSimplex);
  // disabled due to FS#676: vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");

  vtk.addVertexData(new VTKVectorFunction<GridView>("vertex"));
  vtk.addCellData(new VTKVectorFunction<GridView>("cell"));

  char name[256];
  snprintf(name,256,"subsamplingvtktest-%iD-%s-ascii",
           dim, (coerceToSimplex ? "simplex" : "natural"));
  vtk.write(name);

  snprintf(name,256,"subsamplingvtktest-%iD-%s-appendedraw",
           dim, (coerceToSimplex ? "simplex" : "natural"));
  vtk.write(name, Dune::VTK::appendedraw);
}

template<int dim>
void vtkCheck(int* n, double* h)
{
  Dune::FieldVector<double, dim> L(0);
  std::copy(h, h+dim, L.begin());
  Dune::FieldVector<int, dim> s(0);
  std::copy(n, n+dim, s.begin());
  Dune::FieldVector<bool, dim> periodic(false);

  Dune::YaspGrid<dim> g(
#if HAVE_MPI
    MPI_COMM_WORLD,
#endif
    L, s, periodic, 0);
  if(g.comm().rank() == 0)
    std::cout << std::endl
              << "subsamplingVTKCheck dim=" << dim
              << std::endl
              << std::endl;

  g.globalRefine(1);

  doWrite( g.leafView(), false);
  doWrite( g.levelView( 0 ), false);
  doWrite( g.levelView( g.maxLevel() ), false);

  doWrite( g.leafView(), true);
  doWrite( g.levelView( 0 ), true);
  doWrite( g.levelView( g.maxLevel() ), true);
}

int main(int argc, char **argv)
{
  try {

    Dune::MPIHelper::instance(argc, argv);

    int n[] = { 5, 5, 5, 5 };
    double h[] = { 1.0, 2.0, 3.0, 4.0 };

    vtkCheck<1>(n,h);
    vtkCheck<2>(n,h);
    vtkCheck<3>(n,h);

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
