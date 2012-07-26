// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id:$

#if HAVE_CONFIG_H
#include "config.h" // autoconf defines, needed by the dune headers
#endif

#include <algorithm>
#include <iostream>
#include <ostream>
#include <vector>

#include <unistd.h>

#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

const char* VTKDataMode(Dune::VTK::DataMode dm)
{
  switch(dm)
  {
  case Dune::VTK::conforming :
    return "conforming";
  case Dune::VTK::nonconforming :
    return "nonconforming";
  }
  return "";
}

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
  virtual double evaluate (int comp, const Entity& e,
                           const Dune::FieldVector<DT,n>& xi) const
  {
    return comp*0.1;
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
void doWrite( const GridView &gridView, Dune :: VTK :: DataMode dm )
{
  enum { dim = GridView :: dimension };

  const typename GridView :: IndexSet &is = gridView.indexSet();
  std::vector<int> vertexdata(is.size(dim),dim);
  std::vector<int> celldata(is.size(0),0);

  Dune :: VTKWriter< GridView > vtk( gridView, dm );
  vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");

  vtk.addVertexData(new VTKVectorFunction<GridView>("vertex"));
  vtk.addCellData(new VTKVectorFunction<GridView>("cell"));

  char name[256];
  snprintf(name,256,"vtktest-%iD-%s-ascii", dim, VTKDataMode(dm));
  vtk.write(name);

  snprintf(name,256,"vtktest-%iD-%s-base64", dim, VTKDataMode(dm));
  vtk.write(name, Dune::VTK::base64);

  snprintf(name,256,"vtktest-%iD-%s-appendedraw", dim, VTKDataMode(dm));
  vtk.write(name, Dune::VTK::appendedraw);

  snprintf(name,256,"vtktest-%iD-%s-appendedbase64", dim, VTKDataMode(dm));
  vtk.write(name, Dune::VTK::appendedbase64);
}

template<int dim>
void vtkCheck(const Dune::MPIHelper &mpiHelper, int* n, double* h)
{
  const Dune :: PartitionIteratorType VTK_Partition
    = Dune :: InteriorBorder_Partition;
  if(mpiHelper.rank() == 0)
    std::cout << std::endl << "vtkCheck dim=" << dim << std::endl << std::endl;

  typedef Dune::YaspGrid<dim> Grid;
  Dune::FieldVector<typename Grid::ctype, dim> L(0);
  std::copy(h, h+dim, L.begin());
  Dune::FieldVector<int, dim> s(0);
  std::copy(n, n+dim, s.begin());

  Dune::YaspGrid<dim> g(mpiHelper.getCommunicator(), L, s,
                        Dune::FieldVector<bool, dim>(false), 0);
  g.globalRefine(1);

  doWrite( g.template leafView< VTK_Partition >(), Dune::VTK::conforming );
  doWrite( g.template leafView< VTK_Partition >(), Dune::VTK::nonconforming );
  doWrite( g.template levelView< VTK_Partition >( 0 ),
           Dune::VTK::conforming );
  doWrite( g.template levelView< VTK_Partition >( 0 ),
           Dune::VTK::nonconforming );
  doWrite( g.template levelView< VTK_Partition >( g.maxLevel() ),
           Dune::VTK::conforming );
  doWrite( g.template levelView< VTK_Partition >( g.maxLevel() ),
           Dune::VTK::nonconforming );
}

int main(int argc, char **argv)
{
  try {
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if(mpiHelper.rank() == 0)
      std::cout << "vtktest: MPI_Comm_size == " << mpiHelper.size()
                << std::endl;

    int n[] = { 5, 5, 5, 5 };
    double h[] = { 1.0, 2.0, 3.0, 4.0 };

    vtkCheck<1>(mpiHelper,n,h);
    vtkCheck<2>(mpiHelper,n,h);
    vtkCheck<3>(mpiHelper,n,h);

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
