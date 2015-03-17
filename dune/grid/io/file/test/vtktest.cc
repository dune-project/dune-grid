// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h" // autoconf defines, needed by the dune headers
#endif

#include <algorithm>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/test/checkvtkfile.hh>
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
  typedef typename GridView :: Grid :: ctype DT;
  typedef typename GridView :: template Codim< 0 > :: Entity Entity;
  const std::string type_;
public:
  /** @brief Make a new VTKVectorFunction
   *
   * @param type_ Type of the function for use in its name (hint: "cell" or
   *              "vertex")
   */
  VTKVectorFunction(std::string type)
    : type_(type)
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
    return type_ + "-vector-" + std::to_string(ncomps()) + "D";
  }

};

// accumulate exit status
void acc(int &accresult, int result)
{
  if(accresult == 0 || (accresult == 77 && result != 0))
    accresult = result;
}

struct Acc
{
  int operator()(int v1, int v2) const
  {
    acc(v1, v2);
    return v1;
  }
};

template< class GridView >
int doWrite( const GridView &gridView, Dune :: VTK :: DataMode dm )
{
  enum { dim = GridView :: dimension };

  const typename GridView :: IndexSet &is = gridView.indexSet();
  std::vector<int> vertexdata(is.size(dim),dim);
  std::vector<int> celldata(is.size(0),0);

  Dune :: VTKWriter< GridView > vtk( gridView, dm );
  vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");

  vtk.addVertexData(std::make_shared< VTKVectorFunction<GridView> >("vertex"));
  vtk.addCellData(std::make_shared< VTKVectorFunction<GridView> >("cell"));

  int result = 0;
  std::string name;
  std::string prefix =
    "vtktest-" + std::to_string(dim) + "D-" + VTKDataMode(dm);
  int rank = gridView.comm().rank();

  name = vtk.write(prefix + "-ascii");
  if(rank == 0) acc(result, checkVTKFile(name));

  name = vtk.write(prefix + "-base64", Dune::VTK::base64);
  if(rank == 0) acc(result, checkVTKFile(name));

  name = vtk.write(prefix + "-appendedraw", Dune::VTK::appendedraw);
  if(rank == 0) acc(result, checkVTKFile(name));

  name = vtk.write(prefix + "-appendedbase64", Dune::VTK::appendedbase64);
  if(rank == 0) acc(result, checkVTKFile(name));

  return result;
}

template<int dim>
int vtkCheck(const Dune::MPIHelper &mpiHelper, int* n, double* h)
{
  if(mpiHelper.rank() == 0)
    std::cout << std::endl << "vtkCheck dim=" << dim << std::endl << std::endl;

  typedef Dune::YaspGrid<dim> Grid;
  Dune::FieldVector<typename Grid::ctype, dim> L(0);
  std::copy(h, h+dim, L.begin());
  Dune::array<int, dim> s;
  std::copy(n, n+dim, s.begin());

  Dune::YaspGrid<dim> g(L, s, std::bitset<dim>(), 0);

  g.globalRefine(1);

  int result = 0;

  acc(result, doWrite( g.template leafGridView(), Dune::VTK::conforming ));
  acc(result, doWrite( g.template leafGridView(), Dune::VTK::nonconforming ));
  acc(result, doWrite( g.template levelGridView( 0 ),
                       Dune::VTK::conforming ));
  acc(result, doWrite( g.template levelGridView( 0 ),
                       Dune::VTK::nonconforming ));
  acc(result, doWrite( g.template levelGridView( g.maxLevel() ),
                       Dune::VTK::conforming ));
  acc(result, doWrite( g.template levelGridView( g.maxLevel() ),
                       Dune::VTK::nonconforming ));

  return result;
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

    int result = 0; // pass by default

    acc(result, vtkCheck<1>(mpiHelper,n,h));
    acc(result, vtkCheck<2>(mpiHelper,n,h));
    acc(result, vtkCheck<3>(mpiHelper,n,h));

    mpiHelper.getCollectiveCommunication().allreduce<Acc>(&result, 1);

    return result;

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

}
