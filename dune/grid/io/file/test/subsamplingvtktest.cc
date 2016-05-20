// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h" // autoconf defines, needed by the dune headers
#endif

#include <algorithm>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/test/checkvtkfile.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

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
    Dune::FieldVector<DT,n> global = e.geometry().global( xi );
    return global.two_norm2();
  }

  // get name
  virtual std::string name () const
  {
    std::ostringstream os;
    os << type_ << "-vector-" << ncomps() << "D";
    return os.str();
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
int doWrite( const GridView &gridView, bool coerceToSimplex)
{
  enum { dim = GridView :: dimension };

  Dune :: SubsamplingVTKWriter< GridView > vtk( gridView, 1, coerceToSimplex);

  // disabled due to FS#676:
  // const typename GridView :: IndexSet &is = gridView.indexSet();
  // std::vector<int> vertexdata(is.size(dim),dim);
  // std::vector<int> celldata(is.size(0),0);

  // vtk.addVertexData(vertexdata,"vertexData");
  // vtk.addCellData(celldata,"cellData");

  vtk.addVertexData(std::make_shared< VTKVectorFunction<GridView> >("vertex"));
  vtk.addCellData(std::make_shared< VTKVectorFunction<GridView> >("cell"));

  int result = 0;
  std::string name;
  std::ostringstream prefix;
  prefix << "subsamplingvtktest-" << dim << "D-"
         << (coerceToSimplex ? "simplex" : "natural");
  int rank = gridView.comm().rank();

  name = vtk.write(prefix.str() + "-ascii");
  if(rank == 0) acc(result, checkVTKFile(name));

  name = vtk.write(prefix.str() + "-appendedraw", Dune::VTK::appendedraw);
  if(rank == 0) acc(result, checkVTKFile(name));

  return result;
}

template<int dim>
int vtkCheck(const std::array<int, dim>& elements,
              const Dune::FieldVector<double, dim>& upperRight)
{
  Dune::YaspGrid<dim> g(upperRight, elements);

  if(g.comm().rank() == 0)
    std::cout << std::endl
              << "subsamplingVTKCheck dim=" << dim << std::endl
              << std::endl;

  g.globalRefine(1);

  int result = 0;

  acc(result, doWrite( g.leafGridView(), false));
  acc(result, doWrite( g.levelGridView( 0 ), false));
  acc(result, doWrite( g.levelGridView( g.maxLevel() ), false));

  acc(result, doWrite( g.leafGridView(), true));
  acc(result, doWrite( g.levelGridView( 0 ), true));
  acc(result, doWrite( g.levelGridView( g.maxLevel() ), true));

  return result;
}

int main(int argc, char **argv)
{
  try {

    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if(mpiHelper.rank() == 0)
      std::cout << "subsamplingvtktest: MPI_Comm_size == " << mpiHelper.size()
                << std::endl;

    int result = 0; // pass by default
    using Dune::make_array;

    acc(result, vtkCheck<1>(make_array(5), {1.0}));
    acc(result, vtkCheck<2>(make_array(5,5), {1.0, 2.0}));
    acc(result, vtkCheck<3>(make_array(5,5,5), {1.0, 2.0, 3.0}));

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
