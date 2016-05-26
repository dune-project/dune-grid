// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <memory>
#include <vector>
#include <unistd.h>

// dune headers
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

std::string VTKDataMode(Dune::VTK::DataMode dm)
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
  double time_;
public:
  VTKVectorFunction() : time_(0) {}
  void setTime(double time) {
    time_ = time;
  }
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
    return comp*0.1*sin(time_*2.*M_PI);
  }

  // get name
  virtual std::string name () const
  {
    char _name[256];
    snprintf(_name, 256, "vector-%iD", ncomps());
    return std::string(_name);
  }
};

template< class GridView >
void doWrite( const GridView &gridView, Dune::VTK::DataMode dm )
{
  enum { dim = GridView :: dimension };

  const typename GridView :: IndexSet &is = gridView.indexSet();
  std::vector<int> vertexdata(is.size(dim),dim);
  std::vector<int> celldata(is.size(0),0);

  std::stringstream name;
  name << "vtktest-" << dim << "D-" << VTKDataMode(dm);

  auto vtkWriter = std::make_shared<Dune::VTKWriter<GridView> >(gridView, dm);
  Dune :: VTKSequenceWriter< GridView > vtk( vtkWriter, name.str(), ".", "" );

  vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");
  auto vectordata = std::make_shared<VTKVectorFunction<GridView> >();
  vtk.addVertexData(vectordata);
  double time = 0;
  while (time<1) {
    vectordata->setTime(time);
    vtk.write(time);
    time += 0.1;
  }
}

template<int dim>
void vtkCheck(const std::array<int,dim>& n,
              const Dune::FieldVector<double,dim>& h)
{
  std::cout << std::endl << "vtkSequenceCheck dim=" << dim << std::endl << std::endl;
  Dune::YaspGrid<dim> g(h, n);
  g.globalRefine(1);

  doWrite( g.leafGridView(), Dune::VTK::conforming );
  doWrite( g.leafGridView(), Dune::VTK::nonconforming );
  doWrite( g.levelGridView( 0 ), Dune::VTK::conforming );
  doWrite( g.levelGridView( 0 ), Dune::VTK::nonconforming );
  doWrite( g.levelGridView( g.maxLevel() ), Dune::VTK::conforming );
  doWrite( g.levelGridView( g.maxLevel() ), Dune::VTK::nonconforming );
}

int main(int argc, char **argv)
{
  try {

    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if(mpiHelper.rank() == 0)
      std::cout << "subsamplingvtktest: MPI_Comm_size == " << mpiHelper.size()
                << std::endl;

    {
      std::array<int,1> n = { { 5 } };
      Dune::FieldVector<double,1> h = { 1.0 };
      vtkCheck<1>(n,h);
    }
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
