// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h" // autoconf defines, needed by the dune headers

// dune headers
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

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
  virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
  {
    Dune::FieldVector<DT,n> global = e.geometry().global( xi );
    return global.two_norm2();
  }

  // get name
  virtual std::string name () const
  {
    return type_ + "-vector-" + std::to_string(ncomps()) + "D";
  }


};

template< class GridView >
void doWrite( const GridView &gridView, bool coerceToSimplex)
{
  enum { dim = GridView :: dimension };

  Dune :: SubsamplingVTKWriter< GridView > vtk( gridView, 1, coerceToSimplex);
  // disabled due to FS#676:
  // vtk.addVertexData(vertexdata,"vertexData");
  // vtk.addCellData(celldata,"cellData");

  vtk.addVertexData(std::make_shared< VTKVectorFunction<GridView> >("vertex"));
  vtk.addCellData(std::make_shared< VTKVectorFunction<GridView> >("cell"));

  std::string name = "subsamplingvtktest-" + std::to_string(dim) + "D-" + (coerceToSimplex ? "simplex" : "natural") + "-ascii";
  vtk.write(name);

  name = "subsamplingvtktest-" + std::to_string(dim) + "D-" + (coerceToSimplex ? "simplex" : "natural") + "-appendedraw";
  vtk.write(name, Dune::VTK::appendedraw);
}

template<int dim>
void vtkCheck(const Dune::array<int, dim>& elements,
              const Dune::FieldVector<double, dim>& upperRight)
{
  Dune::YaspGrid<dim> g(upperRight, elements);

  if(g.comm().rank() == 0)
    std::cout << std::endl
              << "subsamplingVTKCheck dim=" << dim
              << std::endl
              << std::endl;

  g.globalRefine(1);

  doWrite( g.leafGridView(), false);
  doWrite( g.levelGridView( 0 ), false);
  doWrite( g.levelGridView( g.maxLevel() ), false);

  doWrite( g.leafGridView(), true);
  doWrite( g.levelGridView( 0 ), true);
  doWrite( g.levelGridView( g.maxLevel() ), true);
}

int main(int argc, char **argv)
{
  try {

    Dune::MPIHelper::instance(argc, argv);

    vtkCheck<1>({5}, {1.0});
    vtkCheck<2>({5,5}, {1.0, 2.0});
    vtkCheck<3>({5,5,5}, {1.0, 2.0, 3.0});

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
