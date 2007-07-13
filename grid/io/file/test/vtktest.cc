// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id:$

#include "config.h" // autoconf defines, needed by the dune headers

// dune headers
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <vector>
#include <unistd.h>

const char* VTKDataMode(Dune::VTKOptions::DataMode dm)
{
  switch(dm)
  {
  case Dune::VTKOptions::conforming :
    return "conforming";
  case Dune::VTKOptions::nonconforming :
    return "nonconforming";
  }
  return "";
}

template<class G, class IS>
class VTKVectorFuction : public Dune::VTKWriter<G,IS>::VTKFunction
{
  // extract types
  enum {n=G::dimension};
  enum {w=G::dimensionworld};
  typedef typename G::ctype DT;
  typedef typename G::Traits::template Codim<0>::Entity Entity;
public:
  //! return number of components
  virtual int ncomps () const { return n; };

  //! evaluate single component comp in the entity e at local coordinates xi
  /*! Evaluate the function in an entity at local coordinates.
     @param[in]  comp   number of component to be evaluated
     @param[in]  e      reference to grid entity of codimension 0
     @param[in]  xi     point in local coordinates of the reference element of e
     \return            value of the component
   */
  virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
  {
    return comp*0.1;
  }

  // get name
  virtual std::string name () const
  {
    char _name[256];
    snprintf(_name, 256, "vector-%iD", ncomps());
    return std::string(_name);
  };


};

template<class G, class IS>
void doWrite(G & g, IS & is, Dune::VTKOptions::DataMode dm)
{
  enum { dim = G::dimension };

  Dune::VTKWriter<G, IS> vtk(g,is,dm);
  std::vector<int> vertexdata(is.size(dim),dim);
  std::vector<int> celldata(is.size(0),0);
  vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");

  VTKVectorFuction<G, IS> * vectordata = new VTKVectorFuction<G, IS>;
  vtk.addVertexData(vectordata);

  char name[256];
  snprintf(name,256,"vtktest-%iD-%s-ascii", dim, VTKDataMode(dm));
  vtk.write(name);

  snprintf(name,256,"vtktest-%iD-%s-binary", dim, VTKDataMode(dm));
  vtk.write(name, Dune::VTKOptions::binaryappended);
}

template<int dim>
void vtkCheck(int* n, double* h)
{
  std::cout << std::endl << "vtkCheck dim=" << dim << std::endl << std::endl;
  Dune::SGrid<dim,dim> g(n, h);
  g.globalRefine(2);

  doWrite(g,g.leafIndexSet(),Dune::VTKOptions::conforming);
  doWrite(g,g.leafIndexSet(),Dune::VTKOptions::nonconforming);
  doWrite(g,g.levelIndexSet(0),Dune::VTKOptions::conforming);
  doWrite(g,g.levelIndexSet(0),Dune::VTKOptions::nonconforming);
  doWrite(g,g.levelIndexSet(g.maxLevel()),Dune::VTKOptions::conforming);
  doWrite(g,g.levelIndexSet(g.maxLevel()),Dune::VTKOptions::nonconforming);
}

int main(int argc, char **argv)
{
  try {

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
