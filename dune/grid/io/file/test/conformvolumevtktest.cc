// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/volumewriter.hh>
#include <dune/grid/yaspgrid.hh>

template< class GridView >
class VTKVectorFunction
  : public Dune::VTK::ConformingVolumeWriter<GridView>::VTKFunction
{
  // extract types
  constexpr static int n = GridView :: dimension;
  constexpr static int w = GridView :: dimensionworld;
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
  virtual double evaluate (int comp, [[maybe_unused]] const Entity& e,
                           [[maybe_unused]] const Dune::FieldVector<DT,n>& xi) const
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
void doWrite( const GridView &gridView )
{
  constexpr static int dim = GridView :: dimension;

  const typename GridView :: IndexSet &is = gridView.indexSet();
  std::vector<int> vertexdata(is.size(dim),dim);
  std::vector<int> celldata(is.size(0),0);

  Dune::VTK::ConformingVolumeWriter< GridView > vtk( gridView );
  vtk.addVertexData(vertexdata,"vertexData");
  vtk.addCellData(celldata,"cellData");

  vtk.addVertexData(std::make_shared< VTKVectorFunction<GridView> >("vertex"));
  vtk.addCellData(std::make_shared< VTKVectorFunction<GridView> >("cell"));

  char name[256];
  snprintf(name,256,"conformvolumevtktest-%iD-ascii", dim);
  vtk.write(name, Dune::VTK::ascii);

  snprintf(name,256,"conformvolumevtktest-%iD-base64", dim);
  vtk.write(name, Dune::VTK::base64);

  snprintf(name,256,"conformvolumevtktest-%iD-appendedraw", dim);
  vtk.write(name, Dune::VTK::appendedraw);

  snprintf(name,256,"conformvolumevtktest-%iD-appendedbase64", dim);
  vtk.write(name, Dune::VTK::appendedbase64);
}

template<unsigned int dim>
void vtkCheck(Dune::FieldVector<double,dim>& length, std::array<int,dim> elements)
{
  std::cout << "vtkCheck dim=" << dim << std::endl;
  Dune::YaspGrid<dim> grid(length, elements);
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
      const unsigned int dim = 1;
      Dune::FieldVector<double,dim> length(1.0);
      std::array<int,dim> elements;
      std::fill(elements.begin(), elements.end(), 5);
      vtkCheck<dim>(length, elements);
    }

    {
      const unsigned int dim = 2;
      Dune::FieldVector<double,dim> length(1.0);
      length[1] = 2.0;
      std::array<int,dim> elements;
      std::fill(elements.begin(), elements.end(), 5);
      vtkCheck<dim>(length, elements);
    }

    {
      const unsigned int dim = 3;
      Dune::FieldVector<double,dim> length(1.0);
      length[1] = 2.0;
      length[2] = 3.0;
      std::array<int,dim> elements;
      std::fill(elements.begin(), elements.end(), 5);
      vtkCheck<dim>(length, elements);
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
