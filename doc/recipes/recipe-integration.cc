// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file COPYING in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
/**
 * \page recipe-integration Integrate a function on a grid
 *
 * A fundamental task in any PDE solver is integration of a function over a domain given by a grid:
 * \f[
 *  I = \int_\Omega u(x)\,dx = \sum_{e\in E^0} \int_e u(x)\,dx
 * \f]
 * Here we demonstrate how this can be done for a user-defined function in global coordinates.
 *
 * As always we need to set up a grid first. Here we use the \ref Dune::YaspGrid class in
 * four dimensions in order to demonstrate that even dimension larger than three
 * can be used.
 * \snippet recipe-integration.cc set up grid
 *
 * The Dune grid interface employs the class templates Dune::FieldVector and Dune::FieldMatrix to
 * for small vectors and matrices of compile-time known
 * size in many places, for example to represent positions and Jacobians of certain maps.
 * Here is a small tutorial in the usage of these classes:
 * \snippet recipe-integration.cc small vectors and matrices
 *
 * Next we define the function to integrate as a generic lambda function. We expect
 * the argument x to be an instance of type Dune::FieldVector:
 * \snippet recipe-integration.cc a function to integrate
 *
 * First we demonstrate how to compute an approximation to the
 * integral by using the simple midpoint rule. This means traversing all elements,
 * evalute the function at the barycenter of the element, multiplying with the
 * volume of the element and summing up:
 * \snippet recipe-integration.cc integration with midpoint rule
 *
 * A more accurate approximation of the integral for sufficiently smooth functions
 * can be computed by quadrature rules. These are available in Dune for all the different
 * element types and various integration orders. This simply requires an additional
 * loop over the quadrature points within an element:
 * \snippet recipe-integration.cc integration with quadrature rule
 *
 * For integrating the divergence of a vector field \f$f\f$ we might use the following formula:
 * \f[
 *  I = \int_\Omega \nabla\cdot f(x)\,dx = \int_{\partial \Omega} f(x)\cdot n(x) \,ds = \sum_{e\in E^0} \int_{\partial e\cap \partial \Omega} f(x)\cdot n(x) \,ds.
 * \f]
 * This is implemented by the following snippet illustrating the use of intersections:
 * \snippet recipe-integration.cc integrating a flux
 *
 * Full example code: @ref recipe-integration.cc
 * \example recipe-integration.cc
 * See explanation at @ref recipe-integration
 */


// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include <dune/grid/yaspgrid.hh>


int main(int argc, char** argv)
{
  // Maybe initialize Mpi
  [[maybe_unused]] Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // [set up grid]
  const int dim = 4;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> len; for (auto& l:len) l=1.0;
  std::array<int,dim> cells; for (auto& c : cells) c=5;
  Grid grid(len,cells);
  //! [set up grid]

  // [small vectors and matrices]
  Dune::FieldVector<double,4> x({1,2,3,4}); // make a vector
  auto y(x); // copy constructor
  y *= 1.0/3.0; // scaling
  [[maybe_unused]] auto s = x*y; // scalar product
  [[maybe_unused]] auto norm = x.two_norm(); // Euclidean norm
  Dune::FieldMatrix<double,4,4> A({{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}); // make a matrix
  A.mv(x,y); // matvec: y = Ax
  A.usmv(0.5,x,y); // axpy: y += 0.5*Ax
  //! [small vectors and matrices]

  // [a function to integrate]
  auto u = [](const auto& x){return std::exp(x.two_norm());};
  //! [a function to integrate]

  // [integration with midpoint rule]
  double integral=0.0;
  auto gv = grid.leafGridView(); // extract the grid view
  for (const auto& e : elements(gv))
    integral += u(e.geometry().center())*e.geometry().volume();
  std::cout << "integral = " << integral << std::endl;
  //! [integration with midpoint rule]

  // [integration with quadrature rule]
  double integral2 = 0.0;
  using QR = Dune::QuadratureRules<Grid::ctype,dim>;
  for (const auto& e : elements(gv))
    {
      auto geo = e.geometry();
      auto quadrature = QR::rule(geo.type(),5);
      for (const auto& qp : quadrature)
        integral2 += u(geo.global(qp.position()))
          *geo.integrationElement(qp.position())*qp.weight();
    }
  std::cout << "integral2 = " << integral2 << std::endl;
  //! [integration with quadrature rule]

  // [integrating a flux]
  auto f = [](const auto& x){return x;};
  double divergence=0.0;
  for (const auto& i : elements(gv)) {
    for (const auto& I : intersections(gv,i))
      if (!I.neighbor())
        {
          auto geoI = I.geometry();
          divergence += f(geoI.center())*I.centerUnitOuterNormal()*geoI.volume();
        }
  }
  std::cout << "divergence = " << divergence << std::endl;
  //! [integrating a flux]
}
