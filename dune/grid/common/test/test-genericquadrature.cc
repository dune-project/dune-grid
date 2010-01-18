// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <limits>
#include <iostream>

#include <config.h>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/quadraturerules/gaussquadrature.hh>

bool success = true;

/*
   This is a simple accuracy test on the reference element. It integrates
   x^p and y^p with the quadrature rule of order p, which should give
   an exact result.
 */

/*
   Exact (analytical) solution on different reference elements.
 */

template <class ctype, int dim>
ctype analyticSolution (Dune::GeometryType t, int p, int x) {
  using Dune::GeometryType;
  ctype exact=0;
  switch (t.basicType())
  {
  case GeometryType::cube :
    exact=1.0/(p+1);
    break;
  case GeometryType::simplex :
    /* 1/(prod(k=1..dim,(p+k)) */
    exact = 1.0;
    for (int k=1; k<=dim; k++) exact*=(p+k);
    exact = 1.0/exact;
    break;
  case GeometryType::prism :
    switch(x) {
    case 0 :
      exact=1.0/((p+2)*(p+1));
      break;
    case 1 :
      exact=1.0/((p+2)*(p+1));
      break;
    case 2 :
      exact=1.0/(2*(p+1));
      break;
    };
    break;
  case GeometryType::pyramid :
    switch(x) {
    case 0 :
    case 1 :
      exact=1.0/((p+3)*(p+1));
      break;
    case 2 :
      exact=2.0/((p+1)*(p+2)*(p+3));
      break;
    };
    break;
  default :
    DUNE_THROW(Dune::NotImplemented, __func__ << " for " << t);
  };
  return exact;
};

template <class ctype, int dim>
void checkQuadrature(Dune::GeometryType t, int maxorder)
{
  using namespace Dune;
  typedef GenericGeometry::GaussQuadratureProvider<dim,double> QuadratureProvider;
  typedef typename QuadratureProvider::Object Quadrature;

  for (int p=0; p<=maxorder; ++p)
  {
    const Quadrature &qr = *QuadratureProvider::create(t,p);
    FieldVector<ctype,dim> integral(0);
    for (typename Quadrature::Iterator qp=qr.begin(); qp!=qr.end(); ++qp)
    {
      // pos of integration point
      FieldVector<ctype,dim> const& x = qp->point();
      ctype weight = qp->weight();

      for (int d=0; d<dim; d++)
      {
        integral[d] += weight*std::pow(x[d],p);
      }
    }

    ctype maxRelativeError = 0;
    int dir = -1;
    for (int d=0; d<dim; d++)
    {
      ctype exact = analyticSolution<ctype,dim>(t,p,d);
      ctype relativeError = std::abs(integral[d]-exact) /
                            (std::abs(integral[d])+std::abs(exact));
      if (relativeError > maxRelativeError)
      {
        maxRelativeError = relativeError;
        dir = d;
      }
    }
    ctype epsilon = std::pow(2.0,p)*p*std::numeric_limits<double>::epsilon();
    if (p==0)
      epsilon = 2.0*std::numeric_limits<double>::epsilon();
    if (maxRelativeError > epsilon) {
      std::cerr << "Error: Quadrature for " << t << " and order=" << p << " failed" << std::endl;
      for (int d=0; d<dim; d++)
      {
        ctype exact = analyticSolution<ctype,dim>(t,p,d);
        ctype relativeError = std::abs(integral[d]-exact) /
                              (std::abs(integral[d])+std::abs(exact));
        std::cerr << "       relative error " << relativeError << " in direction " << d << " (exact = " << exact << " numerical = " << integral[d] << ")" << std::endl;
      }
      success = false;
    }
  }
}

template<class ctype, int dim>
void checkWeights(Dune::GeometryType t, int maxorder)
{
  typedef Dune::GenericGeometry::GaussQuadratureProvider<dim,double> QuadratureProvider;
  typedef typename QuadratureProvider::Object Quadrature;
  typedef typename Quadrature::Iterator QuadIterator;
  for (int p=0; p<=maxorder; ++p)
  {
    double volume = 0;
    // Quadratures
    const Quadrature &quad = *QuadratureProvider::create(t,p);
    QuadIterator qp = quad.begin();
    QuadIterator qend = quad.end();
    for (; qp!=qend; ++qp)
    {
      volume += qp->weight();
    }
    if (std::abs(volume -
                 Dune::GenericReferenceElements<ctype, dim>::general(t).volume())
        > 4*dim*(p ? p : 1)*std::numeric_limits<double>::epsilon())
    {
      std::cerr << "Error: Quadrature for " << t << " and order=" << p
                << " does not sum to volume of RefElem" << std::endl;
      std::cerr << "\tSums to " << volume << "( RefElem.volume() = "
                << Dune::GenericReferenceElements<ctype, dim>::general(t).volume()
                << ")" << "(difference " << volume -
      Dune::GenericReferenceElements<ctype, dim>::general(t).volume()
                << ")" << std::endl;
      success = false;
    }
  }
}

int main ()
{
  try {
    Dune::GeometryType cube0d(Dune::GeometryType::cube,0);
    Dune::GeometryType cube1d(Dune::GeometryType::cube,1);
    Dune::GeometryType cube2d(Dune::GeometryType::cube,2);
    Dune::GeometryType cube3d(Dune::GeometryType::cube,3);
    Dune::GeometryType cube4d(Dune::GeometryType::cube,4);

    Dune::GeometryType simplex0d(Dune::GeometryType::simplex,0);
    Dune::GeometryType simplex1d(Dune::GeometryType::simplex,1);
    Dune::GeometryType simplex2d(Dune::GeometryType::simplex,2);
    Dune::GeometryType simplex3d(Dune::GeometryType::simplex,3);
    Dune::GeometryType simplex4d(Dune::GeometryType::simplex,4);

    Dune::GeometryType prism3d(Dune::GeometryType::prism,3);
    Dune::GeometryType pyramid3d(Dune::GeometryType::pyramid,3);

    checkWeights<double, 0>(cube0d,8);
    checkWeights<double, 1>(cube1d,8);
    checkWeights<double, 2>(cube2d,8);
    checkWeights<double, 3>(cube3d,8);
    checkWeights<double, 4>(cube4d,8);

    checkWeights<double, 0>(simplex0d,8);
    checkWeights<double, 1>(simplex1d,8);
    checkWeights<double, 2>(simplex2d,8);
    checkWeights<double, 3>(simplex3d,8);

    checkWeights<double, 3>(prism3d,8);
    checkWeights<double, 3>(pyramid3d,8);

    checkQuadrature<double, 0>(cube0d,8);
    checkQuadrature<double, 1>(cube1d,8);
    checkQuadrature<double, 2>(cube2d,8);
    checkQuadrature<double, 3>(cube3d,8);
    checkQuadrature<double, 4>(cube4d,8);

    checkQuadrature<double, 0>(simplex0d,8);
    checkQuadrature<double, 1>(simplex1d,8);
    checkQuadrature<double, 2>(simplex2d,8);
    checkQuadrature<double, 3>(simplex3d,8);
    checkQuadrature<double, 4>(simplex4d,8);

    checkQuadrature<double, 3>(prism3d,8);
    checkQuadrature<double, 3>(pyramid3d,8);

  }
  catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 1;
  }

  return success ? 0 : 1;
}
