// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <limits>
#include <iostream>

#include <config.h>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/common/genericreferenceelements.hh>

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
void checkQuadrature(Dune::GeometryType t)
{
  using namespace Dune;

  for (int p=0; ; ++p)
    try {
      QuadratureRule<ctype,dim> const& qr =
        QuadratureRules<ctype,dim>::rule(t,p);
      FieldVector<ctype,dim> integral(0);
      for (typename QuadratureRule<ctype,dim>::const_iterator
           qp=qr.begin(); qp!=qr.end(); ++qp)
      {
        // pos of integration point
        FieldVector<ctype,dim> const& x = qp->position();
        ctype weight = qp->weight();

        for (int d=0; d<dim; d++)
        {
          integral[d] += weight*std::pow(x[d],double(p));
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
      ctype epsilon = std::pow(2.0,double(p))*p*std::numeric_limits<double>::epsilon();
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
    catch (Dune::QuadratureOrderOutOfRange & e) {
      std::cout << "tested integration for " << t << std::endl;
      break;
    }
}

template<class ctype, int dim>
void checkWeights(Dune::GeometryType t, int p)
{
  double volume = 0;
  // Quadratures
  typedef Dune::QuadratureRule<ctype, dim> Quad;
  typedef typename Quad::iterator QuadIterator;
  const Quad & quad = Dune::QuadratureRules<ctype,dim>::rule(t, p);
  if (quad.type() != t || quad.order() < p) {
    std::cerr << "Error: Type mismatch! Requested Quadrature for " << t
              << " and order=" << p << "." << std::endl
              << "\tGot Quadrature for " << quad.type() << " and order="
              << quad.order() << std::endl;
    success = false;
    return;
  }
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

template<class ctype, int dim>
void checkWeights(Dune::GeometryType t)
{
  int maxorder;
  for (int i=0;; i++)
  {
    try {
      checkWeights<ctype,dim>(t, i);
    }
    catch (Dune::QuadratureOrderOutOfRange & e) {
      maxorder = i-1;
      break;
    }
  }
  for (int i=maxorder+1;; i++)
  {
    try {
      checkWeights<ctype,dim>(t, i);
    }
    catch (Dune::QuadratureOrderOutOfRange & e) {
      if (i > maxorder+1)
      {
        std::cout << "Error: " << t << " allows higher order in the second run." << std::endl;
        std::cout << "       " << maxorder << " in the first run, "
                  << i-1 << " in the second run." << std::endl;
      }
      std::cout << "tested weights for " << t << " up to max order = " << maxorder << std::endl;
      break;
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
    // Dune::GeometryType cube4d(Dune::GeometryType::cube,4);
    // Dune::GeometryType cube5d(Dune::GeometryType::cube,5);
    // Dune::GeometryType cube6d(Dune::GeometryType::cube,6);

    Dune::GeometryType simplex0d(Dune::GeometryType::simplex,0);
    Dune::GeometryType simplex1d(Dune::GeometryType::simplex,1);
    Dune::GeometryType simplex2d(Dune::GeometryType::simplex,2);
    Dune::GeometryType simplex3d(Dune::GeometryType::simplex,3);

    Dune::GeometryType prism3d(Dune::GeometryType::prism,3);
    Dune::GeometryType pyramid3d(Dune::GeometryType::pyramid,3);

    checkWeights<double, 0>(cube0d);
    checkWeights<double, 1>(cube1d);
    checkWeights<double, 2>(cube2d);
    checkWeights<double, 3>(cube3d);
    // checkWeights<double, 4>(cube4d);
    // checkWeights<double, 5>(cube5d);
    // checkWeights<double, 6>(cube6d);

    checkWeights<double, 0>(simplex0d);
    checkWeights<double, 1>(simplex1d);
    checkWeights<double, 2>(simplex2d);
    checkWeights<double, 3>(simplex3d);

    checkWeights<double, 3>(prism3d);
    checkWeights<double, 3>(pyramid3d);

    checkQuadrature<double, 0>(cube0d);
    checkQuadrature<double, 1>(cube1d);
    checkQuadrature<double, 2>(cube2d);
    checkQuadrature<double, 3>(cube3d);
    // checkQuadrature<double, 4>(cube4d);
    // checkQuadrature<double, 5>(cube5d);
    // checkQuadrature<double, 6>(cube6d);

    checkQuadrature<double, 0>(simplex0d);
    checkQuadrature<double, 1>(simplex1d);
    checkQuadrature<double, 2>(simplex2d);
    checkQuadrature<double, 3>(simplex3d);
    checkQuadrature<double, 3>(prism3d);
    checkQuadrature<double, 3>(pyramid3d);
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
