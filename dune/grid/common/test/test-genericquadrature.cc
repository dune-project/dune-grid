// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

#include <limits>
#include <iostream>

#include <config.h>

#include <dune/common/misc.hh>

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
ctype analyticalSolution (Dune::GeometryType t, int p, int direction )
{
  using Dune::GeometryType;
  ctype exact = 0;
  switch (t.basicType())
  {
  case GeometryType::cube:
    exact=1.0/(p+1);
    break;

  case GeometryType::simplex:
    /* 1/(prod(k=1..dim,(p+k)) */
    exact = ctype( 1 );
    for( int k = 1; k <= dim; ++k )
      exact *= p+k;
    exact = ctype( 1 ) / exact;
    break;

  case GeometryType::prism:
  {
    const int pdim = (dim > 0 ? dim-1 : 0);
    if( direction < dim-1 )
    {
      GeometryType nt( GeometryType::simplex, dim-1 );
      if( dim > 0 )
        exact = analyticalSolution< ctype, pdim >( nt, p, direction );
      else
        exact = ctype( 1 );
    }
    else
      exact = ctype( 1 ) / ctype( Dune::Factorial< pdim >::factorial * (p+1));
    break;
  }

  case GeometryType::pyramid:
    switch( direction )
    {
    case 0:
    case 1:
      exact=1.0/((p+3)*(p+1));
      break;
    case 2:
      exact=2.0/((p+1)*(p+2)*(p+3));
      break;
    };
    break;
  default:
    DUNE_THROW(Dune::NotImplemented, __func__ << " for " << t);
  };
  return exact;
};

template<class Quadrature>
void checkQuadrature(const Quadrature &quad)
{
  using namespace Dune;
  typedef typename Quadrature::Field ctype;
  const unsigned int dim = Quadrature::dimension;
  const unsigned int p = quad.order();
  const Dune::GeometryType& t = quad.type();
  FieldVector<ctype,dim> integral(0);
  for (typename Quadrature::Iterator qp=quad.begin(); qp!=quad.end(); ++qp)
  {
    // pos of integration point
    const FieldVector< ctype, dim > &x = qp->position();
    const ctype weight = qp->weight();

    for (unsigned int d=0; d<dim; d++)
      integral[d] += weight*std::pow(x[d],double(p);
  }
  
  ctype maxRelativeError = 0;
  int dir = -1;
  for( unsigned int d=0; d<dim; d++ )
  {
    ctype exact = analyticalSolution<ctype,dim>(t,p,d);
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
    for (unsigned int d=0; d<dim; d++)
    {
      ctype exact = analyticalSolution<ctype,dim>(t,p,d);
      ctype relativeError = std::abs(integral[d]-exact) /
        (std::abs(integral[d])+std::abs(exact));
      std::cerr << "       relative error " << relativeError << " in direction " << d << " (exact = " << exact << " numerical = " << integral[d] << ")" << std::endl;
    }
    success = false;
  }
}  

template<class Quadrature>
void checkWeights(const Quadrature &quad)
{
  typedef typename Quadrature::Field ctype;
  const unsigned int dim = Quadrature::dimension;
  const unsigned int p = quad.order();
  const Dune::GeometryType& t = quad.type();
  typedef typename Quadrature::Iterator QuadIterator;
  double volume = 0;
  QuadIterator qp = quad.begin();
  QuadIterator qend = quad.end();
  for (; qp!=qend; ++qp)
  {
    volume += qp->weight();
  }
  if (std::abs(volume -
               Dune::GenericReferenceElements<ctype, dim>::general(t).volume())
      > 4*dim*(p?p:1)*std::numeric_limits<double>::epsilon())
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

template<class CF, int dim>
void check( const Dune::GeometryType::BasicType &btype, unsigned int maxOrder )
{
  typedef Dune::GenericGeometry::GaussPoints<CF> OneDPoints;
  typedef Dune::GenericGeometry::GenericQuadratureFactory<dim,double,OneDPoints> QuadratureProvider;
  typedef typename QuadratureProvider::Object Quadrature;
  for (unsigned int p=0;p<=maxOrder; ++p)
  {
    const Quadrature &quad = *QuadratureProvider::create(Dune::GeometryType(btype,dim),p);
    checkWeights(quad);
    checkQuadrature(quad);
    QuadratureProvider::release(&quad);
  }
  if (dim>0 && (dim>2 || 
               btype==Dune::GeometryType::cube ||
               btype==Dune::GeometryType::simplex) )
    check<CF,(dim==0)?0:dim-1>(btype,maxOrder);
}

int main ()
{
  try {
    check<double,4>(Dune::GeometryType::cube,30);
    check<double,4>(Dune::GeometryType::simplex,55);
    check<double,4>(Dune::GeometryType::prism,55);
    check<double,3>(Dune::GeometryType::pyramid,55);
 }
  catch( const Dune::Exception &e )
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 1;
  }
  
  return success ? 0:1;
}
