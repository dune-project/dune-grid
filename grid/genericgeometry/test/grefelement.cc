// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/genericreferenceelements.hh>

using namespace Dune;


template<class T, int dim>
void testWithGeometryType(Dune::GeometryType gt)
{
  typedef Dune::GenericGeometry::MapNumberingProvider<dim> Numbering;
  const unsigned int gtId = Dune::GenericGeometry::topologyId(gt);


  const ReferenceElement<T, dim> &ref = ReferenceElements<T,dim>::general(gt);
  const GenericReferenceElement<T, dim> &gref = GenericReferenceElements<T,dim>::general(gt);

  std::cout << "********************* testing for " << gt << std::endl;

  // test volume
  if (ref.volume() != gref.volume())
    std::cout << "ref.volume() != gref.volume()" << std::endl;

  for(int codim=0; codim<=dim; ++codim)
  {
    std::cout << "codim " << codim << std::endl;

    //test size
    if (ref.size(codim) != gref.size(codim))
      std::cout << "ref.size(codim) != gref.size(codim)" << std::endl;

    for(int i=0; i<ref.size(codim); ++i)
    {
      int gi = Numbering::dune2generic(gtId, i, codim);

      std::cout << "i     " << i << std::endl;
      std::cout << "gi    " << gi << std::endl;

      //test position
      if (ref.position(i, codim) != gref.position(gi, codim))
        std::cout << "ref.position(i, codim) != gref.position(gi, codim)" << std::endl;

      //test codim
      if (ref.type(i, codim) != gref.type(gi, codim))
        std::cout << "ref.type(i, codim) != gref.type(gi, codim)" << std::endl;

      for(int c=codim; c<=dim; ++c)
      {
        if (ref.size(i, codim, c) != gref.size(gi, codim, c))
          std::cout << "ref.size(i, codim, c) != gref.size(gi, codim, c)" << std::endl;
      }
    }
  }
}


int main ( int argc, char **argv )
{
  const GeometryType type( GeometryType::simplex, 3 );
  const GenericReferenceElement< double, 3 > &refElement
    = GenericReferenceElements< double, 3 >::general( type );

  std::cout << refElement.size( 0 ) << std::endl;

  typedef Dune::GeometryType GT;

  testWithGeometryType<double,1>(GT(GT::simplex, 1));
  testWithGeometryType<double,2>(GT(GT::simplex, 2));
  testWithGeometryType<double,3>(GT(GT::simplex, 3));

  testWithGeometryType<double,1>(GT(GT::cube, 1));
  testWithGeometryType<double,2>(GT(GT::cube, 2));
  testWithGeometryType<double,3>(GT(GT::cube, 3));

  testWithGeometryType<double,3>(GT(GT::prism, 3));

  testWithGeometryType<double,3>(GT(GT::pyramid, 3));


  return 0;
}
