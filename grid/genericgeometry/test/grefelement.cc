// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/genericreferenceelements.hh>

using namespace Dune;

int main ( int argc, char **argv )
{
  const GeometryType type( GeometryType::simplex, 3 );
  const GenericReferenceElement< double, 3 > &refElement
    = GenericReferenceElements< double, 3 >::general( type );

  std::cout << refElement.size( 0 ) << std::endl;

  return 0;
}
