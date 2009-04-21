// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune
{

  template< typename ctype, int dim >
  GenericReferenceElementContainer< ctype, dim > GenericReferenceElements< ctype, dim>::general;

  // we use an empty namespace to initialize the singleton,
  // so that this code is hidden
  namespace
  {

    template< class ctype >
    struct InitReferenceElements
    {
      template< int dim >
      struct Dimension
      {
        static void apply ( void *&v )
        {
          v = &GenericReferenceElements< ctype, dim >::general;
        }
      };

      template< int dim >
      static void *apply ()
      {
        void *v;
        GenericGeometry::ForLoop< Dimension, 0, dim >::apply( v );
        return v;
      }
    };

    void initReferenceElements ()
    {
      InitReferenceElements< double >::apply< 3 >();
      InitReferenceElements< float >::apply< 3 >();
    }

  }

}
