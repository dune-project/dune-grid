// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GENERICGEOMETRY_CODIMTABLE_HH
#define DUNE_GENERICGEOMETRY_CODIMTABLE_HH

#include <dune/common/typetraits.hh>
#include <dune/common/tupleutility.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< template< int > class Element, int dim >
    class CodimTable
    {
      friend class CodimTable< Element, dim+1 >;

      typedef typename PushBackTuple<
          typename CodimTable< Element, dim-1 >::ElementTuple,
          Element< dim > >::type ElementTuple;

      ElementTuple map_;

    public:

      template< int codim >
      const Element< codim > &
      operator[] ( const integral_constant< int, codim > codimVariable ) const
      {
        return Dune::get<codim>(map_);
      }

      template< int codim >
      Element< codim > &
      operator[] ( const integral_constant< int, codim > codimVariable )
      {
        return Dune::get<codim>(map_);
      }
    };


    template< template< int > class Element>
    class CodimTable< Element, -1 >
    {
      friend class CodimTable< Element, 0 >;
      typedef typename Dune::tuple<> ElementTuple;
    };

  }

}

#endif
