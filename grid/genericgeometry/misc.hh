// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MISC_HH
#define DUNE_GENERICGEOMETRY_MISC_HH

namespace Dune
{

  namespace GenericGeometry
  {

    template< unsigned int n >
    struct Faculty
    {
      enum { value = n * Faculty< n-1 > :: value };
    };

    template<>
    struct Faculty< 0 >
    {
      enum { value = 1 };
    };



    template< bool condition, class TypeTrue, class TypeFalse >
    struct TypeIf;

    template< class TypeTrue, class TypeFalse >
    struct TypeIf< false, TypeTrue, TypeFalse >
    {
      typedef TypeFalse type;
    };

    template< class TypeTrue, class TypeFalse >
    struct TypeIf< true, TypeTrue, TypeFalse >
    {
      typedef TypeTrue type;
    };

  }

}

#endif
