// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_MISC_HH
#define DUNE_ALBERTA_MISC_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/genericgeometry/misc.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // Exceptions
  // ----------

  class AlbertaError
    : public Exception
  {};

  class AlbertaIOError
    : public IOError
  {};



  namespace Alberta
  {

    // Import Types
    // ------------

    using GenericGeometry::ForLoop;

    typedef ALBERTA REAL Real;
    typedef ALBERTA REAL_D GlobalVector;

    typedef ALBERTA MESH Mesh;
    typedef ALBERTA MACRO_EL MacroElement;
    typedef ALBERTA EL Element;
#if DUNE_ALBERTA_VERSION < 0x200
    typedef ALBERTA BOUNDARY Boundary;
#endif

    typedef ALBERTA FE_SPACE DofSpace;



    // NumSubEntities
    // --------------

    template< int dim, int codim >
    struct NumSubEntities;

    template< int dim >
    struct NumSubEntities< dim, 0 >
    {
      static const int value = 1;
    };

    template< int dim >
    struct NumSubEntities< dim, dim >
    {
      static const int value = dim+1;
    };

    template<>
    struct NumSubEntities< 2, 1 >
    {
      static const int value = 3;
    };

    template<>
    struct NumSubEntities< 3, 1 >
    {
      static const int value = 4;
    };

    template<>
    struct NumSubEntities< 3, 2 >
    {
      static const int value = 6;
    };



    // CodimType
    // ---------

    template< int dim, int codim >
    struct CodimType;

    template< int dim >
    struct CodimType< dim, 0 >
    {
      static const int value = CENTER;
    };

    template< int dim >
    struct CodimType< dim, dim >
    {
      static const int value = VERTEX;
    };

    template<>
    struct CodimType< 2, 1 >
    {
      static const int value = EDGE;
    };

    template<>
    struct CodimType< 3, 1 >
    {
      static const int value = FACE;
    };

    template<>
    struct CodimType< 3, 2 >
    {
      static const int value = EDGE;
    };

  }

}

#endif // #if HAVE_ALBERTA

#endif
