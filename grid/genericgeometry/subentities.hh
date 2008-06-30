// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBENTITIES_HH
#define DUNE_GENERICGEOMETRY_SUBENTITIES_HH

#include <dune/common/static_assert.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/geometrytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Geometry, int codim >
    class NumSubEntities;

    template< class Geometry, int codim, unsigned int i >
    class SubGeometry;



    // NumSubEntities
    // --------------

    template< int codim >
    class NumSubEntities< Point, codim >
    {
    public:
      enum { value = (codim == 0) ? 1 : 0 };
    };


    template< class BaseGeometry, int codim >
    class NumSubEntities< Prism< BaseGeometry >, codim >
    {
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

    public:
      enum { value = n + 2*m };
    };

    template< class BaseGeometry, int codim >
    struct NumSubEntities< Pyramid< BaseGeometry >, codim >
    {
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

      enum { dimension = Pyramid< BaseGeometry > :: dimension };

    public:
      enum { value = (codim == dimension ? m+1 : m+n) };
    };



    // SubGeometry
    // -----------

    template< int codim, unsigned int i >
    class SubGeometry< Point, codim, i >
    {
    public:
      typedef void type;
    };

    template<>
    class SubGeometry< Point, 0, 0 >
    {
    public:
      typedef Point type;
    };

    template< class BaseGeometry, int codim, unsigned int i >
    class SubGeometry< Prism< BaseGeometry >, codim, i >
    {
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };

      enum { j = (i < n+m ? i-n : i-(n+m)) };

    public:
      typedef typename TypeIf
      < (i < n),
      Prism< typename SubGeometry< BaseGeometry, codim, i > :: type >,
      typename SubGeometry< BaseGeometry, codim-1, j > :: type
      > :: type type;
    };

    template< class BaseGeometry, int codim, unsigned int i >
    class SubGeometry< Pyramid< BaseGeometry >, codim, i >
    {
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

      enum { dimension = Pyramid< BaseGeometry > :: dimension };

      typedef typename TypeIf
      < (codim == dimension),
          Point,
          Pyramid< typename SubGeometry< BaseGeometry, codim, i-m > :: type >
      > :: type pyramid_type;

    public:
      typedef typename TypeIf
      < (i < m),
      typename SubGeometry< BaseGeometry, codim-1, i > :: type,
      pyramid_type
      > :: type type;
    };

  }

}

#endif
