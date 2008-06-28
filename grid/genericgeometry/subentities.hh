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

    template< class Geometry, int codim, unsigned int i, int subcodim >
    class NumSubSubEntities;

    template< class Geometry, int codim, unsigned int i >
    class SubGeometry;



    // NumSubEntities
    // --------------

    template< int codim >
    class NumSubEntities< Point, codim >
    {
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

    public:
      enum { value = m + n };
    };



    // NumSubSubEntities
    // -----------------

    template< int codim, unsigned int i, int subcodim >
    class NumSubSubEntities< Point, codim, i, subcodim >
    {
      typedef Point Geometry;
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index." );

    public:
      enum { value = (codim+subcodim == 0) ? 1 : 0 };
    };

    template< class BaseGeometry, int codim, unsigned int i, int subcodim >
    class NumSubSubEntities< Prism< BaseGeometry >, codim, i, subcodim >
    {
      typedef Prism< BaseGeometry > Geometry;
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index." );

      enum { n = NumSubEntities< BaseGeometry, codim > :: value };
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

    public:
      enum
      {
        value = (i >= n)
                ? NumSubSubEntities< BaseGeometry, codim-1, (i-n), subcodim >
                :: value
                : NumSubSubEntities< BaseGeometry, codim, i, subcodim >
                :: value
                + 2 * NumSubSubEntities< BaseGeometry, codim, i, subcodim-1 >
                :: value
      };
    };

    template< class BaseGeometry, int codim, unsigned int i, int subcodim >
    class NumSubSubEntities< Pyramid< BaseGeometry >, codim, i, subcodim >
    {
      typedef Pyramid< BaseGeometry > Geometry;
      dune_static_assert( (i < NumSubEntities< Geometry, codim > :: value),
                          "Invalid subentity index." );

      enum { n = NumSubEntities< BaseGeometry, codim > :: value };
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

    public:
      enum
      {
        value = (i < m)
                ? NumSubSubEntities< BaseGeometry, codim-1, i, subcodim >
                :: value
                : NumSubSubEntities< BaseGeometry, codim, i-m, subcodim >
                :: value
                + NumSubSubEntities< BaseGeometry, codim, i-m, subcodim-1 >
                :: value
      };
    };



    // PrismSubGeometry
    // ----------------

    template< class BaseGeometry, int dim, int codim, int i >
    struct PrismSubGeometry
    {
      enum { n = NumSubEntities< BaseGeometry, codim > :: value };
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

      typedef typename TypeIf
      < (i < n),
      Prism< typename SubGeometry< BaseGeometry, codim, i > :: type >,
      typename SubGeometry< BaseGeometry, codim-1, (i-n)%m > :: type
      > :: type type;
    };

    template< class BaseGeometry, int dim, int i >
    struct PrismSubGeometry< BaseGeometry, dim, 0, i >
    {
      typedef Prism< BaseGeometry > type;
    };

    template< class BaseGeometry, int dim, int i >
    struct PrismSubGeometry< BaseGeometry, dim, dim, i >
    {
      typedef Point type;
    };



    // PyramidSubGeometry
    // ------------------

    template< class BaseGeometry, int dim, int codim, int i >
    struct PyramidSubGeometry
    {
      enum { m = NumSubEntities< BaseGeometry, codim-1 > :: value };

      typedef typename TypeIf
      < (i < m),
      typename SubGeometry< BaseGeometry, codim-1, i > :: type,
      Pyramid< typename SubGeometry< BaseGeometry, codim, i-m > :: type >
      > :: type type;
    };

    template< class BaseGeometry, int dim, int i >
    struct PyramidSubGeometry< BaseGeometry, dim, 0, i >
    {
      typedef Pyramid< BaseGeometry > type;
    };

    template< class BaseGeometry, int dim, int i >
    struct PyramidSubGeometry< BaseGeometry, dim, dim, i >
    {
      typedef Point type;
    };



    // SubGeometry
    // -----------

    template< int codim, unsigned int i >
    class SubGeometry< Point, codim, i >
    {
      dune_static_assert( (i < NumSubEntities< Point, codim > :: value),
                          "Invalid subentity index." );

    public:
      typedef Point type;
    };

    template< class BaseGeometry, int codim, unsigned int i >
    class SubGeometry< Prism< BaseGeometry >, codim, i >
    {
      dune_static_assert( (i < NumSubEntities< Prism< BaseGeometry >, codim > :: value),
                          "Invalid subentity index." );

    public:
      typedef typename PrismSubGeometry
      < BaseGeometry, codim, Prism< BaseGeometry > :: dimension, i > :: type type;
    };

    template< class BaseGeometry, int codim, unsigned int i >
    class SubGeometry< Pyramid< BaseGeometry >, codim, i >
    {
      dune_static_assert( (i < NumSubEntities< Pyramid< BaseGeometry >, codim > :: value),
                          "Invalid subentity index." );

    public:
      typedef typename PyramidSubGeometry
      < BaseGeometry, codim, Pyramid< BaseGeometry > :: dimension, i > :: type type;
    };

  }

}

#endif
