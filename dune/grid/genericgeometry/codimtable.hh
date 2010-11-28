// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GENERICGEOMETRY_CODIMTABLE_HH
#define DUNE_GENERICGEOMETRY_CODIMTABLE_HH

#include <dune/common/typetraits.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< template< int > class Element, int codim >
    class CodimTableStorage
    {
      CodimTableStorage< Element, codim - 1 > map_;
      Element< codim > element_;

    public:
      CodimTableStorage ()
        : map_(),
          element_()
      {}

      CodimTableStorage ( const CodimTableStorage &other )
        : map_( other.map_ ),
          element_( other.element_ )
      {}

      const CodimTableStorage &operator= ( const CodimTableStorage &other )
      {
        map_ = other.map_;
        element_ = other.element_;
        return *this;
      }

      template< int cd >
      const Element< cd > &
      operator[] ( const integral_constant< int, cd > codimVariable ) const
      {
        return map_[ codimVariable ];
      }

      template< int cd >
      Element< cd > &
      operator[] ( const integral_constant< int, cd > &codimVariable )
      {
        return map_[ codimVariable ];
      }

      const Element< codim > &
      operator[] ( const integral_constant< int, codim > &codimVariable ) const
      {
        return element_;
      }

      Element< codim > &
      operator[] ( const integral_constant< int, codim > &codimVariable )
      {
        return element_;
      }
    };



    template< template< int > class Element >
    class CodimTableStorage< Element, 0 >
    {
      Element< 0 > element_;

    public:
      CodimTableStorage ()
        : element_()
      {}

      CodimTableStorage ( const CodimTableStorage &other )
        : element_( other.element_ )
      {}

      CodimTableStorage &operator= ( const CodimTableStorage &other )
      {
        element_ = other.element_;
        return *this;
      }

      const Element< 0 > &
      operator[] ( const integral_constant< int, 0 > codimVaraible ) const
      {
        return element_;
      }

      Element< 0 > &
      operator[] ( const integral_constant< int, 0 > codimVaraible )
      {
        return element_;
      }
    };



    template< template< int > class Element, int dim >
    class CodimTable
    {
      CodimTableStorage< Element, dim > map_;

    public:
      CodimTable ()
      {}

      CodimTable ( const CodimTable &other )
        : map_( other.map_ )
      {}

      const CodimTable &operator= ( const CodimTable &other )
      {
        map_ = other.map_;
        return *this;
      }

      template< int codim >
      const Element< codim > &
      operator[] ( const integral_constant< int, codim > codimVariable ) const
      {
        return map_[ codimVariable ];
      }

      template< int codim >
      Element< codim > &
      operator[] ( const integral_constant< int, codim > codimVariable )
      {
        return map_[ codimVariable ];
      }
    };

  }

}

#endif
