// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_GRID_OBJECT_HH
#define DUNE_PYTHON_GRID_OBJECT_HH

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {
      template< class GridObject >
      inline static const decltype( std::declval< const GridObject & >().gridView() )&
      gridView ( const GridObject &gridObject, PriorityTag< 3 > )
      {
        return gridObject.gridView();
      }
      template< class GridObject >
      inline static const decltype( std::declval< const GridObject & >().entitySet().gridView() )&
      gridView ( const GridObject &gridObject, PriorityTag< 2 > )
      {
        return gridObject.entitySet().gridView();
      }
      template< class GridObject >
      inline static const decltype( std::declval< const typename GridObject::GridPartType::GridViewType & >() )&
      gridView ( const GridObject &gridObject, PriorityTag< 1 > )
      {
        return gridObject.gridPart();
      }
    } // namespace detail


    // gridView
    // --------

    template< class GridObject >
    inline static const auto& gridView ( const GridObject &gridObject )
    {
      return detail::gridView( gridObject, PriorityTag< 42 >() );
    }

    namespace detail
    {

      namespace GridObjectTraits
      {

        using Dune::Python::gridView;


        template< class GridObject >
        using GridView = std::decay_t< decltype( gridView( std::declval< const GridObject & >() ) ) >;


        template< class GridObject >
        typename GridObject::EntitySet::Element element ( const GridObject &, PriorityTag< 1 > );

        template< class GridObject >
        typename GridView< GridObject >::template Codim< 0 >::Entity element ( const GridObject &, PriorityTag< 0 > );


        template< class GridObject >
        typename GridObject::EntitySet::LocalCoordinate localCoordinate ( const GridObject &, PriorityTag< 1 > );

        template< class GridObject >
        typename GridView< GridObject >::template Codim< 0 >::Geometry::LocalCoordinate localCoordinate ( const GridObject &, PriorityTag< 0 > );

      } // namespace GridObjectTraits

    } // namespace detail



    // GridObjectTraits
    // ----------------

    template< class GridObject >
    struct GridObjectTraits
    {
      typedef decltype( detail::GridObjectTraits::element( std::declval< const GridObject & >(), PriorityTag< 42 >() ) ) Element;

      typedef decltype( detail::GridObjectTraits::localCoordinate( std::declval< const GridObject & >(), PriorityTag< 42 >() ) ) LocalCoordinate;
    };

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_OBJECT_HH
