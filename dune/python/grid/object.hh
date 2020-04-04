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
      inline static typename GridObject::GridPartType::GridViewType
      gridView ( const GridObject &gridObject, PriorityTag< 2 > )
      {
        typedef typename GridObject::GridPartType::GridViewType GridView;
        return static_cast< GridView >( gridObject.gridPart() );
      }

      template< class GridObject >
      inline static decltype( std::declval< const GridObject & >().entitySet().gridView() )
      gridView ( const GridObject &gridObject, PriorityTag< 1 > )
      {
        return gridObject.entitySet().gridView();
      }

      template< class GridObject >
      inline static decltype( std::declval< const GridObject & >().gridView() )
      gridView ( const GridObject &gridObject, PriorityTag< 0 > )
      {
        return gridObject.gridView();
      }

    } // namespace detail



    // gridView
    // --------

    template< class GridObject >
    inline static decltype( auto ) gridView ( const GridObject &gridObject )
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
