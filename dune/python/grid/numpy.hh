// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_GRID_NUMPY_HH
#define DUNE_PYTHON_GRID_NUMPY_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <array>
#include <map>
#include <vector>

#include <dune/common/ftraits.hh>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/virtualrefinement.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/python/common/getdimension.hh>
#include <dune/python/grid/object.hh>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    // External Forward Declarations
    // -----------------------------

    template< class GridFunction >
    struct GridFunctionTraits;



    // coordinates
    // -----------

    template< class GridView, class Mapper >
    inline static pybind11::array_t< typename GridView::ctype > coordinates ( const GridView &gridView, const Mapper &mapper )
    {
      typedef typename GridView::ctype ctype;

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( mapper.size() ), static_cast< std::size_t >( GridView::dimensionworld ) };
      const std::vector< std::size_t > stride{ GridView::dimensionworld * sizeof( ctype ), sizeof( ctype ) };
      pybind11::array_t< ctype > coords( pybind11::buffer_info( nullptr, sizeof( ctype ), pybind11::format_descriptor< ctype >::value, 2, shape, stride ) );

      pybind11::buffer_info info = coords.request( true );
      for( const auto &vertex : vertices( gridView, Partitions::all ) )
      {
        typename Mapper::Index index;
        if( !mapper.contains( vertex, index ) )
          continue;

        const auto x = vertex.geometry().center();
        std::copy( x.begin(), x.end(), static_cast< ctype * >( info.ptr ) + GridView::dimensionworld * index );
      }

      return coords;
    }

    template< class GridView >
    inline static pybind11::array_t< typename GridView::ctype > coordinates ( const GridView &gridView )
    {
      MultipleCodimMultipleGeomTypeMapper< GridView > mapper( gridView, mcmgVertexLayout() );
      return coordinates( gridView, mapper );
    }



    // flatCopy
    // --------

    template< class In, class T, class = std::enable_if_t< std::is_convertible< In, T >::value > >
    T *flatCopy ( const In &in, T *out )
    {
      *out = in;
      return ++out;
    }

    template< class In, class T, class = decltype( std::declval< const In & >().begin() ), class = std::enable_if_t< !std::is_convertible< In, T >::value > >
    T *flatCopy ( const In &in, T *out )
    {
      for( auto it = in.begin(), end = in.end(); it != end; ++it )
        out = flatCopy( *it, out );
      return out;
    }



    // makeNumPyArray
    // --------------

    template< class T, class In >
    pybind11::array_t< T > makeNumPyArray ( const In &in, const std::vector< std::size_t > &shape )
    {
      const std::size_t dim = shape.size();

      std::size_t size = sizeof( T );
      std::vector< std::size_t > stride( dim );
      for( std::size_t i = dim; i-- > 0; )
      {
        stride[ i ] = size;
        size *= shape[ i ];
      }

      pybind11::array_t< T > result( pybind11::buffer_info( nullptr, sizeof( T ), pybind11::format_descriptor< T >::value, dim, shape, stride ) );

      pybind11::buffer_info info = result.request( true );
      flatCopy( in, static_cast< T * >( info.ptr ) );

      return result;
    }



    // tessellate
    // ----------
    template< class GridView, unsigned int partitions >
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tessellate ( const GridView &gridView, RefinementIntervals intervals, PartitionSet< partitions > ps )
    {
      typedef typename GridView::ctype ctype;

      const std::size_t dimGrid = GridView::dimension;
      const std::size_t dimWorld = GridView::dimensionworld;

      std::vector< FieldVector< ctype, dimWorld > > coords;
      std::vector< std::array< int, dimGrid+1 > > simplices;

      for( const auto &element : elements( gridView, ps ) )
      {
        const auto &refinement = buildRefinement< dimGrid, double >( element.type(), GeometryTypes::simplex( dimGrid ) );
        const std::size_t offset = coords.size();

        // get coordinates
        const auto geometry = element.geometry();
        for( auto it = refinement.vBegin( intervals ), end = refinement.vEnd( intervals ); it != end; ++it )
          coords.push_back( geometry.global( it.coords() ) );

        // get simplices
        for( auto it = refinement.eBegin( intervals ), end = refinement.eEnd( intervals ); it != end; ++it )
        {
          std::array< int, dimGrid+1 > simplex;
          auto indices = it.vertexIndices();
          assert( indices.size() == simplex.size() );
          std::transform( indices.begin(), indices.end(), simplex.begin(), [ offset ] ( std::size_t i ) { return (i + offset); } );
          simplices.push_back( simplex );
        }
      }

      return std::make_pair( makeNumPyArray< ctype >( coords, { coords.size(), dimWorld } ), makeNumPyArray< int >( simplices, { simplices.size(), dimGrid+1 } ) );
    }

    template< class GridView, unsigned int partitions >
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tessellate ( const GridView &gridView, int level, PartitionSet< partitions > ps )
    {
      return tessellate( gridView, refinementLevels( level ), ps );
    }

    template< class GridView, unsigned int partitions >
    inline static std::vector< pybind11::array_t< typename GridView::ctype > >
    polygons ( const GridView &gridView, PartitionSet< partitions > ps )
    {
      typedef typename GridView::ctype ctype;

      const std::size_t dimWorld = GridView::dimensionworld;

      std::map< std::size_t, std::vector< std::vector< FieldVector< ctype, dimWorld > > > > coords;

      for( const auto &element : elements( gridView, ps ) )
      {
        // get coordinates
        const auto geometry = element.geometry();
        const std::size_t corners = geometry.corners();
        std::vector< FieldVector< ctype, dimWorld > > poly( corners );
        for( std::size_t i = 0; i != corners; ++i )
          poly[ i ] = geometry.corner( i );

        // Martin: This seems to be limited to two spatial dimensions.
        if( element.type().isCube() )
          std::swap( poly[ 0 ], poly[ 1 ] );

        coords[ corners ].push_back( poly );
      }

      std::vector< pybind11::array_t< typename GridView::ctype > > ret;
      for( const auto &entry : coords )
        ret.push_back( makeNumPyArray< ctype >( entry.second, { entry.second.size(), entry.first, dimWorld } ) );
      return ret;
    }

    template< class GridFunction, unsigned int partitions >
    inline static
    // std::pair<
    //    std::vector< pybind11::array_t< typename GridFunction::GridView::ctype > >,
    //    pybind11::array_t< typename GridFunction::GridView::ctype >
    //    >
    auto polygonData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      typedef typename GridFunction::GridView GridView;
      typedef typename GridView::ctype ctype;

      const std::size_t dimWorld = GridView::dimensionworld;
      typedef typename GridFunctionTraits< GridFunction >::Range Range;

      std::map< std::size_t, std::pair<
           std::vector< std::vector< FieldVector< ctype, dimWorld > > >,
           std::vector< Range >
           > > coords;

      const auto &gv = gridView( gridFunction );
      auto lf = localFunction( gridFunction );

      for( const auto &element : elements( gv, ps ) )
      {
        // get coordinates
        const auto geometry = element.geometry();
        const std::size_t corners = geometry.corners();
        std::vector< FieldVector< ctype, dimWorld > > poly( corners );
        for( std::size_t i = 0; i != corners; ++i )
          poly[ i ] = geometry.corner( i );

        // Martin: This seems to be limited to two spatial dimensions.
        if( element.type().isCube() )
          std::swap( poly[ 0 ], poly[ 1 ] );

        lf.bind(element);
        coords[ corners ].first.push_back( poly );
        // for polygons we can't use a reference element but one could use
        // the following for elements with type not none:
        // auto ref = ReferenceElements< typename GridView::ctype, dimGrid >::general( element.type() );
        // coords[ corners ].second.push_back( lf( ref.position(0,0) ) );
        coords[ corners ].second.push_back( lf( geometry.local( geometry.center() )  ) );
        lf.unbind();
      }

      std::vector< pybind11::array_t< typename GridView::ctype > > ret;
      std::vector< pybind11::array_t< typename GridView::ctype > > values;
      for( const auto &entry : coords )
      {
        ret.push_back( makeNumPyArray< ctype >( entry.second.first, { entry.second.first.size(), entry.first, dimWorld } ) );
        values.push_back( Python::makeNumPyArray< typename FieldTraits< Range >::field_type >( entry.second.second, { entry.second.second.size(), GetDimension<Range>::value } ) );
      }
      return std::make_pair(ret, values);
    }

    template< class GridView, unsigned int partitions >
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tessellate ( const GridView &gridView, PartitionSet< partitions > ps )
    {
      return tessellate( gridView, 0, ps );
    }

    template< class GridView >
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tessellate ( const GridView &gridView, int level = 0 )
    {
      return tessellate( gridView, level, Partitions::all );
    }

    template< class GridView >
    inline static std::vector< pybind11::array_t< typename GridView::ctype > >
    polygons ( const GridView &gridView )
    {
      return polygons( gridView, Partitions::all );
    }

    template< class GridFunction >
    inline static
    // std::pair<
    //    std::vector< pybind11::array_t< typename GridFunction::GridView::ctype > >,
    //    pybind11::array_t< typename GridFunction::GridView::ctype >
    //    >
    auto polygonData ( const GridFunction &gridFunction )
    {
      return polygonData( gridFunction, Partitions::all );
    }

    // pointData
    // ---------

    template< class GridFunction, unsigned int partitions >
    inline static auto pointData ( const GridFunction &gridFunction, int level, PartitionSet< partitions > ps )
    {
      typedef typename GridFunctionTraits< GridFunction >::Element Element;
      typedef typename GridFunctionTraits< GridFunction >::LocalCoordinate LocalCoordinate;
      typedef typename GridFunctionTraits< GridFunction >::Range Range;

      typedef typename FieldTraits< LocalCoordinate >::field_type ctype;

      const auto &gv = gridView( gridFunction );

      auto lf = localFunction( gridFunction );
      std::vector< Range > values;
      auto refLevel = refinementLevels(level);
      for( const auto &element : elements( gv, ps ) )
      {
        lf.bind( element );
        const auto &refinement = buildRefinement< Element::mydimension, ctype >( element.type(), GeometryTypes::simplex( Element::dimension ) );
        for( auto it = refinement.vBegin( refLevel ), end = refinement.vEnd( refLevel ); it != end; ++it )
          values.push_back( lf( it.coords() ) );
        lf.unbind();
      }

      return Python::makeNumPyArray< typename FieldTraits< Range >::field_type >( values, { values.size(), GetDimension<Range>::value } );
    }

    template< class GridFunction, unsigned int partitions >
    inline static auto pointData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      return pointData( gridFunction, 0, ps );
    }

    template< class GridFunction >
    inline static auto pointData ( const GridFunction &gridFunction, int level = 0 )
    {
      return pointData( gridFunction, level, Partitions::all );
    }



    // cellData
    // --------

    template< class GridFunction, unsigned int partitions >
    inline static auto cellData ( const GridFunction &gridFunction, int level, PartitionSet< partitions > ps )
    {
      typedef typename GridFunctionTraits< GridFunction >::Element Element;
      typedef typename GridFunctionTraits< GridFunction >::LocalCoordinate LocalCoordinate;
      typedef typename GridFunctionTraits< GridFunction >::Range Range;

      typedef typename FieldTraits< LocalCoordinate >::field_type ctype;

      const auto &gv = gridView( gridFunction );

      auto lf = localFunction( gridFunction );
      std::vector< Range > values;
      auto refLevel = refinementLevels(0);
      for( const Element &element : entities( gv, Dune::Codim< Element::codimension >(), ps ) )
      {
        lf.bind( element );
        const auto &refinement = buildRefinement< Element::mydimension, ctype >( element.type(), GeometryTypes::simplex( Element::mydimension ) );
        for( auto it = refinement.eBegin( refLevel ), end = refinement.eEnd( refLevel ); it != end; ++it )
          values.push_back( lf( it.coords() ) );
        lf.unbind();
      }

      return Python::makeNumPyArray< typename FieldTraits< Range >::field_type >( values, { values.size(), GetDimension<Range>::value } );
    }

    template< class GridFunction, unsigned int partitions >
    inline static auto cellData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      return cellData( gridFunction, 0, ps );
    }

    template< class GridFunction >
    inline static auto cellData ( const GridFunction &gridFunction, int level = 0 )
    {
      return cellData( gridFunction, level, Partitions::all );
    }


    // deprecated tesselate (note spelling) functions
    template< class GridView, unsigned int partitions >
    [[deprecated("use 'tessellate' (note spelling)")]]
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tesselate ( const GridView &gridView, RefinementIntervals intervals, PartitionSet< partitions > ps )
    { return tessellate( gridView, intervals, ps); }
    template< class GridView, unsigned int partitions >
    [[deprecated("use 'tessellate' (note spelling)")]]
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tesselate ( const GridView &gridView, int level, PartitionSet< partitions > ps )
    { return tessellate( gridView, level, ps); }
    template< class GridView, unsigned int partitions >
    [[deprecated("use 'tessellate' (note spelling)")]]
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tesselate ( const GridView &gridView, PartitionSet< partitions > ps )
    { return tessellate( gridView, ps); }
    template< class GridView >
    [[deprecated("use 'tessellate' (note spelling)")]]
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tesselate ( const GridView &gridView, int level = 0 )
    { tessellate(gridView, level); }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_NUMPY_HH
