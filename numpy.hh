#ifndef DUNE_PYTHON_GRID_NUMPY_HH
#define DUNE_PYTHON_GRID_NUMPY_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <array>
#include <vector>
#include <map>

#include <dune/common/ftraits.hh>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/virtualrefinement.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/python/common/getdimension.hh>
#include <dune/python/grid/object.hh>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace CorePy
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
    tesselate ( const GridView &gridView, int level, PartitionSet< partitions > ps )
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
        for( auto it = refinement.vBegin( level ), end = refinement.vEnd( level ); it != end; ++it )
          coords.push_back( geometry.global( it.coords() ) );

        // get simplices
        for( auto it = refinement.eBegin( level ), end = refinement.eEnd( level ); it != end; ++it )
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
    inline static std::vector< pybind11::array_t< typename GridView::ctype > >
    polygons ( const GridView &gridView, PartitionSet< partitions > ps )
    {
      typedef typename GridView::ctype ctype;

      const std::size_t dimGrid = GridView::dimension;
      const std::size_t dimWorld = GridView::dimensionworld;

      std::map< GeometryType, std::vector< std::vector< FieldVector< ctype, dimWorld > > > > coords;

      for( const auto &element : elements( gridView, ps ) )
      {
        const auto &refinement = buildRefinement< dimGrid, double >( element.type(), element.type() );
        std::vector<FieldVector<ctype,dimWorld>> poly;

        // get coordinates
        const auto geometry = element.geometry();
        for( auto it = refinement.vBegin( 0 ), end = refinement.vEnd( 0 ); it != end; ++it )
          poly.push_back( geometry.global( it.coords() ) );

        // Martin: This seems to be limited to two spatial dimensions.
        if( element.type().isCube() )
          std::swap( poly[ 0 ], poly[ 1 ] );

        coords[ element.type() ].push_back( poly );
      }

      std::vector< pybind11::array_t< typename GridView::ctype > > ret;
      // Martin: This will not work for true polygons, due to difference vertex sizes.
      //         A solution could be a NumPy array per vertex count instead of per geometry type.
      for( const auto &entry : coords )
        ret.push_back( makeNumPyArray< ctype >( entry.second, { entry.second.size(), entry.second[ 0 ].size(), dimWorld } ) );
      return ret;
    }

    template< class GridView, unsigned int partitions >
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tesselate ( const GridView &gridView, PartitionSet< partitions > ps )
    {
      return tesselate( gridView, 0, ps );
    }

    template< class GridView >
    inline static std::pair< pybind11::array_t< typename GridView::ctype >, pybind11::array_t< int > >
    tesselate ( const GridView &gridView, int level = 0 )
    {
      return tesselate( gridView, level, Partitions::all );
    }

    template< class GridView >
    inline static std::vector< pybind11::array_t< typename GridView::ctype > >
    polygons ( const GridView &gridView )
    {
      return polygons( gridView, Partitions::all );
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
      for( const auto &element : elements( gv, ps ) )
      {
        lf.bind( element );
        const auto &refinement = buildRefinement< Element::mydimension, ctype >( element.type(), GeometryTypes::simplex( Element::dimension ) );
        for( auto it = refinement.vBegin( level ), end = refinement.vEnd( level ); it != end; ++it )
          values.push_back( lf( it.coords() ) );
        lf.unbind();
      }

      return CorePy::makeNumPyArray< typename FieldTraits< Range >::field_type >( values, { values.size(), GetDimension<Range>::value } );
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
      for( const Element &element : entities( gv, Dune::Codim< Element::codimension >(), ps ) )
      {
        lf.bind( element );
        const auto &refinement = buildRefinement< Element::mydimension, ctype >( element.type(), GeometryTypes::simplex( Element::mydimension ) );
        for( auto it = refinement.eBegin( 0 ), end = refinement.eEnd( 0 ); it != end; ++it )
          values.push_back( lf( it.coords() ) );
        lf.unbind();
      }

      return CorePy::makeNumPyArray< typename FieldTraits< Range >::field_type >( values, { values.size(), GetDimension<Range>::value } );
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

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_NUMPY_HH
