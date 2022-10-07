// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_FACTORY_HH
#define DUNE_PYTHON_GRID_FACTORY_HH

#include <cstddef>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dune/common/visibility.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/gridview.hh>

#include <dune/python/pybind11/extensions.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dune
{

  namespace Python
  {

    // BoundarySegment
    // ---------------

    template< int dim, int dimworld >
    struct DUNE_PRIVATE BoundarySegment final
      : public Dune::BoundarySegment< dim, dimworld >
    {
      explicit BoundarySegment ( pybind11::function parametrization )
        : parametrization_( std::move( parametrization ) )
      {}

      FieldVector< double, dimworld > operator() ( const FieldVector< double, dim-1 > &local ) const
      {
        return parametrization_( local ).template cast< FieldVector< double, dimworld > >();
      }

    private:
      pybind11::function parametrization_;
    };



    namespace detail
    {

      namespace GridFactory
      {

        // insertVertices
        // --------------

        template< class T, class Grid >
        inline static void insertVertices ( const pybind11::buffer_info &info, pybind11::format_descriptor< T >, Dune::GridFactory< Grid > &factory )
        {
          const int dimWorld = Grid::dimensionworld;

          for( int i = 0; i < info.shape[ 0 ]; ++i )
          {
            const std::size_t offset = i * (info.strides[ 0 ] / sizeof( T ));
            FieldVector< typename Grid::ctype, dimWorld > x;
            for( int j = 0; j < dimWorld; ++j )
              x[ j ] = static_cast< T * >( info.ptr )[ offset + j * (info.strides[ 1 ] / sizeof( T )) ];
            factory.insertVertex( x );
          }
        }


        template< class Grid >
        inline static void insertVertices ( const pybind11::buffer_info &info, Dune::GridFactory< Grid > &factory )
        {
          const int dimWorld = Grid::dimensionworld;
          if( (info.ndim != 2) || (info.shape[ 1 ] != dimWorld) )
            throw std::invalid_argument( "vertex array must be of shape (*, " + std::to_string( dimWorld ) + ")" );
          pybind11::handle_buffer_format(info, [ &info, &factory ] ( auto format ) { insertVertices( info, format, factory ); } );
        }


        template< class Grid >
        inline static void insertVertices ( pybind11::list data, Dune::GridFactory< Grid > &factory )
        {
          typedef FieldVector< typename Grid::ctype, Grid::dimensionworld > GlobalCoordinate;
          for( pybind11::handle item : data )
            factory.insertVertex( item.template cast< GlobalCoordinate >() );
        }


        template< class Grid >
        inline static void insertVertices ( pybind11::detail::item_accessor data, Dune::GridFactory< Grid > &factory )
        {
          try
          {
            insertVertices( data.template cast< pybind11::buffer >().request(), factory );
            return;
          }
          catch( const pybind11::error_already_set & )
          {}

          try
          {
            insertVertices( data.template cast< pybind11::list >(), factory );
            return;
          }
          catch( const pybind11::cast_error & )
          {}

          throw std::invalid_argument( "Incompatible array type for vertices." );
        }



        // insertElements
        // --------------

        template< class T, class Grid >
        inline static std::enable_if_t< std::is_integral< T >::value >
        insertElements ( GeometryType type, const pybind11::buffer_info &info, pybind11::format_descriptor< T >, Dune::GridFactory< Grid > &factory )
        {
          const int dimGrid = Grid::dimension;

          const int numVertices = ReferenceElements< typename Grid::ctype, dimGrid >::general( type ).size( dimGrid );
          if( (info.ndim != 2) || (info.shape[ 1 ] != numVertices) )
          {
            std::ostringstream msg;
            msg << "buffer for geometry type " << type << " must be of shape (*, " << numVertices << ")\n";
            if (info.ndim != 2)
              msg << " shape dimension is not 2 but " << info.ndim;
            else
              msg << " shape dimension is 2 but shape[1] is " << info.shape[ 1 ];
            throw std::invalid_argument( msg.str() );
          }

          std::vector< unsigned int > vertices( numVertices );
          for( int i = 0; i < info.shape[ 0 ]; ++i )
          {
            const std::size_t offset = i * (info.strides[ 0 ] / sizeof( T ));
            for( int j = 0; j < numVertices; ++j )
              vertices[ j ] = static_cast< T * >( info.ptr )[ offset + j * (info.strides[ 1 ] / sizeof( T )) ];
            factory.insertElement( type, vertices );
          }
        }

        template< class T, class Grid >
        inline static std::enable_if_t< !std::is_integral< T >::value >
        insertElements ( GeometryType type, const pybind11::buffer_info &info, pybind11::format_descriptor< T >, Dune::GridFactory< Grid > &factory )
        {
          std::ostringstream msg;
          msg << "Incompatible buffer format in array for geometry type " << type << ": '" << info.format << "'.";
          throw std::invalid_argument( msg.str() );
        }


        template< class Grid >
        inline static void insertElements ( GeometryType type, const pybind11::buffer_info &info, Dune::GridFactory< Grid > &factory )
        {
          pybind11::handle_buffer_format( info, [ type, &info, &factory ] ( auto format ) { insertElements( type, info, format, factory ); } );
        }


        template< class Grid >
        inline static void insertElements ( GeometryType type, const pybind11::list &data, Dune::GridFactory< Grid > &factory )
        {
          for( pybind11::handle item : data )
            factory.insertElement( type, item.template cast< std::vector< unsigned int > >() );
        }


        template< class Grid >
        inline static void insertElements ( GeometryType type, pybind11::detail::item_accessor data, Dune::GridFactory< Grid > &factory )
        {
          try
          {
            insertElements( type, data.template cast< pybind11::buffer >().request(), factory );
            return;
          }
          catch( const pybind11::error_already_set & )
          {}

          try
          {
            insertElements( type, data.template cast< pybind11::list >(), factory );
            return;
          }
          catch( const pybind11::cast_error & )
          {}

          std::ostringstream msg;
          msg << "Invalid array type for geometry type " << type << ".";
          throw std::invalid_argument( msg.str() );
        }


        template< class Grid >
        inline static void insertElements ( pybind11::list data, Dune::GridFactory< Grid > &factory )
        {
          for( pybind11::handle hItem : data )
          {
            auto item = pybind11::reinterpret_borrow< pybind11::tuple >( hItem );
            if( item.size() == 2 )
              factory.insertElement( pybind11::cast< GeometryType >( item[ 0 ] ), pybind11::cast< std::vector< unsigned int > >( item[ 1 ] ) );
            else
              throw std::invalid_argument( "Element tuple must be of length 2." );
          }
        }

        template< class Grid >
        inline static void insertElements ( pybind11::detail::item_accessor data, Dune::GridFactory< Grid > &factory )
        {
          insertElements( pybind11::cast< pybind11::list >( data ), factory );
        }



        // insertBoundaries
        // ----------------

        template< class Grid >
        inline static void insertBoundaries ( pybind11::list data, Dune::GridFactory< Grid > &factory )
        {
          const int dimGrid = Grid::dimension;
          const int dimWorld = Grid::dimensionworld;

          for( pybind11::handle item : data )
          {
            auto vertices = item.template cast< std::vector< unsigned int > >();

            pybind11::function f;
            try
            {
              f = item.template cast< pybind11::function >();
            }
            catch( const pybind11::cast_error & )
            {}

            if( f )
              factory.insertBoundarySegment( vertices, std::make_shared< BoundarySegment< dimGrid, dimWorld > >( f ) );
            else
              factory.insertBoundarySegment( vertices );
          }
        }

        template< class Grid >
        inline static void insertBoundaries ( pybind11::detail::item_accessor data, Dune::GridFactory< Grid > &factory )
        {
          insertBoundaries( pybind11::cast< pybind11::list >( data ), factory );
        }

      } // namespace GridFactory

    } // namespace detail



    // fillGridFactory
    // ---------------

    template< class Grid >
    inline void fillGridFactory ( const pybind11::dict &dict, Dune::GridFactory< Grid > &factory )
    {
      const int dimGrid = Grid::dimension;

      if( dict.contains( "vertices" ) )
        detail::GridFactory::insertVertices( dict[ "vertices" ], factory );
      else
        throw std::invalid_argument( "Missing Key: 'vertices'" );

      if( dict.contains( "elements" ) )
        detail::GridFactory::insertElements( dict[ "elements" ], factory );

      if( dict.contains( "lines" ) )
        detail::GridFactory::insertElements( GeometryTypes::line, dict[ "lines" ], factory );
      if( dict.contains( "line" ) )
        detail::GridFactory::insertElements( GeometryTypes::line, dict[ "line" ], factory );

      if( dict.contains( "triangles" ) )
        detail::GridFactory::insertElements( GeometryTypes::triangle, dict[ "triangles" ], factory );
      if( dict.contains( "triangle" ) )
        detail::GridFactory::insertElements( GeometryTypes::triangle, dict[ "triangles" ], factory );

      if( dict.contains( "tetrahedra" ) )
        detail::GridFactory::insertElements( GeometryTypes::tetrahedron, dict[ "tetrahedra" ], factory );
      if( dict.contains( "simplices" ) )
        detail::GridFactory::insertElements( GeometryTypes::simplex( dimGrid ), dict[ "simplices" ], factory );
      if( dict.contains( "tetra" ) )
      {
        std::cout << "reading tetras\n";
        detail::GridFactory::insertElements( GeometryTypes::simplex( dimGrid ), dict[ "tetra" ], factory );
      }

      if( dict.contains( "quadrilaterals" ) )
        detail::GridFactory::insertElements( GeometryTypes::quadrilateral, dict[ "quadrilaterals" ], factory );
      if( dict.contains( "hexahedra" ) )
        detail::GridFactory::insertElements( GeometryTypes::hexahedron, dict[ "hexahedra" ], factory );
      if( dict.contains( "cubes" ) )
        detail::GridFactory::insertElements( GeometryTypes::cube( dimGrid ), dict[ "cubes" ], factory );

      if( dict.contains( "prisms" ) )
        detail::GridFactory::insertElements( GeometryTypes::prism, dict[ "prisms" ], factory );
      if( dict.contains( "pyramid" ) )
        detail::GridFactory::insertElements( GeometryTypes::prism, dict[ "pyramids" ], factory );

      if( dict.contains( "boundaries" ) )
        detail::GridFactory::insertBoundaries( dict[ "boundaries" ], factory );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_FACTORY_HH
