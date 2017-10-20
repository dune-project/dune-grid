// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_GRIDPART_HH
#define DUNE_PYTHON_GRID_GRIDPART_HH

#include <array>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/python/grid/entity.hh>
#include <dune/python/grid/function.hh>
#include <dune/python/grid/indexset.hh>
#include <dune/python/grid/intersection.hh>
#include <dune/python/grid/mapper.hh>
#include <dune/python/grid/numpy.hh>
#include <dune/python/grid/range.hh>
#include <dune/python/grid/vtk.hh>
#include <dune/python/grid/capabilities.hh>

#include <dune/python/pybind11/pybind11.h>

#include <iostream>

namespace Dune
{

  namespace Python
  {

    // registerPyGridViewIterators
    // ---------------------------

    template< class GridView, int codim >
    inline static auto registerPyGridViewIterator ( pybind11::handle scope, PriorityTag< 1 > )
      -> std::enable_if_t< Capabilities::canIterate< typename GridView::Grid, codim >::value, std::function< pybind11::object ( const GridView & ) > >
    {
      typedef PyGridViewIterator<GridView,codim> Iterator;
      auto entry = insertClass<Iterator>(scope, "EntityIterator_"+std::to_string(codim),
           GenerateTypeName("PyGridViewIterator",MetaType<GridView>(),codim),
           IncludeFiles{"dune/python/grid/range.hh"} );
      if (entry.second)
        registerPyIterator( scope, entry.first );

      return [] ( const GridView &gridView ) -> pybind11::object {
          return pybind11::cast( Iterator( gridView.template begin< codim >(), gridView.template end< codim >() ) );
        };
    }

    template< class GridView, int codim >
    inline static std::function< pybind11::object ( const GridView & ) > registerPyGridViewIterator ( pybind11::handle scope, PriorityTag< 0 > )
    {
      return [] ( const GridView &gridView ) -> pybind11::object {
          throw std::invalid_argument( "Iterators for codimension " + std::to_string( codim ) + " are not implemented." );
        };
    }

    template< class GridView, int codim >
    inline static std::function< pybind11::object ( const GridView & ) > registerPyGridViewIterator ( pybind11::handle scope )
    {
      return registerPyGridViewIterator< GridView, codim >( scope, PriorityTag< 42 >() );
    }

    template< class GridView, int... codim >
    inline static std::array< std::function< pybind11::object( const GridView & ) >, sizeof...( codim ) >
    registerPyGridViewIterators( pybind11::handle scope, std::integer_sequence< int, codim... > )
    {
      return { registerPyGridViewIterator< GridView, codim >( scope )... };
    }



    // registerPyGridViewParIterator
    // -----------------------------

    template< class GridView, int codim, int partition >
    inline static auto registerPyGridViewParIterator( pybind11::handle scope, PriorityTag< 1 > )
      -> std::enable_if_t< Capabilities::canIterate< typename GridView::Grid, codim >::value, std::function< pybind11::object ( const GridView & ) > >
    {
      typedef PyGridViewParIterator< GridView, codim, partition > Iterator;
      auto entry = insertClass<Iterator>(scope, "EntityIterator_"+std::to_string(codim)+std::to_string((int)partition),
           GenerateTypeName("PyGridViewParIterator",MetaType<GridView>(),codim,partition),
           IncludeFiles{"dune/python/grid/range.hh"} );
      if (entry.second)
        registerPyIterator< Iterator >( scope, entry.first );
      return [] ( const GridView &gridView ) -> pybind11::object {
          auto begin = gridView.template begin<codim, static_cast<Dune::PartitionIteratorType>(partition)>();
          auto end   = gridView.template end  <codim, static_cast<Dune::PartitionIteratorType>(partition)>();
          return pybind11::cast(PyGridViewParIterator<GridView, codim, partition>(begin, end));
      };
    }

    template< class GridView, int codim, int partition >
    inline static std::function< pybind11::object ( const GridView & ) > registerPyGridViewParIterator( pybind11::handle scope, PriorityTag< 0 > )
    {
      return [] ( const GridView &gridView ) -> pybind11::object {
          throw std::invalid_argument( "Iterators for codimension " + std::to_string( codim ) + " are not implemented." );
        };
    }

    template< class GridView, int codim, int partition >
    inline static std::function< pybind11::object ( const GridView & ) > registerPyGridViewParIterator( pybind11::handle scope )
    {
      return registerPyGridViewParIterator< GridView, codim, partition >( scope, PriorityTag< 42 >() );
    }

    template<class GridView, int codim, int... partitions>
    auto registerPyGridViewParIterators_(pybind11::handle scope, std::integer_sequence<int, partitions...>)
    {
      std::array<std::function<pybind11::object(const GridView&)>, sizeof...(partitions)>
        pyGridViewParIterators =
          {registerPyGridViewParIterator<GridView, codim, partitions>(scope)...};

      return pyGridViewParIterators;
    }

    template<class GridView, int... codim>
    auto registerPyGridViewParIterators(pybind11::handle scope, std::integer_sequence<int, codim...>)
    {
      constexpr unsigned int nPartitionType = 6; // see Dune::PartitionIteratorType enum

      std::array<
        std::array< std::function<pybind11::object(const GridView&)>, nPartitionType >,
        sizeof...(codim)
      >
        pyGridViewParIterators =
          {registerPyGridViewParIterators_<GridView, codim>(scope, std::make_integer_sequence<int, nPartitionType>())...};

      return pyGridViewParIterators;
    }



    // ProxyDataHandle for parallel communication
    // ------------------------------------------

    struct DUNE_PRIVATE ProxyDataHandle
      : public Dune::CommDataHandleIF< ProxyDataHandle, double >
    {
      ProxyDataHandle ( pybind11::object dataHandle )
        : contains_( method( dataHandle, "contains" ) ), fixedsize_( method( dataHandle, "fixedsize" ) ),
          size_( method( dataHandle, "size" ) ), gather_( method( dataHandle, "gather" ) ), scatter_( method( dataHandle, "scatter" ) )
      {}

      bool contains ( int dim, int codim ) const { return contains_( dim, codim ).cast< bool >(); }
      bool fixedsize ( int dim, int codim ) const { return fixedsize_( dim, codim ).cast< bool >(); }

      template< class Entity >
      std::size_t size ( const Entity &entity ) const
      {
        return size_( entity ).template cast< std::size_t >();
      }

      template< class Buffer, class Entity >
      void gather ( Buffer &buffer, const Entity &entity ) const
      {
        pybind11::list data = gather_( entity );
        for( const pybind11::handle &x : data )
          buffer.write( x.template cast< double >() );
      }

      template< class Buffer, class Entity >
      void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
      {
        pybind11::list data;
        for( std::size_t i = 0; i < n; ++i )
        {
          double x;
          buffer.read( x );
          data.append( pybind11::cast( x ) );
        }
        scatter_( entity, data );
      }

    private:
      pybind11::object method ( pybind11::handle dataHandle, const char *name )
      {
        pybind11::object method = dataHandle.attr( name );
        if( !method )
          throw std::invalid_argument( std::string( "Method \"" ) + name + std::string( "\" missing from data handle." ) );
        return method;
      }

      pybind11::object contains_, fixedsize_;
      pybind11::object size_, gather_, scatter_;
    };


    namespace detail
    {

      // registerGridViewConstructorFromGrid
      // -----------------------------------

      template< class GridView, class... options, std::enable_if_t< std::is_same< GridView, typename GridView::Grid::LeafGridView >::value, int > = 0 >
      void registerGridViewConstructorFromGrid ( pybind11::class_< GridView, options... > &cls, PriorityTag< 2 > )
      {
        cls.def( pybind11::init( [] ( typename GridView::Grid &grid ) { return new GridView( grid.leafGridView() ); } ), pybind11::keep_alive< 1, 2 >() );
      }

      template< class GridView, class... options, std::enable_if_t< std::is_constructible< GridView, typename GridView::Grid & >::value, int > = 0 >
      void registerGridViewConstructorFromGrid ( pybind11::class_< GridView, options... > &cls, PriorityTag< 1 > )
      {
        cls.def( pybind11::init( [] ( typename GridView::Grid &grid ) { return new GridView( grid ); } ), pybind11::keep_alive< 1, 2 >() );
      }

      template< class GridView, class... options >
      void registerGridViewConstructorFromGrid ( pybind11::class_< GridView, options... > &, PriorityTag< 0 > )
      {}

    } // namespace detail



    // registerGridView
    // ----------------

    template< class GridView, class... options >
    void registerGridView ( pybind11::handle scope, pybind11::class_< GridView, options... > cls )
    {
      using pybind11::operator""_a;

      typedef typename GridView::Grid Grid;
      const int dim = GridView::dimension;

      detail::registerGridViewConstructorFromGrid( cls, PriorityTag< 42 >() );

      cls.attr( "dimGrid" ) = pybind11::int_( static_cast< int >( GridView::dimension ) );
      cls.attr( "dimWorld" ) = pybind11::int_( static_cast< int >( GridView::dimensionworld ) );

      registerGridEntities< GridView >( cls );
      registerGridIntersection< GridView >( cls );

      registerGridViewIndexSet< GridView >( cls );
      registerMultipleCodimMultipleGeomTypeMapper< GridView >( cls );
      typedef MultipleCodimMultipleGeomTypeMapper<GridView> MCMGMapper;
      cls.def("mapper", [](GridView &self, pybind11::function contains)
          { return MCMGMapper(self,
                  [contains](Dune::GeometryType gt, int griddim)
                  { return (unsigned)(contains(gt).cast<int>()); } ); },
          pybind11::keep_alive<0,1>() );

      // STATIC, needed? - visibility problem?
      static const auto pyGridViewIterators = registerPyGridViewIterators< GridView >(cls, std::make_integer_sequence<int, dim+1>() );

      cls.def( "elements", [] ( const GridView &self ) {
            return pyGridViewIterators[ 0 ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "facets", [] ( const GridView &self ) {
            return pyGridViewIterators[ 1 ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "edges", [] ( const GridView &self ) {
            return pyGridViewIterators[ dim-1 ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "vertices", [] ( const GridView &self ) {
            return pyGridViewIterators[ dim ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "entities", [] ( const GridView &self, int codim ) {
            return pyGridViewIterators[ codim ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      // STATIC needed? - visibility needed?
      static const auto pyGridViewParIterators = registerPyGridViewParIterators< GridView >( cls, std::make_integer_sequence<int, dim+1>() );

      cls.def( "elements", [] ( const GridView &self, Dune::PartitionIteratorType p ) {
            return pyGridViewParIterators[ 0 ][ static_cast< std::size_t >( p ) ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "facets", [] ( const GridView &self, Dune::PartitionIteratorType p ) {
            return pyGridViewParIterators[ 1 ][ static_cast< std::size_t >( p ) ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "edges", [] ( const GridView &self, Dune::PartitionIteratorType p ) {
            return pyGridViewParIterators[ dim-1 ][ static_cast< std::size_t >( p ) ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "vertices", [] ( const GridView &self, Dune::PartitionIteratorType p ) {
            return pyGridViewParIterators[ dim ][ static_cast< std::size_t >( p ) ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "entities", [] ( const GridView &self, std::size_t codim, Dune::PartitionIteratorType p ) {
            return pyGridViewParIterators[ codim ][ static_cast< std::size_t >( p ) ]( self );
          }, pybind11::keep_alive< 0, 1 >() );

      // TODO
      auto entry = insertClass< PyIntersectionIterator<GridView> >( cls, "IntersectionIterator",
          GenerateTypeName("PyIntersectionIterator", cls),
          IncludeFiles{"dune/python/range.hh"} );
      if (entry.second)
        registerPyIterator< PyIntersectionIterator< GridView > >( cls, entry.first );
      cls.def( "intersections", [] ( const GridView &self, const typename GridView::template Codim< 0 >::Entity &e ) {
            return PyIntersectionIterator< GridView >( self.ibegin( e ), self.iend( e ) );
          }, pybind11::keep_alive< 0, 1 >() );

      cls.def("__repr__",
          [] (const GridView &gridView) -> std::string {
            return "LeafGrid with " + std::to_string(gridView.indexSet().size(0)) + " elements";
          });

      cls.def_property_readonly( "hierarchicalGrid", [] ( const GridView &self ) -> const Grid & { return self.grid(); } );

      cls.def_property_readonly( "dimension", [] ( const GridView &_ ) { return static_cast< int >( GridView::dimension ); } );
      cls.def_property_readonly( "dimensionworld", [] ( const GridView &_ ) { return static_cast< int >( GridView::dimensionworld ); } );

      cls.def_property_readonly( "conforming", [] ( const GridView &_ ) { return static_cast< bool >( GridView::conforming ); } );

      cls.def( "size", [] ( const GridView &self, int codim ) { return self.size( codim ); } );
      cls.def( "size", [] ( const GridView &self, Dune::GeometryType gt ) { return self.size( gt ); } );

      registerVTKWriter< GridView >( cls );
      cls.def( "vtkWriter", [] ( const GridView &self ) {
          return new VTKWriter< GridView >( self );
        }, pybind11::keep_alive< 0, 1 >() );
      cls.def( "vtkWriter", [] ( const GridView &self, int subsampling ) {
            return new SubsamplingVTKWriter< GridView >( self, subsampling );
          }, pybind11::keep_alive< 0, 1 >(), "subsampling"_a );

      cls.def("overlapSize", &GridView::overlapSize);
      cls.def("ghostSize", &GridView::ghostSize);

      cls.def_property_readonly("indexSet", &GridView::indexSet,
          pybind11::return_value_policy::reference_internal);

      cls.def_property_readonly( "comm", [] ( const GridView &gridView ) { return gridView.grid().comm(); } );

      cls.def( "communicate", [] ( const GridView &gridView, pybind11::object dataHandle, InterfaceType iftype, CommunicationDirection dir ) {
            ProxyDataHandle proxyDataHandle( std::move( dataHandle ) );
            gridView.communicate( proxyDataHandle, iftype, dir );
          });

      cls.def( "coordinates", [] ( const GridView &self ) { return coordinates( self ); } );
      cls.def( "tesselate", [] ( const GridView &self, int level ) { return tesselate( self, level ); }, "level"_a = 0 );
      cls.def( "polygons", [] ( const GridView &self ) { return polygons( self ); } );

      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &cls ] ( auto codim ) {
          cls.def( "contains", [] ( GridView &self, const typename GridView::template Codim< decltype( codim )::value >::Entity &entity ) -> bool {
              return self.contains( entity );
            }, "entity"_a );
        } );

      cls.def( "function", Dune::Python::defGridFunction< GridView >( cls, "GridFunction", std::make_integer_sequence< int, 11 >() ) );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_GRIDPART_HH
