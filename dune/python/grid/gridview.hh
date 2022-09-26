// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
#include <dune/python/grid/indexset.hh>
#include <dune/python/grid/mapper.hh>
#include <dune/python/grid/intersection.hh>
#include <dune/python/grid/numpy.hh>
#include <dune/python/grid/range.hh>
#include <dune/python/grid/vtk.hh>
#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/numpycommdatahandle.hh>

#include <dune/python/pybind11/pybind11.h>

#include <iostream>

#if HAVE_DUNE_VTK
#include <dune/vtk/function.hh>
#endif

namespace Dune
{

  namespace Python
  {

    // ProxyDataHandle for parallel communication
    // ------------------------------------------

    struct DUNE_PRIVATE ProxyDataHandle
      : public Dune::CommDataHandleIF< ProxyDataHandle, double >
    {
      ProxyDataHandle ( pybind11::object dataHandle )
        : contains_( method( dataHandle, "contains" ) ),
          fixedSze_( method( dataHandle, "fixedSize" ) ),
          size_( method( dataHandle, "size" ) ),
          gather_( method( dataHandle, "gather" ) ),
          scatter_( method( dataHandle, "scatter" ) )
      {}

      bool contains ( int dim, int codim ) const { return contains_( dim, codim ).cast< bool >(); }
      bool fixedSize ( int dim, int codim ) const { return fixedSze_( dim, codim ).cast< bool >(); }

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

      pybind11::object contains_, fixedSze_;
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
    inline static void registerGridView ( pybind11::handle scope, pybind11::class_< GridView, options... > cls )
    {
      typedef typename GridView::Grid Grid;
      typedef PyGridViewIterator< GridView, 0 > PyElementIterator;

      using pybind11::operator""_a;

      pybind11::options opts;
      opts.disable_function_signatures();

      detail::registerGridViewConstructorFromGrid( cls, PriorityTag< 42 >() );

      cls.attr( "dimGrid" ) = pybind11::int_( static_cast< int >( GridView::dimension ) );
      cls.attr( "dimWorld" ) = pybind11::int_( static_cast< int >( GridView::dimensionworld ) );

      registerGridEntities< GridView >( cls );
      registerGridIntersection< GridView >( cls );

      cls.def( "_mapper", [] ( GridView &self, pybind11::object layout ) {
          return makeMultipleCodimMultipleGeomTypeMapper( self, layout );
        }, pybind11::keep_alive< 0, 1 >(), "layout"_a,
        R"doc(
          Set up a mapper to attach data to the grid. The layout argument defines how many
          degrees of freedom to assign to each subentity of a geometry type.

          Args:
              layout:     function, dict, tuple, or list defining the number of indices to reserve
                          for each geometry type.

          If layout is a dict, is must map geometry types to integers. All types not mentioned in
          the dictionary are assumed to be zero.

          If layout is a tuple or a list, it must contain exactly dimension+1 integers, one for
          each codimension in the grid.

          if layout is a function it maps a geometry type to the number of degrees of freedom to
          reserve. Here a return value of 0 or `False` indicates that no data is to be attach, `True` can be used instead of 1.

          Returns:   the mapper
        )doc" );

      // register iterators
      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [] ( auto codim ) {
          registerPyGridViewIterator< GridView, codim >();
        } );

      if( Capabilities::canIterate< Grid, 0 >::value )
        cls.def_property_readonly( "elements", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, 0 >( self ); },
          R"doc(
            sequence of all elements (i.e., entities of codimension 0)
          )doc" );
      if( Capabilities::canIterate< Grid, 1 >::value )
        cls.def_property_readonly( "facets", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, 1 >( self ); },
          R"doc(
            range of all facets (i.e., entities of codimension 1)
          )doc" );
      if( Capabilities::canIterate< Grid, GridView::dimension-1 >::value )
        cls.def_property_readonly( "edges", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, GridView::dimension-1 >( self ); },
          R"doc(
            range of all edges (i.e., entities of dimension 1)
          )doc" );
      if( Capabilities::canIterate< Grid, GridView::dimension >::value )
        cls.def_property_readonly( "vertices", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, GridView::dimension >( self ); },
          R"doc(
            range of all vertices (i.e., entities of dimension 0)
          )doc" );

      std::array< pybind11::object (*) ( pybind11::object ), GridView::dimension+1 > makePyGridViewIterators;
      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &makePyGridViewIterators ] ( auto codim ) {
          makePyGridViewIterators[ codim ] = makePyGridViewIterator< GridView, codim >;
        } );
      cls.def( "entities", [ makePyGridViewIterators ] ( pybind11::object self, int codim ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          return makePyGridViewIterators[ codim ]( self );
        }, "codim"_a,
        R"doc(
          get range of entities for a codimension

          Args:
              codim:    Codimension to obtain range of entities for
        )doc" );

      registerPyIntersectionIterator< GridView >();
      cls.def( "intersections", [] ( const GridView &self, const typename GridView::template Codim< 0 >::Entity &e ) {
          return PyIntersectionIterator< GridView >( self.ibegin( e ), self.iend( e ) );
        }, pybind11::keep_alive< 0, 1 >(), "element"_a,
        R"doc(
          get range of all codim-1 intersections for an element

          Args:
              element:    Element of obtain intersections for
        )doc" );

      registerPyBoundaryIntersectionIterator< GridView, PyElementIterator >();
      cls.def_property_readonly( "boundaryIntersections", [] ( const GridView &self ) {
          return PyBoundaryIntersectionIterator< GridView, PyElementIterator >( self, PyElementIterator( self.template begin< 0 >(), self.template end< 0 >() ) );
        }, pybind11::keep_alive< 0, 1 >(),
        R"doc(
          range of all codim-1 boundary intersections of the grid
        )doc" );

      // register partitions

      registerGridViewPartition< GridView, Interior_Partition >();
      registerGridViewPartition< GridView, InteriorBorder_Partition >();
      registerGridViewPartition< GridView, Overlap_Partition >();
      registerGridViewPartition< GridView, OverlapFront_Partition >();
      registerGridViewPartition< GridView, All_Partition >();
      registerGridViewPartition< GridView, Ghost_Partition >();

      cls.def_property_readonly( "interiorPartition", [] ( pybind11::object self ) {
          return GridViewPartition< GridView, Interior_Partition >( self );
        }, R"doc(
          Partition containing only interior entities.
        )doc" );
      cls.def_property_readonly( "interiorBorderPartition", [] ( pybind11::object self ) {
          return GridViewPartition< GridView, InteriorBorder_Partition >( self );
        }, R"doc(
          Partition containing only interior and border entities.
        )doc" );
      cls.def_property_readonly( "overlapPartition", [] ( pybind11::object self ) {
          return GridViewPartition< GridView, Overlap_Partition >( self );
        }, R"doc(
          Partition containing only interior, border and overlap entities.
        )doc" );
      cls.def_property_readonly( "overlapFrontPartition", [] ( pybind11::object self ) {
          return GridViewPartition< GridView, OverlapFront_Partition >( self );
        }, R"doc(
          Partition containing only interior, border, overlap, and front entities.
        )doc" );
      cls.def_property_readonly( "allPartition", [] ( pybind11::object self ) {
          return GridViewPartition< GridView, All_Partition >( self );
        }, R"doc(
          Partition containing all entities.
        )doc" );
      cls.def_property_readonly( "ghostPartition", [] ( pybind11::object self ) {
          return GridViewPartition< GridView, Ghost_Partition >( self );
        }, R"doc(
          Partition containing only ghost entities.
        )doc" );

      cls.def("__repr__",
          [] (const GridView &gridView) -> std::string {
            return "LeafGrid with " + std::to_string(gridView.indexSet().size(0)) + " elements";
          });

      cls.def_property_readonly( "hierarchicalGrid", [] ( const GridView &self ) -> const Grid & { return self.grid(); },
        R"doc(
          associated hierarchical grid
        )doc" );

      cls.def_property_readonly_static( "dimension", [] ( pybind11::object ) { return int(GridView::dimension); } );
      cls.def_property_readonly_static( "dimensionworld", [] ( pybind11::object ) { return int(GridView::dimensionworld); } );

      cls.def( "size", [] ( const GridView &self, int codim ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          return self.size( codim );
        }, "codim"_a,
        R"doc(
          Args:
              codim:     required codimension
          Returns:       number of subentities of codimension `codim`
        )doc" );
      cls.def( "size", [] ( const GridView &self, Dune::GeometryType gt ) {
          if( (gt.dim() < 0) || (gt.dim() > GridView::dimension) )
            throw pybind11::value_error( "Invalid geometry type (dimension must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          return self.size( gt );
        }, "gt"_a,
        R"doc(
          Args:
              gt:        a geometry type
          Returns:       number of subentities of the given geometry type
        )doc" );

      registerVTKWriter< GridView >( cls );
      cls.def( "vtkWriter", [] ( const GridView &self, const bool nonconforming = false ) {
          const VTK::DataMode dm = nonconforming ? VTK::nonconforming : VTK::conforming;
          return new VTKWriter< GridView >( self, dm );
        }, pybind11::keep_alive< 0, 1 >() );
      cls.def( "vtkWriter", [] ( const GridView &self, int subsampling ) {
            return new SubsamplingVTKWriter< GridView >( self,
                  Dune::refinementIntervals(1<<subsampling) );
          }, pybind11::keep_alive< 0, 1 >(), "subsampling"_a );

      cls.def( "overlapSize", [] ( const GridView &self, int codim ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          return self.overlapSize( codim );
        }, "codim"_a );
      cls.def( "ghostSize", [] ( const GridView &self, int codim ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          return self.ghostSize( codim );
        }, "codim"_a );

      cls.def_property_readonly( "_indexSet", [] ( const GridView &self ) -> const typename GridView::IndexSet & {
          return self.indexSet();
        }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
        R"doc(
          index set for the grid
        )doc" );

      cls.def_property_readonly( "comm", [] ( const GridView &gridView ) -> const typename Grid::CollectiveCommunication & {
          return gridView.grid().comm();
        }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
        R"doc(
          collective communication for the grid

          Note: For collective (or global) operations, all processes in this
                collective communication must call the corresponding method.
        )doc" );

      typedef NumPyCommDataHandle< MultipleCodimMultipleGeomTypeMapper< GridView >, double, std::function< double ( double, double ) > > CommDataHandle;
      cls.def( "communicate", [] ( const GridView &gridView, CommDataHandle &dataHandle, InterfaceType iftype, CommunicationDirection dir ) {
            gridView.communicate( dataHandle, iftype, dir );
          } );
      cls.def( "communicate", [] ( const GridView &gridView, pybind11::object dataHandle, InterfaceType iftype, CommunicationDirection dir ) {
            ProxyDataHandle proxyDataHandle( std::move( dataHandle ) );
            gridView.communicate( proxyDataHandle, iftype, dir );
          });

      // export grid capabilities

      cls.def_property_readonly( "conforming", [] ( pybind11::object ) { return static_cast< bool >( GridView::conforming ); } );

      if( Capabilities::hasSingleGeometryType< Grid >::v )
        cls.def_property_readonly_static( "type", [] ( pybind11::object ) {
            return GeometryType( Capabilities::hasSingleGeometryType< Grid >::topologyId, Grid::dimension );
          } );

      cls.def_property_readonly_static( "isCartesian", [] ( pybind11::object ) { return Capabilities::isCartesian< Grid >::v; } );
      cls.def_property_readonly_static( "canCommunicate", [] ( pybind11::object ) {
          pybind11::tuple canCommunicate( Grid::dimension+1 );
          Hybrid::forEach( std::make_integer_sequence< int, Grid::dimension+1 >(), [ &canCommunicate ] ( auto codim ) {
              canCommunicate[ codim ] = pybind11::cast( bool( Capabilities::canCommunicate< Grid, codim >::v ) );
            } );
          return canCommunicate;
        } );

      cls.def_property_readonly_static( "threadSafe", [] ( pybind11::object ) { return Capabilities::viewThreadSafe< Grid >::v; } );

      // export utility methods

      cls.def( "coordinates", [] ( const GridView &self ) { return coordinates( self ); },
        R"doc(
          Returns: `numpy` array with the coordinates of all vertices in the grid in
                   the format `[ [x_1,y_1], [x_2,y_2], ..., [x_N,y_N] ]` for example
                   in 2d.
        )doc" );
      cls.def( "tessellate", [] ( const GridView &self, int level ) { return tessellate( self, level ); }, "level"_a = 0,
        R"doc(
          Generated a possibly refined tessellation using only simplices.

          Args:
              level: virtual refinement level to use to generate the tessellation

          Returns: (coordinates,simplices) where coordinates is a `numpy` array
                   of the vertex coodinates
                   (e.g. in 2d `[ [x_1,y_1], [x_2,y_2], ..., [x_N,y_N] ]` )
                   and simplices is a `numpy` array of the vertices of the simplices
                   (e.g. in 2d `[s_11,s_12,s_13], [s_21,s_22,s_23], ..., [s_N1,s_N2,s_N3] ]` )

        )doc" );
      cls.def( "polygons", [] ( const GridView &self ) { return polygons( self ); },
        R"doc(
          Store the grid in numpy arrays.

          Returns: coordinate array storing the vertex coordinate of each polygon
                   in the grid.
        )doc" );

      cls.def( "contains", [] ( GridView &self, pybind11::object entity ) {
          bool found = false, contained = false;
          Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &self, entity, &found, &contained ] ( auto codim ) {
              typedef typename GridView::template Codim< decltype( codim )::value >::Entity Entity;
              if( pybind11::isinstance< Entity >( entity ) )
              {
                found = true;
                contained = self.contains( pybind11::cast< const Entity & >( entity ) );
              }
            } );
          if( found )
            return contained;
          else
            throw pybind11::value_error( "Argument 'entity' is not a valid entity for this grid." );
        }, "entity"_a,
        R"doc(
          Check whether an entity is contained in the grid instance

          Args:
              entity:   entity to check

          Note:
          - The entity must be contained in the corresponding hierarchical grid.
        )doc" );

#if HAVE_DUNE_VTK
      using VirtualizedGF = Dune::Vtk::Function<GridView>;
      auto vgfClass = Python::insertClass<VirtualizedGF>(scope,"VtkFunction",
          Python::GenerateTypeName("Dune::Vtk::Function", MetaType<GridView>()),
          Python::IncludeFiles{"dune/vtk/function.hh"});
      if( vgfClass.second )
      {
        vgfClass.first.def("name",[](VirtualizedGF &self) { return self.name(); });
      }
#endif
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_GRIDPART_HH
