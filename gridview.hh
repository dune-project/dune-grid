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

#include <dune/python/common/logger.hh>
#include <dune/python/grid/entity.hh>
#include <dune/python/grid/function.hh>
#include <dune/python/grid/indexset.hh>
#include <dune/python/grid/intersection.hh>
#include <dune/python/grid/mapper.hh>
#include <dune/python/grid/numpy.hh>
#include <dune/python/grid/range.hh>
#include <dune/python/grid/vtk.hh>
#include <dune/python/grid/capabilities.hh>
#include <dune/python/utility/numpycommdatahandle.hh>

#include <dune/python/pybind11/pybind11.h>

#include <iostream>

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
    inline static void registerGridView ( pybind11::handle scope, pybind11::class_< GridView, options... > cls )
    {
      using pybind11::operator""_a;

      typedef typename GridView::Grid Grid;
      typedef PyGridViewIterator< GridView, 0 > PyElementIterator;

      detail::registerGridViewConstructorFromGrid( cls, PriorityTag< 42 >() );

      cls.attr( "dimGrid" ) = pybind11::int_( static_cast< int >( GridView::dimension ) );
      cls.attr( "dimWorld" ) = pybind11::int_( static_cast< int >( GridView::dimensionworld ) );

      registerGridEntities< GridView >( cls );
      registerGridIntersection< GridView >( cls );

      registerGridViewIndexSet< GridView >( cls );
      registerMultipleCodimMultipleGeomTypeMapper< GridView >( cls );

      typedef MultipleCodimMultipleGeomTypeMapper< GridView > MCMGMapper;
      cls.def( "mapper", [] ( GridView &self, pybind11::function layout ) {
          return new MCMGMapper( self, [ layout ] ( Dune::GeometryType gt, int griddim ) { return static_cast< unsigned int >( pybind11::cast< int >( layout( gt ) )); } );
        }, pybind11::keep_alive< 0, 1 >(), "layout"_a,
        R"doc(
          Set up a mapper to attach data to the given grid. The layout argument
          defines the set of subentities (given by geometry type)
          to attach the data to and how many degrees of freedom to use per subentity.

          Args:
              layout:     function taken a geometry type and returning the number
                          of dof to attach to the entities of that geometry type.
                          0 or `False`: do not attach any, `True` can be used instead of 1.

          Returns:   the mapper
        )doc" );

      // register iterators

      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &cls ] ( auto &&codim ) {
          registerPyGridViewIterator< GridView, codim >();
        } );

      if( Capabilities::canIterate< Grid, 0 >::value )
        cls.def_property_readonly( "elements", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, 0 >( self ); } );
      if( Capabilities::canIterate< Grid, 1 >::value )
        cls.def_property_readonly( "facets", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, 1 >( self ); } );
      if( Capabilities::canIterate< Grid, GridView::dimension-1 >::value )
        cls.def_property_readonly( "edges", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, GridView::dimension-1 >( self ); } );
      if( Capabilities::canIterate< Grid, GridView::dimension >::value )
        cls.def_property_readonly( "vertices", [] ( pybind11::object self ) { return makePyGridViewIterator< GridView, GridView::dimension >( self ); } );

      std::array< pybind11::object (*) ( pybind11::object ), GridView::dimension+1 > makePyGridViewIterators;
      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &makePyGridViewIterators ] ( auto &&codim ) {
          makePyGridViewIterators[ codim ] = makePyGridViewIterator< GridView, codim >;
        } );
      cls.def( "entities", [ makePyGridViewIterators ] ( pybind11::object self, int codim ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          return makePyGridViewIterators[ codim ]( self );
        }, "codim"_a );
      Logger logger( "dune.grid" );
      cls.def( "entities", [ logger ] ( pybind11::object self, int codim, PartitionIteratorType pitype ) {
          if( (codim < 0) || (codim > GridView::dimension) )
            throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
          logger.error( "The entities no longer has a partition type argument." );
          switch( pitype )
          {
          case Interior_Partition:
            logger.error( "The parallel versions can be obtained from the corresponding partition, i.e., interiorPartition." );
            break;

          case InteriorBorder_Partition:
            logger.error( "The parallel versions can be obtained from the corresponding partition, i.e., interiorBorderPartition." );
            break;

          case Overlap_Partition:
            logger.error( "The parallel versions can be obtained from the corresponding partition, i.e., overlapPartition." );
            break;

          case OverlapFront_Partition:
            logger.error( "The parallel versions can be obtained from the corresponding partition, i.e., overlapFrontPartition." );
            break;

          case All_Partition:
            logger.error( "The parallel versions can be obtained from the corresponding partition, i.e., allPartition." );
            break;

          case Ghost_Partition:
            logger.error( "The parallel versions can be obtained from the corresponding partition, i.e., ghostPartition." );
            break;
          }
          throw pybind11::value_error( "The entities no longer has a partition type argument." );
        } );

      registerPyIntersectionIterator< GridView >();
      cls.def( "intersections", [] ( const GridView &self, const typename GridView::template Codim< 0 >::Entity &e ) {
          return PyIntersectionIterator< GridView >( self.ibegin( e ), self.iend( e ) );
        }, pybind11::keep_alive< 0, 1 >() );

      registerPyBoundaryIntersectionIterator< GridView, PyElementIterator >();
      cls.def_property_readonly( "boundaryIntersections", [] ( const GridView &self ) {
          return PyBoundaryIntersectionIterator< GridView, PyElementIterator >( self, PyElementIterator( self.template begin< 0 >(), self.template end< 0 >() ) );
        }, pybind11::keep_alive< 0, 1 >() );

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

      cls.def_property_readonly( "hierarchicalGrid", [] ( const GridView &self ) -> const Grid & { return self.grid(); } );

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
      cls.def( "vtkWriter", [] ( const GridView &self ) {
          return new VTKWriter< GridView >( self );
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

      cls.def_property_readonly("indexSet", &GridView::indexSet,
          pybind11::return_value_policy::reference_internal);

      cls.def_property_readonly( "comm", [] ( const GridView &gridView ) { return gridView.grid().comm(); } );

      cls.def( "communicate", [] ( const GridView &gridView,
                                   NumPyCommDataHandle<MCMGMapper,double,std::function<double(double,double)>> &dataHandle, InterfaceType iftype, CommunicationDirection dir ) {
            gridView.communicate( dataHandle, iftype, dir );
          });
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
          Hybrid::forEach( std::make_integer_sequence< int, Grid::dimension+1 >(), [ &canCommunicate ] ( auto &&codim ) {
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
      cls.def( "tesselate", [] ( const GridView &self, int level ) { return tesselate( self, level ); }, "level"_a = 0,
        R"doc(
          Generated a possibly refined tesselation using only simplices.

          Args:
              level: virtual refinement level to use to generate the tesselation

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
