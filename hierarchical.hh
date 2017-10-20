// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_HIERARCHICAL_HH
#define DUNE_PYTHON_GRID_HIERARCHICAL_HH

#include <array>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <type_traits>

#include <dune/common/iteratorrange.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/python/common/common.hh>

#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/factory.hh>
#include <dune/python/grid/gridview.hh>

#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dune
{

  namespace Python
  {

    // GridModificationListener
    // ------------------------

    template< class Grid >
    struct GridModificationListener
    {
      virtual ~GridModificationListener ()  {}

      virtual void preModification ( const Grid &grid ) {}
      virtual void postModification ( const Grid &grid ) {}
    };



    namespace detail
    {

      template< class Grid >
      using GridModificationListeners = std::multimap< const Grid *, GridModificationListener< Grid > * >;

      template< class Grid >
      inline GridModificationListeners< Grid > &gridModificationListeners ()
      {
        // TODO: add the listeners to dune.grid?
        static GridModificationListeners< Grid > listeners;
        return listeners;
      }

      template< class Grid >
      inline static IteratorRange< typename GridModificationListeners< Grid >::iterator > gridModificationListeners ( const Grid &grid )
      {
        typedef typename GridModificationListeners< Grid >::iterator Iterator;
        std::pair< Iterator, Iterator > range = gridModificationListeners< Grid >().equal_range( &grid );
        return IteratorRange< Iterator >( range.first, range.second );
      }

      template< class Grid >
      inline static void addGridModificationListener ( const Grid &grid, GridModificationListener< Grid > *listener, pybind11::handle nurse = pybind11::handle() )
      {
        // add listener to list of listeners
        auto pos = gridModificationListeners< Grid >().insert( std::make_pair( &grid, listener ) );

        // if nurse is given, ensure the listener dies with the nurse
        if( nurse )
        {
          pybind11::cpp_function remove( [ pos ] ( pybind11::handle weakref ) {
              delete pos->second;
              gridModificationListeners< Grid >().erase( pos );
              weakref.dec_ref();
            } );
          pybind11::weakref( nurse, remove ).release();
        }
      }

    } // namespace detail



    // registerHierarchicalGridPicklingSupport
    // ---------------------------------------

    template< class Grid, class... options >
    inline static std::enable_if_t< Capabilities::hasBackupRestoreFacilities< Grid >::v >
    registerHierarchicalGridPicklingSupport ( pybind11::class_< Grid, options... > cls, PriorityTag< 1 > )
    {
      cls.def( pybind11::pickle( [] ( const Grid &self ) -> pybind11::bytes {
          std::ostringstream stream;
          BackupRestoreFacility< Grid >::backup( self, stream );
          return stream.str();
        }, [] ( pybind11::bytes state ) -> std::unique_ptr< Grid > {
          std::istringstream stream( state );
          return std::unique_ptr< Grid >( BackupRestoreFacility< Grid >::restore( stream ) );
        } ) );
    }

    template< class Grid, class... options >
    inline static void
    registerHierarchicalGridPicklingSupport ( pybind11::class_< Grid, options... > cls, PriorityTag< 0 > )
    {}

    template< class Grid, class... options >
    inline static void registerHierarchicalGridPicklingSupport ( pybind11::class_< Grid, options... > cls )
    {
      registerHierarchicalGridPicklingSupport( cls, PriorityTag< 42 >() );
    }



    // readDGF
    // -------

    template< class Grid >
    inline static std::unique_ptr< Grid > readDGF ( const std::string &fileName )
    {
      DGFGridFactory< Grid > dgfFactory( fileName );
      std::unique_ptr< Grid > grid( dgfFactory.grid() );
      grid->loadBalance();
      return grid;
    }

    template< class Grid >
    inline static std::unique_ptr< Grid > readDGF ( std::istream &input )
    {
      DGFGridFactory< Grid > dgfFactory( input );
      std::unique_ptr< Grid > grid( dgfFactory.grid() );
      grid->loadBalance();
      return grid;
    }



    // readGmsh
    // --------

    template< class Grid, std::enable_if_t< Capabilities::HasGridFactory< Grid >::value, int > = 0 >
    inline static std::unique_ptr< Grid > readGmsh ( const std::string &fileName )
    {
      Dune::GridFactory< Grid > gridFactory;
      Dune::GmshReader< Grid >::read( gridFactory, fileName, false, false );
      return std::unique_ptr< Grid >( gridFactory.createGrid() );
    }

    template< class Grid, std::enable_if_t< !Capabilities::HasGridFactory< Grid >::value, int > = 0 >
    inline static std::unique_ptr< Grid > readGmsh ( const std::string &fileName )
    {
      throw std::invalid_argument( "Can only read Gmsh files into grids supporting the GridFactory concept." );
    }



    // reader
    // ------

    template< class Grid >
    inline static std::unique_ptr< Grid > reader ( const std::tuple< Reader, std::string > &args )
    {
      switch( std::get< 0 >( args ) )
      {
        case Reader::dgf:
          return readDGF< Grid >( std::get< 1 >( args ) );

        case Reader::dgfString:
          {
            std::istringstream input( std::get< 1 >( args ) );
            return readDGF< Grid >( input );
          }

        case Reader::gmsh:
          return readGmsh< Grid >( std::get< 1 >( args ) );

        default:
          return nullptr;
      }
    }

    template< typename Grid >
    using StructuredReader = std::tuple<
        Reader,
        std::string,
        FieldVector< typename Grid::ctype, Grid::dimensionworld >,
        FieldVector< typename Grid::ctype, Grid::dimensionworld >,
        std::array< unsigned int, Grid::dimension >
      >;

    template< class Grid >
    inline static auto reader ( const StructuredReader< Grid > &args, PriorityTag< 1 > )
      -> std::enable_if_t< Capabilities::HasStructuredGridFactory< Grid >::value, std::unique_ptr< Grid > >
    {
      using std::get;

      if( get< 0 >( args ) != Reader::structured )
        throw pybind11::value_error( "This overload only supports the 'structured' grid reader." );

      if( get< 1 >( args ) == "cube" )
        return StructuredGridFactory< Grid >::createCubeGrid( get< 2 >( args ), get< 3 >( args ), get< 4 >( args ) );
      else if( get< 1 >( args ) == "simplex" )
        return StructuredGridFactory< Grid >::createSimplexGrid( get< 2 >( args ), get< 3 >( args ), get< 4 >( args ) );
      else
        throw pybind11::value_error( "Unknown grid element type." );
    }

    template< class Grid >
    inline static std::unique_ptr< Grid > reader ( const StructuredReader< Grid > &args, PriorityTag< 0 > )
    {
      throw pybind11::value_error( "Can only create structured grids for grids supporting the GridFactory concept." );
    }

    template< class Grid >
    inline static std::unique_ptr< Grid > reader ( const StructuredReader< Grid > &args )
    {
      return reader< Grid >( args, PriorityTag< 42 >() );
    }

    template< class Grid >
    inline static auto reader ( const pybind11::dict &args, PriorityTag< 1 > )
      -> std::enable_if_t< Capabilities::HasGridFactory< Grid >::value, std::unique_ptr< Grid > >
    {
      GridFactory< Grid > factory;
      fillGridFactory( args, factory );
      return std::unique_ptr< Grid >( factory.createGrid() );
    }

    template< class Grid >
    inline static std::unique_ptr< Grid > reader ( const pybind11::dict &args, PriorityTag< 0 > )
    {
      throw pybind11::value_error( "Can only read Python dictionaries into grids supporting the GridFactory concept." );
    }

    template< class Grid >
    inline static std::unique_ptr< Grid > reader ( const pybind11::dict &args )
    {
      return reader< Grid >( args, PriorityTag< 42 >() );
    }



    // registerHierarchicalGrid
    // ------------------------

    template< class Grid, class... options >
    void registerHierarchicalGrid ( pybind11::module module, pybind11::class_< Grid, options... > cls )
    {
      pybind11::module::import( "dune.geometry" );

      auto clsLevelView =
        insertClass<typename Grid::LevelGridView>
        ( module, "LevelGrid", GenerateTypeName(Dune::MetaType<Grid>(), "LevelGridView") );
      if (clsLevelView.second)
        registerGridView( module, clsLevelView.first );

      auto clsLeafView = insertClass<typename Grid::LeafGridView>
        ( module, "LeafGrid", GenerateTypeName(Dune::MetaType<Grid>(), "LeafGridView") );
      if (clsLeafView.second)
        registerGridView( module, clsLeafView.first );

      module.def( "reader", [] ( const std::tuple< Reader, std::string > &args ) { return reader< Grid >( args ); } );
      module.def( "reader", [] ( const std::string &args ) { return reader< Grid >( std::make_tuple( Reader::dgf,args ) ); } );
      module.def( "reader", [] ( const StructuredReader<Grid> &args ) { return reader< Grid >( args ); } );
      module.def( "reader", [] ( const pybind11::dict &args ) { return reader< Grid >( args ); } );

      registerHierarchicalGridPicklingSupport( cls );

      cls.def_property_readonly( "leafView", pybind11::cpp_function( [] ( const Grid &self ) { return self.leafGridView(); }, pybind11::keep_alive< 0, 1 >() ) );
      cls.def( "levelView", [] ( const Grid &self, int level ) { return self.levelGridView( level ); }, pybind11::keep_alive< 0, 1 >() );

      typedef typename Grid::template Codim< 0 >::Entity Element;

      cls.def( "mark", [] ( Grid &grid, const std::function< Marker( const Element &e ) > &marking ) {
          std::pair< int, int > marked;
          for( const Element &element : elements( grid.leafGridView() ) )
          {
            Marker marker = marking( element );
            marked.first += static_cast< int >( marker == Marker::Refine );
            marked.second += static_cast< int >( marker == Marker::Coarsen );
            grid.mark( static_cast< int >( marker ), element );
          }
          return marked;
        } );

      cls.def( "adapt", [] ( Grid &grid ) {
          for( const auto &listener : detail::gridModificationListeners( grid ) )
            listener.second->preModification( grid );
          grid.preAdapt();
          grid.adapt();
          grid.postAdapt();
          for( const auto &listener : detail::gridModificationListeners( grid ) )
            listener.second->postModification( grid );
        } );

      cls.def( "globalRefine", [] ( Grid &grid, int level ) {
          for( const auto &listener : detail::gridModificationListeners( grid ) )
            listener.second->preModification( grid );
          grid.globalRefine( level );
          for( const auto &listener : detail::gridModificationListeners( grid ) )
            listener.second->postModification( grid );
        } );

      cls.def( "loadBalance", [] ( Grid &grid ) {
          for( const auto &listener : detail::gridModificationListeners( grid ) )
            listener.second->preModification( grid );
          grid.loadBalance();
          for( const auto &listener : detail::gridModificationListeners( grid ) )
            listener.second->postModification( grid );
        } );

      cls.def_property_readonly( "refineStepsForHalf", [] ( const Grid & ) { return DGFGridInfo< Grid >::refineStepsForHalf(); } );

      cls.def_property_readonly( "comm", [] ( const Grid &grid ) { return grid.comm(); } );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_HIERARCHICAL_HH
