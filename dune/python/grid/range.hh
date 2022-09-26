// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_RANGE_HH
#define DUNE_PYTHON_GRID_RANGE_HH

#include <string>
#include <utility>

#include <dune/common/visibility.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/python/grid/capabilities.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/common/logger.hh>
#include <dune/python/pybind11/extensions.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    // PyIterator
    // ----------

    template< class Iterator, class Entity >
    struct PyIterator
    {
      PyIterator ( Iterator begin, Iterator end )
        : it_( std::move( begin ) ), end_( std::move( end ) )
      {}

      Entity next ()
      {
        if( it_ == end_ )
          throw pybind11::stop_iteration();

        Entity entity = *it_;
        ++it_;
        return entity;
      }

    private:
      Iterator it_, end_;
    };



    // registerPyIterator
    // ------------------

    template< class Iterator >
    void registerPyIterator( pybind11::handle scope, pybind11::class_< Iterator > cls )
    {
      cls.def( "__iter__", [] ( Iterator &self ) -> Iterator & { return self; } );
      cls.def( "__next__", [] ( Iterator &self ) { return self.next(); } );
    }



    // PyGridViewIterator
    // ------------------

    template< class GridView, int codim >
    using PyGridViewIterator = PyIterator< typename GridView::template Codim< codim >::Iterator, typename GridView::template Codim< codim >::Entity >;



    // PyGridViewPartitionIterator
    // ---------------------------

    template< class GridView, int codim, PartitionIteratorType partition >
    using PyGridViewPartitionIterator = PyIterator< typename GridView::template Codim< codim >::template Partition< partition >::Iterator, typename GridView::template Codim< codim >::Entity >;



    // PyIntersectionIterator
    // ----------------------

    template< class GridView >
    using PyIntersectionIterator = PyIterator< typename GridView::IntersectionIterator, typename GridView::Intersection >;



    // PyBoundaryIntersectionIterator
    // ------------------------------

    template< class GridView, class PyElementIterator >
    struct PyBoundaryIntersectionIterator
    {
      typedef typename GridView::Intersection Intersection;

      PyBoundaryIntersectionIterator ( const GridView &gridView, PyElementIterator it )
        : gridView_( gridView ), elementIt_( std::move( it ) )
      {}

      Intersection next ()
      {
        while( true )
        {
          if( intersectionIt_ != intersectionEnd_ )
          {
            Intersection intersection = *intersectionIt_;
            ++intersectionIt_;
            if( intersection.boundary() )
              return intersection;
          }
          else
          {
            auto element = elementIt_.next();
            while( !element.hasBoundaryIntersections() )
              element = elementIt_.next();

            intersectionIt_ = gridView_.ibegin( element );
            intersectionEnd_ = gridView_.iend( element );
          }
        }
      }

    private:
      const GridView &gridView_;
      PyElementIterator elementIt_;
      typename GridView::IntersectionIterator intersectionIt_, intersectionEnd_;
    };



#if 0
    // PyGridViewPartitionIntersectionIterator
    // ---------------------------------------

    template< class GridView, PartitionIteratorType partition >
    struct PyGridViewPartitionIntersectionIterator
    {
      typedef typename GridView::Intersection Intersection;

      PyGridViewPartitionIntersectionIterator ( const GridView &gridView )
        : gridView_( gridView ), mapper_( gridView_, mcmgElementLayout() ),
          elementIt_( gridView_.template begin< 0, partition >(), gridView_.template end< 0, partition >() )
      {}

      Intersection next ()
      {
        while( true )
        {
          if( intersectionIt_ != intersectionEnd_ )
          {
            Intersection intersection = *intersectionIt_;
            ++intersectionIt_;
            if( !intersection.neighbor() )
              return intersection;

            auto outside = intersection.outside();
            if( !partitionSet< partition >().contains( outside.partitionType() ) || (insideIndex_ < mapper_.index( outside )) )
              return intersection;
          }
          else
          {
            auto element = elementIt_.next();
            insideIndex_ = mapper_.index( element );
            intersectionIt_ = gridView_.ibegin( element );
            intersectionEnd_ = gridView_.iend( element );
          }
        }
      }

    private:
      typedef MultipleCodimMultipleGeomTypeMapper< GridView > Mapper;

      const GridView &gridView_;
      Mapper mapper_;
      PyGridViewPartitionIterator< GridView, 0, partition > elementIt_;
      typename Mapper::Index insideIndex_;
      typename GridView::IntersectionIterator intersectionIt_, intersectionEnd_;
    };
#endif




    // registerPyGridViewIterator
    // --------------------------

    template< class GridView, int codim >
    inline static auto registerPyGridViewIterator ( pybind11::handle scope, PriorityTag< 1 > )
      -> std::enable_if_t< Capabilities::canIterate< typename GridView::Grid, codim >::value >
    {
      typedef PyGridViewIterator< GridView, codim > Iterator;
      auto typeName = GenerateTypeName( "PyGridViewIterator", MetaType< GridView >(), codim );
      auto entry = insertClass< Iterator >( scope, "EntityIterator", typeName, IncludeFiles{ "dune/python/grid/range.hh" } );
      if( entry.second )
      {
        registerPyIterator( scope, entry.first );
        Logger logger( "dune.grid" );
        entry.first.def( "__call__", [ logger ] ( pybind11::object self ) {
            logger.warning( "The methods elements, facets, edges, and vertices have been converted to properties." );
            logger.warning( "Please remove the trailing parentesis." );
            return self;
          } );
        entry.first.def( "__call__", [ logger ] ( pybind11::object self, PartitionIteratorType pitype ) {
            logger.error( "The methods elements, facets, edges, and vertices have been converted to properties." );
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
            throw pybind11::value_error( "The methods elements, facets, edges, and vertices have been converted to properties." );
          } );
      }
    }

    template< class GridView, int codim >
    inline static void registerPyGridViewIterator ( pybind11::handle scope, PriorityTag< 0 > )
    {}

    template< class GridView, int codim >
    inline static void registerPyGridViewIterator ( pybind11::handle scope = {} )
    {
      return registerPyGridViewIterator< GridView, codim >( scope, PriorityTag< 42 >() );
    }



    // makePyGridViewIterator
    // ----------------------

    template< class GridView, int codim >
    inline static auto makePyGridViewIterator ( pybind11::object obj )
      -> std::enable_if_t< Capabilities::canIterate< typename GridView::Grid, codim >::value, pybind11::object >
    {
      const GridView &gridView = pybind11::cast< const GridView & >( obj );
      pybind11::object iterator = pybind11::cast( new PyGridViewIterator< GridView, codim >( gridView.template begin< codim >(), gridView.template end< codim >() ) );
      pybind11::detail::keep_alive_impl( iterator, obj );
      return iterator;
    }

    template< class GridView, int codim >
    inline static auto makePyGridViewIterator ( pybind11::object )
      -> std::enable_if_t< !Capabilities::canIterate< typename GridView::Grid, codim >::value, pybind11::object >
    {
      throw pybind11::value_error( "Iterators for codimension " + std::to_string( codim ) + " are not implemented." );
    }



    // registerPyIntersectionIterator
    // ------------------------------

    template< class GridView >
    inline static void registerPyIntersectionIterator ( pybind11::handle scope = {} )
    {
      typedef PyIntersectionIterator< GridView > Iterator;
      auto typeName = GenerateTypeName( "Dune::Python::PyIntersectionIterator", MetaType< GridView >() );
      auto entry = insertClass< Iterator >( scope, "IntersectionIterator", typeName, IncludeFiles{ "dune/python/range.hh" } );
      if( entry.second )
        registerPyIterator< Iterator >( scope, entry.first );
    }



    // registerPyBoundaryIntersectionIterator
    // --------------------------------------

    template< class GridView, class PyElementIterator >
    inline static void registerPyBoundaryIntersectionIterator ( pybind11::handle scope = {} )
    {
      typedef PyBoundaryIntersectionIterator< GridView, PyElementIterator > Iterator;
      auto typeName = GenerateTypeName( "Dune::Python::PyBoundaryIntersectionIterator", MetaType< GridView >(), MetaType< PyElementIterator >() );
      auto entry = insertClass< Iterator >( scope, "BoundaryIntersectionIterator", typeName, IncludeFiles{ "dune/python/range.hh" } );
      if( entry.second )
        registerPyIterator< Iterator >( scope, entry.first );
    }



#if 0
    // registerPyGridViewPartitionIntersectionIterator
    // -----------------------------------------------

    template< class GridView, PartitionIteratorType partition >
    inline static void registerPyGridViewPartitionIntersectionIterator ( pybind11::handle scope = {} )
    {
      typedef PyGridViewPartitionIntersectionIterator< GridView, partition > Iterator;
      auto typeName = GenerateTypeName( "Dune::Python::PyGridViewPartitionIntersectionIterator", MetaType< GridView >(), static_cast< int >( partition ) );
      auto entry = insertClass< Iterator >( scope, "PartitionIntersectionIterator", typeName, IncludeFiles{ "dune/python/range.hh" } );
      if( entry.second )
        registerPyIterator< Iterator >( scope, entry.first );
    }
#endif



    // GridViewPartition
    // -----------------

    template< class GridView, PartitionIteratorType partition >
    struct DUNE_PRIVATE GridViewPartition
    {
      explicit GridViewPartition ( pybind11::object o )
        : gridView( pybind11::cast< const GridView & >( o ) ), obj( std::move( o ) )
      {}

      const GridView &gridView;
      pybind11::object obj;
    };



    // registerPyGridViewPartitionIterator
    // -----------------------------------

    template< class GridView, int codim, PartitionIteratorType partition >
    inline static auto registerPyGridViewPartitionIterator ( pybind11::handle scope, PriorityTag< 1 > )
      -> std::enable_if_t< Capabilities::canIterate< typename GridView::Grid, codim >::value >
    {
      typedef PyGridViewPartitionIterator< GridView, codim, partition > Iterator;
      auto typeName = GenerateTypeName( "Dune::Python::PyGridViewPartitionIterator", MetaType< GridView >(), codim, static_cast< int >( partition ) );
      auto entry = insertClass< Iterator >( scope, "EntityIterator", typeName, IncludeFiles{ "dune/python/grid/range.hh" } );
      if( entry.second )
        registerPyIterator( scope, entry.first );
    }

    template< class GridView, int codim, PartitionIteratorType partition >
    inline static void registerPyGridViewPartitionIterator ( pybind11::handle scope, PriorityTag< 0 > )
    {}

    template< class GridView, int codim, PartitionIteratorType partition >
    inline static void registerPyGridViewPartitionIterator ( pybind11::handle scope = {} )
    {
      return registerPyGridViewPartitionIterator< GridView, codim, partition >( scope, PriorityTag< 42 >() );
    }



    // makePyGridViewPartitionIterator
    // -------------------------------

    template< class GridView, int codim, PartitionIteratorType partition >
    inline static auto makePyGridViewPartitionIterator ( pybind11::object obj )
      -> std::enable_if_t< Capabilities::canIterate< typename GridView::Grid, codim >::value, pybind11::object >
    {
      typedef GridViewPartition< GridView, partition > Partition;
      const GridView &gridView = pybind11::cast< const Partition & >( obj ).gridView;
      pybind11::object iterator = pybind11::cast( new PyGridViewPartitionIterator< GridView, codim, partition >( gridView.template begin< codim, partition >(), gridView.template end< codim, partition >() ) );
      pybind11::detail::keep_alive_impl( iterator, obj );
      return iterator;
    }

    template< class GridView, int codim, PartitionIteratorType partition >
    inline static auto makePyGridViewPartitionIterator ( pybind11::object )
      -> std::enable_if_t< !Capabilities::canIterate< typename GridView::Grid, codim >::value, pybind11::object >
    {
      throw pybind11::value_error( "Iterators for codimension " + std::to_string( codim ) + " are not implemented." );
    }



    // registerGridViewPartition
    // -------------------------

    template< class GridView, PartitionIteratorType partition >
    inline static void registerGridViewPartition ( pybind11::handle scope = {} )
    {
      typedef GridViewPartition< GridView, partition > Partition;
      typedef typename GridView::Grid Grid;
      typedef PyGridViewPartitionIterator< GridView, 0, partition > PyElementIterator;

      using pybind11::operator""_a;

      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [] ( auto codim ) {
          registerPyGridViewPartitionIterator< GridView, codim, partition >();
        } );

      if( !pybind11::already_registered< Partition >() )
      {
        pybind11::class_< Partition > cls( scope, "Partition" );

        if( Capabilities::canIterate< Grid, 0 >::value )
          cls.def_property_readonly( "elements", [] ( pybind11::object self ) { return makePyGridViewPartitionIterator< GridView, 0, partition >( self ); } );
        if( Capabilities::canIterate< Grid, 1 >::value )
          cls.def_property_readonly( "facets", [] ( pybind11::object self ) { return makePyGridViewPartitionIterator< GridView, 1, partition >( self ); } );
        if( Capabilities::canIterate< Grid, GridView::dimension-1 >::value )
          cls.def_property_readonly( "edges", [] ( pybind11::object self ) { return makePyGridViewPartitionIterator< GridView, GridView::dimension-1, partition >( self ); } );
        if( Capabilities::canIterate< Grid, GridView::dimension >::value )
          cls.def_property_readonly( "vertices", [] ( pybind11::object self ) { return makePyGridViewPartitionIterator< GridView, GridView::dimension, partition >( self ); } );

        std::array< pybind11::object (*) ( pybind11::object ), GridView::dimension+1 > makePyIterators;
        Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &makePyIterators ] ( auto codim ) {
            makePyIterators[ codim ] = makePyGridViewPartitionIterator< GridView, codim, partition >;
          } );
        cls.def( "entities", [ makePyIterators ] ( pybind11::object self, int codim ) {
            if( (codim < 0) || (codim > GridView::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [0, " + std::to_string( GridView::dimension ) + "])." );
            return makePyIterators[ codim ]( self );
          }, "codim"_a );

        registerPyBoundaryIntersectionIterator< GridView, PyElementIterator >();
        cls.def_property_readonly( "boundaryIntersections", [] ( const Partition &self ) {
            const GridView &gv = self.gridView;
            return PyBoundaryIntersectionIterator< GridView, PyElementIterator >( gv, PyElementIterator( gv.template begin< 0, partition >(), gv.template end< 0, partition >() ) );
          }, pybind11::keep_alive< 0, 1 >() );

      }
#if 0
      registerPyGridViewPartitionIntersectionIterator< GridView, partition >();
      cls.def_property_readonly( "intersections", [] ( const Partition &self ) {
          return PyGridViewPartitionIntersectionIterator< GridView, partition >( self.gridView );
        }, pybind11::keep_alive< 0, 1 >() );
#endif
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_RANGE_HH
