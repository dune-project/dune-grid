// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_MAPPER_HH
#define DUNE_PYTHON_GRID_MAPPER_HH

#include <functional>

#include <dune/common/visibility.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/python/grid/commops.hh>
#include <dune/python/grid/range.hh>
#include <dune/python/utility/numpycommdatahandle.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/numpy.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // mapperSubIndex
      // --------------

      template< class Mapper, class Entity >
      inline static pybind11::object mapperSubIndex ( const Mapper &mapper, const Entity &entity, int i, int c )
      {
        if( (c < Entity::codimension) || (c > Entity::dimension) )
          throw pybind11::value_error( "Invalid codimension: " + std::to_string( c ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
        const int size = entity.subEntities( c );
        if( (i < 0) || (i >= size) )
          throw pybind11::value_error( "Invalid index: " + std::to_string( i ) + " (must be in [0, " + std::to_string( size ) + "))." );
        typename Mapper::Index index;
        return (mapper.contains( entity, i, c, index ) ? pybind11::cast( index ) : pybind11::none());
      }



      // registerMapperSubIndex
      // ----------------------

      template< class Entity, class Mapper, class... options >
      inline static std::enable_if_t< Entity::codimension == 0 >
      registerMapperSubIndex ( pybind11::class_< Mapper, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        pybind11::options opts;
        opts.disable_function_signatures();

        cls.def( "subIndex", [] ( const Mapper &self, const Entity &entity, int i, int c ) {
            return detail::mapperSubIndex( self, entity, i, c );
          }, "entity"_a, "i"_a, "c"_a,
          R"doc(
            Args:
                entity:     an entity of the grid
                i:          number of subentity of `entity`
                c:          codimension of the subentity

            Returns: the starting index of the dofs attached to the `ith` subentity
                      of codimension `c` of the given entity.
          )doc" );

        cls.def( "subIndex", [] ( const Mapper &self, const Entity &entity, std::tuple< int, int > subEntity ) {
            return detail::mapperSubIndex( self, entity, std::get< 0 >( subEntity ), std::get< 1 >( subEntity ) );
          }, "entity"_a, "subEntity"_a,
          R"doc(
            Args:
                entity:     an entity of the grid
                (i,c):      number and codimension of a subentity of `entity`

            Returns: the starting index of the dofs attached to the `ith` subentity
                      of codimension `c` of the given entity.
          )doc" );

        cls.def( "subIndices", [] ( const Mapper &self, const Entity &entity, int codim ) {
            if( (codim < Entity::codimension) || (codim > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            const int size = entity.subEntities( codim );
            pybind11::tuple subIndices( size );
            for( int i = 0; i < size; ++i )
            {
              typename Mapper::Index index;
              subIndices[ i ] = (self.contains( entity, i, codim, index ) ? pybind11::cast( index ) : pybind11::none());
            }
            return subIndices;
          }, "entity"_a, "codim"_a,
          R"doc(
            Args:
                entity:     an entity of the grid
                codim:      codimension of a subentity of `entity`

            Returns: indices of dofs attached to all subentity
                     of codimension `codim` of the given entity.
          )doc" );
      }

      template< class Entity, class Mapper, class... options >
      inline static void registerMapperSubIndex ( pybind11::class_< Mapper, options... > cls, PriorityTag< 0 > )
      {}

      template< class Entity, class Mapper, class... options >
      inline static void registerMapperSubIndex ( pybind11::class_< Mapper, options... > cls )
      {
        registerMapperSubIndex< Entity >( cls, PriorityTag< 42 >() );
      }

      template <class Mapper>
      void mapperCommunicate(const Mapper &mapper, CommOp commOp,
                             std::vector<pybind11::array_t<double>> data,
                             InterfaceType iftype, CommunicationDirection dir )
      {
        switch (commOp)
        {
        case add:
          {
          auto dataHandle = numPyCommDataHandle( mapper, data,
              [] (double local, double remote) {return remote+local;});
          mapper.gridView().communicate( dataHandle, iftype, dir );
          }
        case set:
          {
          auto dataHandle = numPyCommDataHandle( mapper, data,
              [] (double local, double remote) {return remote;});
          mapper.gridView().communicate( dataHandle, iftype, dir );
          }
        }
      }
    } // namespace detail

    // registerMapperCommunicate
    // -------------------------
    template < Dune::PartitionIteratorType fromType,
               Dune::PartitionIteratorType toType,
               Dune::InterfaceType interfaceType,
               class Mapper >
    inline static void registerMapperCommunicate(pybind11::class_<Mapper> cls)
    {
      using pybind11::operator""_a;
      typedef typename Mapper::GridView GridView;
      cls.def( "communicate",
        [] ( const Mapper &self,
          GridViewPartition< GridView, fromType > from,
          GridViewPartition< GridView, toType > to,
          detail::CommOp commOp, pybind11::args args ) {
            std::vector<pybind11::array_t<double>> data(args.size());
            for (int i=0;i<args.size();++i)
              data[i] = pybind11::array_t<double>(args[i]);
            mapperCommunicate(self,commOp, data, interfaceType, ForwardCommunication);
          }, "from"_a, "to"_a, "commOp"_a);
      cls.def( "communicate",
        [] ( const Mapper &self,
          GridViewPartition< GridView, fromType > from,
          GridViewPartition< GridView, toType > to,
          pybind11::function operation, pybind11::args args ) {
            std::vector<pybind11::array_t<double>> data(args.size());
            for (int i=0;i<args.size();++i)
              data[i] = pybind11::array_t<double>(args[i]);
            auto dataHandle = numPyCommDataHandle( self, data, operation.template cast<std::function<double(double,double)>>());
            self.gridView().communicate( dataHandle, interfaceType, ForwardCommunication);
          }, "from"_a, "to"_a, "operation"_a );
      if (fromType != toType)
      {
        cls.def( "communicate",
          [] ( const Mapper &self,
            GridViewPartition< GridView, fromType > from,
            GridViewPartition< GridView, toType > to,
            detail::CommOp commOp, pybind11::args args ) {
              std::vector<pybind11::array_t<double>> data(args.size());
              for (int i=0;i<args.size();++i)
                data[i] = pybind11::array_t<double>(args[i]);
              mapperCommunicate(self,commOp, data, interfaceType, BackwardCommunication);
            }, "from"_a, "to"_a, "commOp"_a );
        cls.def( "communicate",
          [] ( const Mapper &self,
            GridViewPartition< GridView, fromType > from,
            GridViewPartition< GridView, toType > to,
            pybind11::function operation, pybind11::args args ) {
              std::vector<pybind11::array_t<double>> data(args.size());
              for (int i=0;i<args.size();++i)
                data[i] = pybind11::array_t<double>(args[i]);
              auto dataHandle = numPyCommDataHandle( self, data, operation.template cast<std::function<double(double,double)>>());
              self.gridView().communicate( dataHandle, interfaceType, BackwardCommunication);
            }, "from"_a, "to"_a, "operation"_a );
      }
    }

    // registerMapper
    // --------------

    template< class GridView, class Mapper, class... options >
    inline static void registerMapper ( pybind11::class_< Mapper, options... > cls )
    {
      using pybind11::operator""_a;

      pybind11::options opts;
      opts.disable_function_signatures();

      cls.def( "__len__", [] ( const Mapper &self ) { return self.size(); },
        R"doc(
          Returns the maximum index plus 1 returns by the mapper.
          Use to set the number of entries in the dof vector.

          Note: corresponds to the `size()` method on the MCMGMapper.
        )doc" );
      cls.def_property_readonly( "size", [] ( const Mapper &self ) { return self.size(); },
        R"doc(
          Returns the maximum index plus 1 returns by the mapper.
          Use to set the number of entries in the dof vector.
        )doc" );

      cls.def( "types", [] ( const Mapper &self, int codim ) {
          pybind11::list types;
          for( GeometryType type : self.types( codim ) )
            types.append( pybind11::cast( type ) );
          return types;
        }, "codim"_a,
        R"doc(
          Args:
              codim:      codimension

          Returns: list of geometry types of the given codimension with attached data
        )doc" );

      Hybrid::forEach( std::make_integer_sequence< int, GridView::dimension+1 >(), [ &cls ] ( auto &&codim ) {
          typedef typename GridView::template Codim< codim >::Entity Entity;

          using pybind11::operator""_a;

          cls.def( "index", [] ( const Mapper &self, const Entity &e ) -> pybind11::object {
              typename Mapper::Index index;
              if( self.contains( e, index ) )
                return pybind11::cast( index );
              else
                return pybind11::none();
            }, "entity"_a,
            R"doc(
              Args:
                  entity:      entity in the grid

              Returns:  the first index for dofs stored on the given entity
                        if no dofs are stored on this entity the return value is `None`
            )doc" );
          cls.def( "contains", [] ( const Mapper &self, const Entity &entity ) {
              typename Mapper::Index index;
              return self.contains( entity, index );
            }, "entity"_a,
            R"doc(
              Args:
                  entity:      entity in the grid

              Returns:  `True` if the data is stored on this entity
            )doc" );
          detail::registerMapperSubIndex< Entity >( cls );
        } );

      cls.def( "update", [] ( Mapper &self ) { self.update(); },
        R"doc(
          Update the mapper after a grid modification.

          Note: after grid modification the mapper is not valid anymore and needs to
                be updated!
        )doc" );

      cls.def( "__call__", [] ( const Mapper &self, const typename GridView::template Codim< 0 >::Entity &element ) {
          // need a cache gt(cdim=0) -> nof indices then we could store directly in retArray
          std::vector< typename Mapper::Index > indices;
          for( int c = GridView::dimension; c>=0; --c )
            for( auto se : range( element.subEntities( c ) ) )
            {
              const auto &i = self.indices( element, se, c );
              indices.insert( indices.end(), i.begin(), i.end() );
            }
          pybind11::array_t< std::size_t > retArray( indices.size() );
          auto y = retArray.template mutable_unchecked< 1 >();
          std::size_t idx = 0;
          for ( auto i : indices )
            y[idx++] = i;
          return retArray;
        }, "element"_a,
        R"doc(
          Args:
              element:      codim zero entity in the grid

          Returns: list of all degrees of freedom attached to the given and
                   all subentities of the given entity.
        )doc" );


      registerMapperCommunicate<InteriorBorder_Partition,InteriorBorder_Partition,
          InteriorBorder_InteriorBorder_Interface>(cls);
      registerMapperCommunicate<InteriorBorder_Partition,All_Partition,
          InteriorBorder_All_Interface>(cls);
      registerMapperCommunicate<Overlap_Partition,OverlapFront_Partition,
          Overlap_OverlapFront_Interface>(cls);
      registerMapperCommunicate<Overlap_Partition,All_Partition,
          Overlap_All_Interface>(cls);
      registerMapperCommunicate<All_Partition,All_Partition,
          All_All_Interface>(cls);
      cls.def( "communicate",
        [] ( const Mapper &self,
          pybind11::object from, pybind11::object to, pybind11::object, pybind11::args args ) {
            throw pybind11::value_error("Combination of 'to' and 'from' partition not legal");
          }, "from"_a, "to"_a, "commOp"_a,
        R"doc(
          Communicate data stored in `numpy` arrays attached using the given mapper

          Args:
              from:       partition used on sending processor
              to:         partition used on receiving processor
              commOp:     function (local,remote)->result defining how to combine received with locally stored data on receiving process
              *args:      data vectors to communicate

          Returns: None

          Note:
              - At the moment only double numpy vectors are supported
              - Not all combinations of `from` and `to` partitions are valid.
                Possible combinations are
                - from=interiorBorder and to=interiorBorder|all
                - from=overlap and to=overlapFront|all
                - from=overlapFront and to=overlap
                - from=all and to=interiorBorder|overlap|all,
        )doc" );
    }

    // registerMultipleCodimMultipleGeomTypeMapper
    // -------------------------------------------

    template<typename GridView>
    auto registerMultipleCodimMultipleGeomTypeMapper(pybind11::handle scope)
    {
      typedef MultipleCodimMultipleGeomTypeMapper<GridView> MCMGMapper;
      auto cls = insertClass<MCMGMapper>(scope, "MultipleCodimMultipleGeomTypeMapper",
          GenerateTypeName("Dune::MultipleCodimMultipleGeomTypeMapper", MetaType<GridView>()),
          IncludeFiles{"dune/grid/common/mcmgmapper.hh","dune/python/grid/mapper.hh"}).first;
      registerMapper<GridView>(cls);
      cls.def( pybind11::init( [] ( const GridView& gridView, pybind11::function contains ) {
            return MCMGMapper(gridView,
                  [contains](Dune::GeometryType gt, int griddim)
                  { return (unsigned)(contains(gt).cast<int>()); }
                ); } ),
          pybind11::keep_alive<1, 2>() );
      return cls;
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_MAPPER_HH
