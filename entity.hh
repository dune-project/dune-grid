// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_ENTITY_HH
#define DUNE_PYTHON_GRID_ENTITY_HH

#include <string>
#include <tuple>
#include <utility>

#include <dune/common/visibility.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/grid/geometry.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/operators.h>

namespace Dune
{

  namespace Python
  {

    // PyHierarchicIterator
    // --------------------

    template<class Entity>
    struct PyHierarchicIterator
    {
      typedef typename Entity::HierarchicIterator HierarchicIterator;

      PyHierarchicIterator(const Entity& e, int maxLevel)
        : it_(e.hbegin(maxLevel)), end_(e.hend(maxLevel))
      {}

      auto next()
      {
        if (it_ == end_)
          throw pybind11::stop_iteration();

        return *it_++;
      }

    private:
      HierarchicIterator it_;
      HierarchicIterator end_;
    };


    // registerPyHierarchicIterator
    // ----------------------------

    template<class Entity>
    auto registerPyHierarchicIterator(pybind11::handle scope)
    {
      typedef PyHierarchicIterator<Entity> Iterator;

      auto cls = insertClass<Iterator>(scope, "PyHierarchicIterator",
          GenerateTypeName("PyHierarchicIterator",MetaType<Entity>()),
          IncludeFiles{"dune/python/entity.hh"}).first;
       cls.def("__iter__", [] (Iterator& it) -> Iterator& { return it; });
       cls.def("__next__", &Iterator::next);
    }



    namespace detail
    {

      template< class Entity, int codim >
      auto registerEntitySubEntity_()
      {
        return [](const Entity& entity, int i) {
          return pybind11::cast(entity.template subEntity<codim>(i));
        };
      }

      template< class Entity, int... codim >
      std::array<std::function<pybind11::object(const Entity& e, int i)>, sizeof...(codim)>
      registerEntitySubEntity(std::integer_sequence<int, codim...>)
      {
        return { {(registerEntitySubEntity_<Entity, codim>())...} };
      }

      // registerExtendedEntityInterface
      // -------------------------------

      template< class Entity, std::enable_if_t< Entity::codimension == 0, int > = 0 >
      inline static void registerExtendedEntityInterface ( pybind11::class_< Entity > &cls )
      {
        cls.def_property_readonly( "father", [] ( const Entity &self ) {
            return (self.hasFather() ? pybind11::cast( self.father() ) : pybind11::none());
          } );

        cls.def_property_readonly( "geometryInFather", &Entity::geometryInFather );

        cls.def( "subEntities", &Entity::subEntities );
        cls.def( "isLeaf", &Entity::isLeaf );
        cls.def( "isRegular", &Entity::isRegular );
        cls.def( "isNew", &Entity::isNew );
        cls.def( "mightVanish", &Entity::mightVanish );
        cls.def( "hasBoundaryIntersections", &Entity::hasBoundaryIntersections );

        registerPyHierarchicIterator< Entity >( cls );
        cls.def( "descendants", [] ( const Entity &e, int maxLevel ) {
            return PyHierarchicIterator< Entity >( e, maxLevel );
          }, pybind11::keep_alive< 0, 1 >() );

        // STATIC variable needed? - visibiility problems?
        static const auto subEntity
          = registerEntitySubEntity<Entity>(std::make_integer_sequence<int, Entity::dimension+1>{});
        cls.def( "subEntity",
                 [](const Entity& e, int i, int codim) { return subEntity.at(codim)(e, i); },
                 pybind11::arg("i"), pybind11::arg("codim") );
      }

      template< class Cls >
      inline static void registerExtendedEntityInterface ( Cls &cls )
      {}



      // registerGridEntity
      // ------------------

      template< class Entity, class... options >
      inline static void registerGridEntity ( pybind11::handle scope,
          pybind11::class_< Entity, options... > cls )
      {
        cls.def_property_readonly( "codimension", [] ( const Entity &e) -> int { return e.codimension; } );
        cls.def_property_readonly( "dimension", [] ( const Entity &e ) -> int { return e.dimension; } );
        cls.def_property_readonly( "mydimension", [] ( const Entity &e ) -> int { return e.mydimension; } );

        cls.def_property_readonly( "geometry", &Entity::geometry );
        cls.def_property_readonly( "level", &Entity::level );
        cls.def_property_readonly( "type", &Entity::type );
        cls.def_property_readonly( "partitionType", &Entity::partitionType );
        cls.def_property_readonly( "domain", [](Entity &self) { return referenceElement<double,self.dimension>(self.type()); },
            pybind11::keep_alive<0,1>() );

        cls.def( pybind11::self == pybind11::self );
        cls.def( pybind11::self != pybind11::self );

        registerExtendedEntityInterface( cls );
      }

    } // namespace detail



    // registerGridEntity
    // ------------------

    template< class Grid, int codim >
    auto registerGridEntity ( pybind11::handle scope )
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;

      auto entry = insertClass<Entity>(scope, "Entity_" + std::to_string(codim),
         GenerateTypeName(MetaType<Grid>(), "template Codim<" + std::to_string(codim) + ">::Entity"));
      if (entry.second)
      {
        detail::registerGridEntity< Entity >( scope, entry.first );
        registerGridGeometry<Entity>( entry.first );
      }
      return entry.first;
    }

    // registerGridEntities
    // --------------------

    template< class Grid, int... codim >
    inline static auto registerGridEntities ( pybind11::handle scope, std::integer_sequence< int, codim... > )
    {
      return std::make_tuple( registerGridEntity< Grid, codim >( scope )... );
    }

    template< class Grid >
    auto registerGridEntities ( pybind11::handle scope )
    {
      return registerGridEntities< Grid >( scope, std::make_integer_sequence< int, Grid::dimension+1 >() );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_ENTITY_HH
