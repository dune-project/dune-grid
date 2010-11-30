// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TOPOLOGYFACTORY_HH
#define DUNE_TOPOLOGYFACTORY_HH

#include <vector>
#include <map>

#include <dune/common/fvector.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  /**
   * @brief Provide a factory over the generic topologies
   *
   * This class can be used to dynamically create objects
   * statically bound by there generic topology.
   * The method create returns a pointer to an object depending
   * on the topology id and a key; the dimension corresponding
   * to the topology id is static and is provided by the
   * Traits class. A static method (taking the Topology as template
   * argument is also provided).
   * The Traits class must provide the space dimension
   * the types for the key (Key),
   * the objects returned (Object) and the underlying factory
   * (Factory). This class must have a template method
   * createObject taking a key and returning a pointer to
   * the newly create Object - for destruction call the release
   * method.
   **/
  template <class Traits>
  struct TopologyFactory
  {
    // extract types from Traits class
    static const unsigned int dimension = Traits::dimension;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;
    typedef typename Traits::Factory Factory;
    //! dynamically create objects
    static Object *create(unsigned int topologyId, const Key &key) DUNE_DEPRECATED
    {
      Object *object;
      GenericGeometry::IfTopology< Maker, dimension >::apply( topologyId, key, object );
      return object;
    }
    //! dynamically create objects
    static Object *create(const Dune::GeometryType &gt, const Key &key)
    {
      Object *object;
      GenericGeometry::IfTopology< Maker, dimension >::apply( gt.id(), key, object );
      return object;
    }
    //! statically create objects
    template <class Topology>
    static Object *create(const Key &key)
    {
      return Factory::template createObject<Topology> ( key );
    }
    //! release the object returned by the create methods
    static void release( Object *object)
    {
      delete object;
    }
  private:
    // Internal maker class used in ifTopology helper
    template< class Topology >
    struct Maker
    {
      static void apply ( const Key &key, Object *&object )
      {
        object = create<Topology>( key );
      };
    };
  };


  /** @brief A wrapper for a TopologyFactory providing
   *         singleton storage. Same usage as TopologyFactory
   *         but with empty release method an internal storage.
   **/
  template <class Factory>
  struct TopologySingletonFactory
  {
    static const unsigned int dimension = Factory::dimension;
    typedef typename Factory::Key Key;
    typedef const typename Factory::Object Object;
    //! @copydoc TopologyFactory::create(const unsigned int topologyId,const Key &key)
    static Object *create ( const unsigned int topologyId, const Key &key ) DUNE_DEPRECATED
    {
      assert( topologyId < numTopologies );
      return instance().getObject( topologyId, key );
    }
    //! @copydoc TopologyFactory::create(const Dune::GeometryType &gt,const Key &key)
    static Object *create(const Dune::GeometryType &gt, const Key &key)
    {
      assert( gt.id() < numTopologies );
      return instance().getObject( gt.id(), key );
    }
    //! @copydoc TopologyFactory::create(const Key &key)
    template< class Topology >
    static Object *create ( const Key &key )
    {
      dune_static_assert( (Topology::dimension == dimension),
                          "Topology with incompatible dimension used" );
      return instance().template getObject< Topology >( key );
    }
    //! @copydoc TopologyFactory::release
    static void release ( Object *object )
    {}
  private:
    static TopologySingletonFactory &instance ()
    {
      static TopologySingletonFactory instance;
      return instance;
    }

    static const unsigned int numTopologies = (1 << dimension);
    typedef FieldVector< Object *, numTopologies > Array;
    typedef std::map< Key, Array > Storage;

    TopologySingletonFactory ()
    {}
    ~TopologySingletonFactory ()
    {
      const typename Storage::iterator end = storage_.end();
      for( typename Storage::iterator it = storage_.begin(); it != end; ++it )
      {
        for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        {
          Object *&object = it->second[ topologyId ];
          if( object != 0 )
            Factory::release( object );
          object = 0;
        }
      }
    }

    Object *&find( const unsigned int topologyId, const Key &key )
    {
      typename Storage::iterator it = storage_.find( key );
      if( it == storage_.end() )
        it = storage_.insert( std::make_pair( key, Array( 0 ) ) ).first;
      return it->second[ topologyId ];
    }

    Object *getObject ( const unsigned int topologyId, const Key &key )
    {
      Object *&object = find(topologyId,key);
      if( object == 0 )
        object = Factory::create( topologyId, key );
      return object;
    }
    template< class Topology >
    Object *getObject ( const Key &key )
    {
      Object *&object = find(Topology::id,key);
      if( object == 0 )
        object = Factory::template create< Topology >( key );
      return object;
    }
    Storage storage_;
  };

}

#endif // #ifndef DUNE_TOPOLOGYFACTORY_HH
