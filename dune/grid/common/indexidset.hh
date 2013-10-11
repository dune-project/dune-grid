// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#ifndef DUNE_INDEXIDSET_HH
#define DUNE_INDEXIDSET_HH

#include <iostream>
#include <vector>
#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>
#include <dune/grid/common/grid.hh>


/** @file
        @author Peter Bastian
        @brief Provides base classes for index and id sets
 */

namespace Dune
{

#include <dune/common/bartonnackmanifcheck.hh>

  /** @brief Index Set %Interface base class.

     This class template is used as a base class for all index set implementations.
     It uses the Barton-Nackman trick to ensure conformity to the interface.

     \tparam GridImp Type that is a model of Dune::Grid.
     \tparam IndexSetImp Type that is a model of Dune::IndexSet.
     \tparam IndexTypeImp The type used by IndexSetImp to store the indices

     <H3>Overview</H3>

     An index set provides a map \f[ m : E \to \mathbf{N}\f] where
     \f$E\f$ is a subset of the entities of a grid and \f$\mathbf{N}\f$ is the set of
     natural numbers (including 0).

     We define the subsets
     \f[ E_g^c = \{e\in E \ | \ \textrm{$e$ has codimension $c$ and geometry type $g$} \}.\f]

     The index map \f$m\f$ has the following properties:

     - It is unique within the subsets \f$E_g^c\f$, i.e. for any \f$e,e^\prime\in E_g^c\f$
     we have \f$e\neq e^\prime \rightarrow m(e)\neq m(e^\prime)\f$.
     - It is consecutive and zero-starting within the subsets \f$E_g^c\f$, i.e. we have
     \f$0\leq m(e) < |E_g^c|\f$ for any \f$e\in E_g^c\f$.

     Index sets are used to assign user defined data (e.g. degrees of freedom
     of a discretization) to entities of the grid. For efficiency reasons the prefered
     data structure for user data is the array. In order to access the data from the
     entity, its index (with respect to an index set - there may be several) is evaluated
     and used as an index to an array (or some other container providing random access).

     Usually an index set is not used directly but a Mapper is used to
     compute the array index from the information supplied by an index set.

     It is important to note that the index assigned to an entity may change during
     grid modification (i.e. refinement or dynamic load balancing). The user is reponsible
     for reorganizing the information stored in the external arrays appropriately. In
     order to do this the IdSet concept is supplied.

     <H3>Level index</H3>

     Index set where \f$E\f$ corresponds to all entities of a given grid level.
     All grid implementations provide level indices.

     <H3>Leaf Index</H3>

     Index set where \f$E\f$ corresponds to all entities of the leaf grid.
     All grid implementations provide a leaf index.

     @ingroup IndexIdSets
   */
  template<class GridImp, class IndexSetImp, class IndexTypeImp>
  class IndexSet
  {
    /* We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet. */
    typedef typename remove_const< GridImp >::type::Traits Traits;

  public:
    /** \brief The type used for the indices */
    typedef IndexTypeImp IndexType;

    /** \brief dimension of the grid (maximum allowed codimension) */
    static const int dimension = remove_const< GridImp >::type::dimension;

    //===========================================================
    /** @name Index access from entity
     */
    //@{
    //===========================================================

    /** @brief Map entity to index. The result of calling this method with an entity that is not
            in the index set is undefined.

            \param e Reference to codim cc entity, where cc is the template parameter of the function.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cc>
    IndexType index (const typename remove_const<GridImp>::type::
                     Traits::template Codim<cc>::Entity& e) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template index<cc>(e)));
      return asImp().template index<cc>(e);
    }

    /** @brief Map entity to index. Easier to use than the above because codimension template
            parameter need not be supplied explicitly.
            The result of calling this method with an entity that is not
            in the index set is undefined.

            \param e Reference to codim cc entity. Since
           entity knows its codimension, automatic extraction is possible.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<class EntityType>
    IndexType index (const EntityType& e) const
    {
      enum { cc = EntityType::codimension };
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template index<cc>(e)));
      return asImp().template index<cc>(e);
    }

    /** \brief Map a subentity to an index.
     *
     *  The result of calling this method with an entity that is not in the
     *  index set is undefined.
     *
     *  \tparam  cc  codimension of the entity
     *
     *  \param[in]  e      reference to codimsion cc entity
     *  \param[in]  i      number subentity of e within the codimension
     *  \param[in]  codim  codimension of the subentity we're interested in
     *                     (must satisfy cc <= codim <= dimension)
     *
     *  \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template< int cc >
    IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &e,
                         int i, unsigned int codim ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template subIndex< cc >(e,i,codim)));
      return asImp().template subIndex< cc >(e,i,codim);
    }

    /** \brief Map a subentity to an index.
     *
     *  The result of calling this method with an entity that is not in the
     *  index set is undefined.
     *
     *  \note This method exists for convenience only.
     *        It extracts the codimension from the type of the entity, which can
     *        be guessed by the compiler.
     *
     *  \tparam  Entity  type of entity (must be GridImp::Codim< cc >::Entity
     *                   for some cc)
     *
     *  \param[in]  e      reference to entity
     *  \param[in]  i      number subentity of e within the codimension
     *  \param[in]  codim  codimension of the subentity we're interested in
     *
     *  \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template< class Entity >
    IndexType subIndex ( const Entity &e, int i, unsigned int codim ) const
    {
      static const int cc = Entity::codimension;
      return asImp().template subIndex< cc >( e, i, codim );
    }
    //@}


    //===========================================================
    /** @name Access to entity set
     */
    //@{
    //===========================================================

    /** @brief Return vector with all geometry types of entities in domain of index map.
            Return a vector with all geometry types of a given codimension
            contained in the Entity set \f$E\f$.

       \param[in] codim A valid codimension.
       \return Const reference to a vector of geometry types.
     */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().geomTypes(codim)));
      return asImp().geomTypes(codim);
    }

    /** @brief Return total number of entities of given geometry type in entity set \f$E\f$.

       \param[in] type A valid geometry type.
       \return         number of entities.
     */
    IndexType size (GeometryType type) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().size(type)));
      return asImp().size(type);
    }

    /** @brief Return total number of entities of given codim in the entity set \f$E\f$. This
            is simply a sum over all geometry types.

       \param[in] codim A valid codimension
            \return    number of entities.
     */
    IndexType size (int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().size(codim)));
      return asImp().size(codim);
    }

    /** @brief Return true if the given entity is contained in \f$E\f$.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
     */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().contains(e)));
      return asImp().contains(e);
    }

    // Must be explicitely defined although this class should get a default constructor.
    IndexSet() {}

  private:
    //! Forbid the copy constructor
    IndexSet(const IndexSet&);
    //! Forbid the assignment operator
    IndexSet& operator=(const IndexSet&);

    //!  Barton-Nackman trick
    IndexSetImp& asImp () {return static_cast<IndexSetImp &> (*this);}
    //!  Barton-Nackman trick
    const IndexSetImp& asImp () const {return static_cast<const IndexSetImp &>(*this);}
  };

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION



  /**\brief Provide default implementation of method if IndexSet
         @ingroup GridDevel
   */
  template<class GridImp, class IndexSetImp>
  class IndexSetDefaultImplementation
    : public IndexSet< GridImp, IndexSetImp >
  {
    typedef IndexSet< GridImp, IndexSetImp > Base;
    typedef typename remove_const< GridImp >::type::Traits Traits;

  public:
    /** \brief The type used for the indices */
    typedef typename Base::IndexType IndexType;

    /** \brief dimension of the grid (maximum allowed codimension) */
    static const int dimension = Base::dimension;

    using Base::index;
    using Base::subIndex;

    //===========================================================
    /** @name Index access from entity
     */
    //@{
    //===========================================================

    /** \copydoc Dune::IndexSet::subIndex(const typename Traits::template Codim< cc >::Entity &e,int i,unsigned int codim) const
     *
     *  The default implementation is as follows:
     *  \code
     *  index( *(e.subEntity( i, codim )) );
     *  \endcode
     *  It does only work for cc=0 since the subEntity method is not present otherwise.
     */
    template< int cc >
    IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &e, int i, unsigned int codim ) const
    {
      // this does not work, since subEntity is a template method requiring codim to be
      // a template parameter
      // return index( *(e.subEntity( i, codim )) );
      DUNE_THROW(NotImplemented,"subIndex for entities is not is not implemented");
      return -1;
    }
    //@}

    //===========================================================
    /** @name Access to entity set
     */
    //@{
    //===========================================================

    /** @brief Return total number of entities of given codim in the entity set \f$E\f$. This
            is simply a sum over all geometry types.

       \param[in] codim A valid codimension
                \return    number of entities.
     */
    IndexType size ( const int codim ) const
    {
      IndexType s( 0 );
      const std::vector< GeometryType > &geomTs = Base::geomTypes( codim );
      typedef typename std::vector< GeometryType >::const_iterator Iterator;
      const Iterator end = geomTs.end();
      for( Iterator it = geomTs.begin(); it != end; ++it )
        s += Base::size( *it );
      return s;
    }
    //@{
  };


  /** @brief Id Set %Interface.

     This class template is used as a base class for all id set implementations.
     It uses the Barton-Nackman trick to ensure conformity to the interface.

     \tparam GridImp Type that is a model of Dune::Grid.
     \tparam IdSetImp Type that is a model of Dune::IdSet.
     \tparam IdTypeImp Type used for ids

     <H3>Overview</H3>

     An id set provides a map \f[ m : E \to \mathbf{I}\f] where
     \f$E\f$ is a subset of the entities of a grid and \f$\mathbf{I}\f$ is a discrete
     set of ids. These ids need not be consecutive nor positive.
     However, the ids must be usable as keys for STL associative containers
     (e.g., <tt>std::map</tt>). For debugging purposes, it must also be possible
     to write them into standard C++ streams.
     More precisely, for such a type <tt>Id</tt>, at least the following operators
     have to be provided:
     \code
     bool operator== ( const Id &, const Id & );
     bool operator!= ( const Id &, const Id & );
     bool operator<  ( const Id &, const Id & );

     template< class C, class T >
     std::basic_ostream< C, T > &operator<< ( std::basic_ostream< C, T > &, const Id & );
     \endcode

     The index map \f$m\f$ has the following properties:

     - It is injective, i.e. for any \f$e,e^\prime\in E\f$
     we have \f$e\neq e^\prime \Rightarrow m(e)\neq m(e^\prime)\f$.
     - It is persistent with respect to grid modification, i.e. if there exists an entity \f$e\f$ with
     id \f$i\f$ before grid modification and an entity \f$e^\prime\f$ with id \f$i\f$ after mesh
     modification it is guaranteed that \f$e=e^\prime\f$.

     The set of ids \f$ \mathbf{I} = \{i\ |\ \exists e\in E : m(e)=i\}\f$ used by the
     id set is not necessarily consecutive. In practice the numbers can be quite large, especially
     in a parallel implementation. Therefore the type used to represent the id can be chosen
     by the application.

     <H3>Ids and leaf entities</H3>

     An element is a copy of its father element if it is the only son. This
     concept can be transfered to all higher codimensions because in a nested grid
     structure the entities of any codimension form a set of trees. However, the roots
     of these trees are not necessarily on level 0.
     Thus, we define that an entity is a copy of another entity if it is the only descendant
     of this entity in the refinement tree. This is illustrated in the following figure where,
     for example, vertex w is a copy of vertex v.

     \image html  idlocalref.png "Sharing of ids."
     \image latex idlocalref.eps "Sharing of ids." width=\textwidth

     The copy relation can be trivially extended to be an equivalence relation.
     With respect to ids we define that <EM> all copies of an entity share the same id.</EM>
     In the example of the figure the vertices v and w would have the same id.

     This definition is useful to transfer data related to the leaf grid during grid modification.

     <H3>Global id set</H3>

     A global id set provides ids that are unique over all processes over which the
     grid is distributed.
     All grid implementations provide a global id set.

     <H3>Local id set</H3>

     A local id set provides ids that are unique within one process but two entities
     in different processes may have the same id. Obviously, a global id set is also
     a local id set. A grid implementation may provide an extra local id set for efficiency reasons.
     In sequential grids local and global id set are identical.
     All grid implementations provide a local id set.

     @ingroup IndexIdSets
   */
  template<class GridImp, class IdSetImp, class IdTypeImp>
  class IdSet
  {
  public:
    //! Type used to represent an id.
    typedef IdTypeImp IdType;

    //! Get id of an entity. This method is simpler to use than the one below.
    template<class EntityType>
    IdType id (const EntityType& e) const
    {
      enum { cc = EntityType::codimension };
      return asImp().template id<cc>(e);
    }

    //! Get id of an entity of codim cc. Unhandy because template parameter must be supplied explicitely.
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cc>
    IdType id (const typename remove_const<GridImp>::type::
               Traits::template Codim<cc>::Entity& e) const
    {
      return asImp().template id<cc>(e);
    }

    /** \brief Get id of subentity i of co-dimension codim of a co-dimension 0 entity.
     */
    IdType subId (const typename remove_const<GridImp>::type::
                  Traits::template Codim<0>::Entity& e, int i, unsigned int codim) const
    {
      return asImp().subId(e,i,codim);
    }

    // Default constructor (is not provided automatically because copy constructor is private)
    IdSet() {}

  private:
    //! Forbid the copy constructor
    IdSet(const IdSet&);
    //! Forbid the assignment operator
    IdSet& operator=(const IdSet&);

    //!  Barton-Nackman trick
    IdSetImp& asImp () {return static_cast<IdSetImp &> (*this);}
    //!  Barton-Nackman trick
    const IdSetImp& asImp () const {return static_cast<const IdSetImp &>(*this);}
  };

}

#endif
