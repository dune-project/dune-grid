// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_COMMON_INDEXIDSET_HH
#define DUNE_GRID_COMMON_INDEXIDSET_HH

#include <vector>
#include <dune/common/exceptions.hh>
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
     \tparam TypesImp     iterator range for all geometry types in domain

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
     of a discretization) to entities of the grid. For efficiency reasons the preferred
     data structure for user data is the array. In order to access the data from the
     entity, its index (with respect to an index set - there may be several) is evaluated
     and used as an index to an array (or some other container providing random access).

     Usually an index set is not used directly but a Mapper is used to
     compute the array index from the information supplied by an index set.

     It is important to note that the index assigned to an entity may change during
     grid modification (i.e. refinement or dynamic load balancing). The user is responsible
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
  template< class GridImp, class IndexSetImp, class IndexTypeImp, class TypesImp >
  class IndexSet
  {
    /* We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet. */
    typedef typename std::remove_const< GridImp >::type::Traits Traits;

  public:
    /** \brief Export the type of the entity used as parameter in the index(...) method */
    template <int cc>
    struct Codim
    {
      typedef typename Traits :: template Codim<cc> :: Entity Entity;
    };

    /** \brief The type used for the indices */
    typedef IndexTypeImp IndexType;

    /** \brief iterator range for geometry types in domain */
    typedef TypesImp Types;

    /** \brief dimension of the grid (maximum allowed codimension) */
    static const int dimension = std::remove_const< GridImp >::type::dimension;

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
    template<int cc>
    IndexType index (const typename Traits::template Codim<cc>::Entity& e) const
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
    template<class Entity>
    IndexType index (const Entity& e) const
    {
      constexpr static int cc = Entity::codimension;
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
     *  \param[in]  e      reference to codimension cc entity
     *  \param[in]  i      number subentity of e within the codimension
     *  \param[in]  codim  codimension of the subentity we're interested in
     *
     *  \note The parameter <tt>codim</tt> denotes the codimension with respect
     *        to the grid, i.e., it must satisfy cc <= codim <= dimension.
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
     *  \note The parameter <tt>codim</tt> denotes the codimension with respect
     *        to the grid, i.e., it must satisfy cc <= codim <= dimension.
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

    /**
     * \brief obtain all geometry types of entities in domain
     *
     * This method returns an iterator range (something that behaves like
     * Dune::IteratorRange) visiting all geometry types of codimension codim
     * in the domain of the index map exactly once.
     * The iterator must implement the concept of a forward iterator (in the
     * sense of the STL).
     * The elements in the iterator range are required to be of type
     * Dune::GeometryType.
     *
     * \param[in]  codim  a valid codimension
     *
     * \return iterator range over Const reference to a vector of geometry types.
     */
    Types types ( int codim ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( (asImp().types( codim )) );
      return asImp().types( codim );
    }

    /** @brief Return total number of entities of given geometry type in entity set \f$E\f$.

       \param[in] type A valid geometry type.
       \return    number of entities (type is auto determined by the
                  implementation. std::size_t is the expected return type).
     */
    auto size (GeometryType type) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().size(type)));
      return asImp().size(type);
    }

    /** @brief Return total number of entities of given codim in the entity set \f$E\f$. This
            is simply a sum over all geometry types.

       \param[in] codim A valid codimension
       \return    number of entities (type is auto determined by the
                  implementation. std::size_t is the expected return type).
     */
    auto size (int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().size(codim)));
      return asImp().size(codim);
    }

    /** @brief Return true if the given entity is contained in \f$E\f$.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
     */
    template<class Entity>
    bool contains (const Entity& e) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().contains(e)));
      return asImp().contains(e);
    }

  protected:
    // Must be explicitly defined although this class should get a default constructor.
    IndexSet() = default;

  public:
    //! Forbid the copy constructor
    IndexSet(const IndexSet&) = delete;
    //! Forbid the assignment operator
    IndexSet& operator=(const IndexSet&) = delete;

  private:
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
    typedef typename std::remove_const< GridImp >::type::Traits Traits;

  public:
    /** \brief The type used for the indices */
    typedef typename Base::IndexType IndexType;

    typedef typename Base::Types Types;

    /** \brief dimension of the grid (maximum allowed codimension) */
    static const int dimension = Base::dimension;

    using Base::index;
    using Base::subIndex;

    //===========================================================
    /** @name Access to entity set
     */
    //@{
    //===========================================================

    Types types ( int codim ) const { return asImp().geomTypes( codim ); }

    /** @brief Return total number of entities of given codim in the entity set \f$E\f$. This
            is simply a sum over all geometry types.

       \param[in] codim A valid codimension
       \return    number of entities (type is auto determined by the
                  implementation. std::size_t is the expected return type).
     */
    auto size ( const int codim ) const
    {
      using SizeType = std::decay_t<decltype( Base::size( Dune::GeometryType() ) )>;

      const std::vector< GeometryType > &geomTs = asImp().geomTypes( codim );
      typedef typename std::vector< GeometryType >::const_iterator Iterator;

      const Iterator end = geomTs.end();

      SizeType s ( 0 );
      for( Iterator it = geomTs.begin() ; it != end; ++it )
        s += Base::size( *it );

      return s;
    }
    //@{

  private:
    IndexSetImp &asImp () { return static_cast< IndexSetImp & >( *this );}
    const IndexSetImp &asImp () const { return static_cast< const IndexSetImp & >( *this ); }
  };


  /** @brief Id Set %Interface.

     This class template is used as a base class for all id set implementations.
     It uses the Barton-Nackman trick to ensure conformity to the interface.

     \tparam GridImp Type that is a model of Dune::Grid.
     \tparam IdSetImp Type that is a model of Dune::IdSet.
     \tparam IdTypeImp Type used for ids

     <H3>Overview</H3>

     An id set provides a map \f[ m : E \to \mathbf{I}\f] where
     \f$E\f$ is a subset of the entities of a grid and \f$\mathbf{I}\f$ is a finite
     set of indices.  These indices  are persistent under grid
     modifications, i.e. if there exists an entity \f$e\f$ with index \f$i\f$
     before grid modification and an entity \f$e^\prime\f$ with index \f$i\f$ after grid
     modification it is guaranteed that \f$e=e^\prime\f$.

     In the terminology of the Dune grid interface, these indices
     are called 'ids'.  Ids are typically numbers, but may be other things, too.
     If they are numbers, the ids used in a grid are not necessarily positive,
     and they are usually not consecutive.
     However, the ids must be usable as keys for STL associative containers,
     like a <tt>std::map</tt> or a <tt>std::unordered_map</tt>.  As ids are persistent,
     they are used to keep data on a grid that undergoes modifications.

     <H3>Injectivity properties</H3>

     The id map \f$m\f$ is injective on the entire set of entities of a hierarchical
     grid.  More formally, we have:

     - For any \f$e,e^\prime\in E\f$ we have \f$e\neq e^\prime \Rightarrow m(e)\neq m(e^\prime)\f$.

     An exception to this rule holds for entities that are copies of entities on a lower refinement
     level.
     An element is a copy of its father element if it is the only son. This
     concept can be transferred to all higher codimensions because in a nested grid
     structure the entities of any codimension form a set of trees. However, the roots
     of these trees are not necessarily on level 0.
     Thus, we define that an entity is a copy of another entity if it is the only descendant
     of this entity in the refinement tree. This is illustrated in the following figure where,
     for example, vertex w is a copy of vertex v.

     \image html  idlocalref.png "Sharing of ids."
     \image latex idlocalref.eps "Sharing of ids." width=\textwidth

     The copy relation is an equivalence relation.
     We define that <EM> all copies of an entity share the same id.</EM>
     In the example of the figure the vertices v and w would have the same id.
     This definition is useful to transfer data on the leaf grid during grid modification.

     <b>Note:</b> In "Bastian et al.: A Generic Grid Interface for Adaptive and Parallel Scientific Computing.
     Part I: Abstract Framework, Computing 82 (2-3), 2008, pp. 103-119" it is claimed that the
     id map should be injective only for entities of the same codimension.  This is in
     contrast to the actual %Dune code: in the code, grid implementations are required to provide
     ids that are injective even on sets of entities with different codimensions.

     <H3>Properties of the <tt>IdTypeImp</tt> type</H3>

     The type <tt>IdTypeImp</tt> used for ids will be specified by the grid implementation.
     Conceptually, it can be pretty much anything, but we require that it is
     usable as <tt>key</tt> for the <tt>std::map</tt> and <tt>std::unordered_map</tt>
     containers of the standard library.  In particular, for <tt>std::map</tt> this means
     that <tt>IdTypeImp</tt> implements
     \code
     bool operator<  ( const IdTypeImp &, const IdTypeImp & );
     \endcode
     and that this operator implements a strict weak ordering.  For use with <tt>std::unordered_map</tt>,
     we additionally require
     \code
     bool operator== ( const IdTypeImp &, const IdTypeImp & );
     \endcode
     and that <tt>std::hash<IdType></tt> fulfills the requirements of an stl hash object.

     For convenience we further require that
     \code
     bool operator!= ( const IdTypeImp &, const IdTypeImp & );
     \endcode
     is implemented and returns the negation of <tt>operator==</tt>.

     The <tt>IdTypeImp</tt> must be
     <ul>
      <li> default-constructible </li>
      <li> copy-constructible </li>
      <li> copy-assignable </li>
     </ul>

     Finally, for debugging purposes it must be possible to write objects of type
     <tt>IdTypeImp</tt> into standard C++ streams.  Therefore
     \code
     template< class C, class T >
     std::basic_ostream< C, T > &operator<< ( std::basic_ostream< C, T > &, const IdTypeImp & );
     \endcode
     has to exist and do something reasonable.

     <H3>Global and local id sets</H3>

     Grids are required to provide two types of id sets.  These differ only of the
     grid is distributed over several processes:

     <ul>
      <li>
       A global id set provides ids that are unique over all processes over which the
       grid is distributed.
      </li>
      <li>
     A local id set provides ids that are unique within one process, but two entities
     in different processes may have the same id. Obviously, a global id set is also
     a local id set, but a dedicated local id set implementation may be more efficient.
      </li>
     </ul>

     @ingroup IndexIdSets
   */
  template<class GridImp, class IdSetImp, class IdTypeImp>
  class IdSet
  {
    /* We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet. */
    using Traits = typename std::remove_const< GridImp >::type::Traits;
  public:
    //! Type used to represent an id.
    typedef IdTypeImp IdType;

    /** \brief Export the type of the entity used as parameter in the id(...) method */
    template <int cc>
    struct Codim {
      using Entity = typename Traits::template Codim<cc>::Entity;
    };

    /** \brief dimension of the grid (maximum allowed codimension) */
    static constexpr auto dimension = std::remove_const< GridImp >::type::dimension;

    //! Get id of an entity. This method is simpler to use than the one below.
    template<class Entity>
    IdType id (const Entity& e) const
    {
      constexpr static int cc = Entity::codimension;
      return asImp().template id<cc>(e);
    }

    //! Get id of an entity of codim cc. Unhandy because template parameter must be supplied explicitly.
    template<int cc>
    IdType id (const typename Codim<cc>::Entity& e) const
    {
      return asImp().template id<cc>(e);
    }

    /** \brief Get id of subentity i of co-dimension codim of a co-dimension 0 entity.
     */
    IdType subId (const typename Codim<0>::Entity& e, int i, unsigned int codim) const
    {
      return asImp().subId(e,i,codim);
    }

  protected:
    // Default constructor (is not provided automatically because copy constructor is private)
    IdSet() = default;

  public:
    //! Forbid the copy constructor
    IdSet(const IdSet&) = delete;
    //! Forbid the assignment operator
    IdSet& operator=(const IdSet&) = delete;

  private:
    //!  Barton-Nackman trick
    IdSetImp& asImp () {return static_cast<IdSetImp &> (*this);}
    //!  Barton-Nackman trick
    const IdSetImp& asImp () const {return static_cast<const IdSetImp &>(*this);}
  };

}

#endif   // DUNE_GRID_COMMON_INDEXIDSET_HH
