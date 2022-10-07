// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PERSISTENTCONTAINERINTERFACE_HH
#define DUNE_PERSISTENTCONTAINERINTERFACE_HH

#ifndef HEADERCHECK
#error "This header exists for documentation purposes only and should never be included directly."
#endif

namespace Dune
{

  /** \brief Persistent storage of data on all entities of a grid.
   *
   *  This container allows to store data which is to remain persistent
   *  even during adaptation cycles.
   *  It provides storage for all entities in the hierarchy of a given
   *  codimension (provided dynamically during construction) and behaves much
   *  like an STL container.
   *
   *  The container stores one entry for each entity in the hierarchical grid.
   *  However, it may also store some additional entries, which are not (or no
   *  longer) attached to an entity.
   *
   *  After grid modification the method resize must be called to ensure entries
   *  for each entity in the modified grid.
   *  Accessing newly created entities before calling resize results in
   *  undefined behavior (e.g., a segmentation fault).
   *  To reduce the amount of overallocated entries, the method shrinkToFit
   *  may be called.
   *  It is explicitly possible that the grid adapts any persistent containers
   *  directly during the adaptation process.
   *
   *  The containers are also be persistent over backup / restore of the grid.
   *  After 'shrinkToFit', the entries in the container (and their order) must
   *  match those of a newly created container, even after a backup and restore
   *  of the grid.
   *
   *  There is a default implementation based on std::map but a grid
   *  implementation may provide a specialized implementation.
   *  Grids with a hashable id type can use std::unordered_map to store
   *  the data by simply deriving their PersistentContainer from
   *  Dune::PersistentContainerMap.
   *  For grids providing an id set suitable addressing vector-like storages,
   *  i.e., the id is an integral type and a method size() is provided,
   *  Dune::PersistentContainerVector can be used.
   *
   * \tparam G Grid type
   * \tparam T Container's value type
   */
  template< class G, class T >
  class PersistentContainerInterface
  {
    typedef PersistentContainerInterface< G, T > This;

    struct ImplementationDefined;

  public:
    typedef G Grid;

    typedef T Value;

    typedef ImplementationDefined Size;
    typedef ImplementationDefined ConstIterator;
    typedef ImplementationDefined Iterator;

    // construction

    /** \brief constuctor
     *
     *  \param[in]  grid       reference to the grid to attach data to
     *  \param[in]  codim      codimension to attach data to
     *  \param[in]  value      value for initial entries
     *
     *  \note All implementations must provide at least this constructor.
     */
    PersistentContainerInterface ( Grid &grid, int codim, const Value &value = Value() );

    /** \brief copy constructor */
    PersistentContainerInterface ( const This &other );

    /** \brief assignment operator */
    const This &operator= ( const This &other );

    // element access

    /** \brief access the data associated with an entity
     *
     *  \note The entity must be of the same codimension as the container
     */
    template< class Entity >
    const Value &operator[] ( const Entity &entity ) const;

    /** \brief access the data associated with an entity
     *
     *  \note The entity must be of the same codimension as the container
     */
    template< class Entity >
    Value &operator[] ( const Entity &entity );

    /** \brief access the data associated with a subentity
     *
     *  \note The codimension of the entity must be less or equal to the
     *        codimension of the container.
     */
    template< class Entity >
    const Value &operator() ( const Entity &entity, int subEntity ) const;

    /** \brief access the data associated with a subentity
     *
     *  \note The codimension of the entity must be less or equal to the
     *        codimension of the container.
     */
    template< class Entity >
    Value &operator() ( const Entity &entity, int subEntity );

    // capacity

    /** \brief number of entries in the container
     *
     *  \note The number of entries in the container may be large than the
     *        number of entities in the grid.
     */
    Size size () const;

    /** \brief reserve memory for all entities in the grid
     *
     *  This method will enlarge the container such that there is an
     *  entry for each entity.
     *
     *  After a grid modification, this method must be called before accessing
     *  any data associated to newly created entities.
     *
     *  \note The container might still hold entries for entities that no
     *        longer exist in the grid.
     *        While those entries ca no longer be accessed through an entity,
     *        iterators will still stop at such entries.
     */
    void resize ( const Value &value = Value() );

    /** \brief remove unnecessary entries from container
     *
     *  This method will remove entries from the container that can no longer
     *  be accessed from the grid.
     *
     *  The entries in the container (and their order) must match those of a
     *  newly created container.
     *  This property must be persistent over backup / restore of the grid.
     *
     *  \note This method is merely a hint to the container.
     *        According to the container's internal addressing rules, some
     *        inaccessible entries might remain in the container.
     */
    void shrinkToFit ();

    // modifiers

    /** \brief set all accessible entries to a given value
     *
     *  \note It is undefined, whether the inaccessible values are also
     *        motified.
     */
    void fill ( const Value &value );

    /** \brief exchange the content of the container with another one
     *
     *  \note std::swap is overloaded to refor to this method
     */
    void swap ( This &other );

    // iterators

    /** \brief returns an iterator pointing to the first element of the container
     *
     *  \note Iterators stop at all entries of the container, even if they are
     *        no longer accessible from the grid.
     */
    ConstIterator begin () const;
    /** \brief returns an iterator pointing to the first element of the container
     *
     *  \note Iterators stop at all entries of the container, even if they are
     *        no longer accessible from the grid.
     */
    Iterator begin ();

    /** \brief returns an iterator pointing to the last element of the container
     *
     *  \note Iterators stop at all entries of the container, even if they are
     *        no longer accessible from the grid.
     */
    ConstIterator end () const;
    /** \brief returns an iterator pointing to the last element of the container
     *
     *  \note Iterators stop at all entries of the container, even if they are
     *        no longer accessible from the grid.
     */
    Iterator end ();

    // information

    /** \brief return the codimension, the container attaches data to */
    int codimension () const;
  };

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERINTERFACE_HH
