// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_ENTITYRANGEGENERATORS_HH
#define DUNE_GRID_COMMON_ENTITYRANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

#ifdef DOXYGEN

  /**
   * \addtogroup GIIteration
   * \brief Iterator ranges for entities and intersections to support iteration with range-based for loops.
   *
   * Iterating over sets of entities or intersections is one of the most common operations when
   * writing DUNE code. The grid interface implements this in the standard C++ ways by providing
   * iterators and matching begin() and end() methods, but their usage can be rather unwieldy.
   *
   * This page describes a much simpler alternative based on a C++11 feature called range-based for loop.
   *
   * <h2>Range-based for loop</h2>
   *
   * \note Range-based for loops are availabe in GCC 4.6+, Clang 3.2+ and Intel ICC 13+.
   *
   * A range-based for loop is a short-hand way of iterating over any object that provides the standard
   * begin() and end() methods. It looks like this:
   *
   * \code
   * for (auto&& i : vec)
   *   i *= 2;
   * \endcode
   *
   * This code will multiply all entries of `vec` by 2. The loop head always looks like `for (<type> <variable-name> : <container>)`.
   * You can also specify the exact type of the variable, but it is normally much easier to let the compiler
   * do that for you using `auto`.
   *
   * \note Due to a number of upcoming changes to our grid interface, you should **always**
   *       use `auto&&` to capture the Entity / Intersection type. Otherwise, your code will
   *       break after the release of DUNE 2.4!
   *
   * For those interested in the technical details, the compiler has translated the loop into
   * something resembling this (not quite, but we'll forget about the minor details here):
   *
   * \code
   * for (auto it = vec.begin(),
   *        end = vec.end();
   *      it != end;
   *      ++it)
   * {
   *   auto&& i = *it;
   *   i *= 2;
   * }
   * \endcode
   *
   * For further details, see e.g. http://en.cppreference.com/w/cpp/language/range-for.
   *
   * <h2>Entities</h2>
   *
   * You cannot simply iterate over all entities of a GridView by passing it to a range-based for loop,
   * simply because you cannot (and probably do not want to) iterator over *all* of its entities. Instead,
   * algorithms typically have to iterate over all entities with a particular codimension, e.g. over all
   * grid cells. The functions listed at the top of this page allow you to do just this. Assuming you have
   * a GridView `gv`, you can iterate over its cells and vertices like this:
   *
   * \code
   * // iterate over cells
   * for (auto&& cell : elements(gv))
   * {
   *   std::cout << "Cell " << gv.indexSet().index(cell) << " is centered at "
   *             << cell.geometry().center() << std::endl;
   * }
   * // iterate over vertices
   * for (auto&& vertex : vertices(gv))
   * {
   *   std::cout << "Vertex " << gv.indexSet().index(vertex) << " at "
   *             << vertex.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \note As explained above, **always** use `auto&&` for the type of the Entity!
   *
   * There are also functions for iterating over facets() and edges() as well as generic entities() functions,
   * which allow you to specify a numeric codimension or dimension.
   *
   * If you are using Dune for parallel computations, you probably know that the GridView offers iterators with
   * different \link PartitionIteratorType PartitionIteratorTypes\endlink. Those are also supported with range-based
   * for loops: By default, the functions above will iterate over all entities (Dune::All_Partition), but they
   * accept an additional parameter for selecting a different set of partitions. This parameter is of type
   * Dune::PartitionSet. You can find pre-instantiated objects for all partitions and allowed combinations in the
   * namespace Dune::Partitions. Using those objects, you can iterate over all interior and border vertices like
   * this:
   *
   * \code
   * // use prebuild PartitionSet
   * for (auto&& vertex : vertices(gv,Dune::Partitions::interiorBorder))
   * {
   *    ...
   * }
   * // construct PartitionSet by combining partitions
   * for (auto&& vertex : vertices(gv,Dune::Partitions::interior + Dune::Partitions::border))
   * {
   *    ...
   * }
   * \endcode
   *
   * <h2>Intersections</h2>
   *
   * If you want to iterate over the \link Intersection intersections \endlink of an Entity, you can
   * use the function intersections() to obtain a range that is suitable for a range-based for loop.
   * Intersections are always defined with respect to a GridView; for further information, see the discussion
   * on locally refined grids.
   *
   * As an example, the following code counts the number of boundary intersections for a given GridView `gv`:
   *
   * \code
   * std::size_t count = 0;
   * for (auto&& e : elements(gv))
   * {
   *   do_stuff();
   *   for (auto&& i : intersections(gv,e))
   *   {
   *     if (e.boundary())
   *       ++count;
   *   }
   * }
   * \endcode
   *
   * \{
   */


  // *****************************************************************************************
  // Doxygen documentation
  // *****************************************************************************************
  //
  // In the following, the range generating functions are documented for Doxygen; the actual
  // implementations are further down in this file and hidden from Doxygen.
  // The main reason for this split are the return types of those functions, which either contain
  // long type listings to obtain the iterator type or (in the case of the forwarded functions
  // use the new-style function syntax and calculate the return type using decltype. In both cases,
  // Doxygen generates function signatures that are very confusing to the average user.
  //
  // *****************************************************************************************



  //! \name Common entity ranges
  //! \brief Entity ranges for common entity types. If in doubt, use one of these.
  //! \{


  //! Iterates over all elements / cells (entities with codimension 0) of a GridView.
  /**
   * This functions returns an object representing the range of *elements* or *cells*
   * contained in the GridView gv. The main purpose of this function is to enable iteration
   * over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all cells in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : elements(gv))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This is the default version of the elements() function. It will always iterate over all
   *         elements in the GridView, regardless of their Dune::PartitionType. If you are interested
   *         in cells with specific PartitionType(s), use elements(const GV&,PartitionSet<partitions>)
   *         instead.
   *
   * \sa elements(const GV&,PartitionSet<partitions>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the elements.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV>
  IteratorRange<...> elements(const GV& gv);


  //! Iterates over all facets (entities with codimension 1) of a GridView.
  /**
   * This functions returns an object representing the range of *facets*
   * contained in the GridView gv. The main purpose of this function is to enable iteration
   * over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all facets in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : facets(gv))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This is the default version of the facets() function. It will always iterate over all
   *         elements in the GridView, regardless of their Dune::PartitionType. If you are interested
   *         in cells with specific PartitionType(s), use facets(const GV&,PartitionSet<partitions>)
   *         instead.
   *
   * \sa facets(const GV&,PartitionSet<partitions>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the facets.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV>
  IteratorRange<...> facets(const GV& gv);


  //! Iterates over all edges (entities with dimension 1) of a GridView.
  /**
   * This functions returns an object representing the range of *edges*
   * contained in the GridView gv. The main purpose of this function is to enable iteration
   * over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all edges in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : edges(gv))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This is the default version of the edges() function. It will always iterate over all
   *         elements in the GridView, regardless of their Dune::PartitionType. If you are interested
   *         in cells with specific PartitionType(s), use edges(const GV&,PartitionSet<partitions>)
   *         instead.
   *
   * \sa edges(const GV&,PartitionSet<partitions>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the edges.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV>
  IteratorRange<...> edges(const GV& gv);


  //! Iterates over all vertices (entities with dimension 0) of a GridView.
  /**
   * This functions returns an object representing the range of *vertices*
   * contained in the GridView gv. The main purpose of this function is to enable iteration
   * over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all vertices in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : vertices(gv))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This is the default version of the vertices() function. It will always iterate over all
   *         elements in the GridView, regardless of their Dune::PartitionType. If you are interested
   *         in cells with specific PartitionType(s), use vertices(const GV&,PartitionSet<partitions>)
   *         instead.
   *
   * \sa vertices(const GV&,PartitionSet<partitions>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the vertices.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV>
  IteratorRange<...> vertices(const GV& gv);


  //! \}



  //! \name Intersection Range
  //! \brief Iterator range for \link Intersection intersections \endlink.
  //! \{

  //! Iterates over all \link Intersection Intersections \endlink of an Entity with respect to the given GridView.
  /**
   * This functions returns an object representing the range of `Intersection`s of the Entity
   * e with respect to the GridView gv. The main purpose of this function is to enable
   * iteration over those intersections by means of a range-based for loop:
   *
   * \code
   * // iterate over all intersections of an entity with respect to the LeafGridView
   * auto gv = grid.leafGridView();
   * auto&& entity = ...; // get an entity from somewhere
   * for (auto&& i : intersections(gv,entity))
   * {
   *   std::cout << i.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \relates    GridView
   * \relates    Entity
   * \param gv   the GridView to use for determining interior intesections.
   * \param e    the Entity whose intersections should be iterated over.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, typename Entity>
  IteratorRange<...> intersections(const GV& gv, const Entity& e);


  //! \}


  //! \name Entity ranges with (co)dimension template argument
  //! \brief Entity ranges which allow specifying the codimension / dimension as a numeric template parameter.
  //! \{


  //! Iterates over all entities of a GridView with the given codimension.
  /**
   * This functions returns an object representing the range of entities of codimension
   * cd contained in the GridView gv. The main purpose of this function is to enable iteration
   * over those entities by means of a range-based for loop:
   *
   * \code
   * // iterate over all cells in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : entities(gv,Dune::Codim<0>()))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This function allows you to write loops that are parameterized on the codimension of
   *         the entities. While this allows for extra flexibility (for e.g. some codimension-agnostic
   *         algorithms), it reduces the code readibility. If you don't need this flexibility, consider
   *         using elements(), facets(), edges() or vertices() instead.
   *
   * \remark This is the default version of the entities() function. It will always iterate over all
   *         entities in the GridView, regardless of their Dune::PartitionType. If you are interested
   *         in entities with specific PartitionType(s), use
   *         entities(const GV&,Codim<codim>,PartitionSet<partitions>) instead.
   *
   * \remark If you have to iterate over entities with a specific *dimension*, consider using
   *         entities(const GV&,Dim<dim>) instead to improve the readability of your code.
   *
   * \sa entities(const GV&,Codim<codim>,PartitionSet<partitions>)
   * \sa entities(const GV&,Dim<dim>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the entities.
   * \param cd   a Codim object that is used to specify the codimension of the entities by means
   *             of its template parameter.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, int codim>
  IteratorRange<...> entities(const GV& gv, Codim<codim> cd);


  //! Iterates over all entities of a GridView with the given dimension.
  /**
   * This functions returns an object representing the range of entities of dimension
   * d contained in the GridView gv. The main purpose of this function is to enable iteration
   * over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all edges in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : entities(gv,Dune::Dim<1>()))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode

   * \remark This function allows you to write loops that are parameterized on the codimension of
   *         the entities. While this allows for extra flexibility (for e.g. some codimension-agnostic
   *         algorithms), it reduces the code readibility. If you don't need this flexibility, consider
   *         using elements(), facets(), edges() or vertices() instead.
   *
   * \remark This is the default version of the entities() function. It will always iterate over all
   *         entities in the GridView, regardless of their Dune::PartitionType. If you are interested
   *         in entities with specific PartitionType(s), use
   *         entities(const GV&,Dim<dim>,PartitionSet<partitions>) instead.
   *
   * \remark If you have to iterate over entities with a specific *codimension*, consider using
   *         entities(const GV&,Codim<codim>) instead to improve the readability of your code.
   *
   * \sa entities(const GV&,Dim<dim>,PartitionSet<partitions>)
   * \sa entities(const GV&,Codim<codim>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the entities.
   * \param d    a Dim object that is used to specify the dimension of the entities by means
   *             of its template parameter.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, int dim>
  IteratorRange<...> entities(const GV& gv, Dim<dim> d);

  //! \}



  //! \name Common entity ranges for non-standard parallel partitions
  //! \brief The following Entity ranges make it possible to specify a PartitionSet which is sometimes needed in parallel code.
  //! \{


  //! Iterates over all elements / cells (entities with codimension 0) of a GridView that belong to the given PartitionSet.
  /**
   * This functions returns an object representing the range of *elements* or *cells* contained
   * in the GridView gv which belong to the PartitionSet ps. The main purpose of this function
   * is to enable iteration over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all ghost cells in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : elements(gv,Dune::Partitions::ghost))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \sa elements(const GV&)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the elements.
   * \param ps   a PartitionSet object that is used to specify the set of Dune::PartitionType to which
   *             the elements must belong.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> elements(const GV& gv, PartitionSet<partitions> ps);


  //! Iterates over all facets (entities with codimension 1) of a GridView that belong to the given PartitionSet.
  /**
   * This functions returns an object representing the range of *facets* contained
   * in the GridView gv which belong to the PartitionSet ps. The main purpose of this function
   * is to enable iteration over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all interior and border facets in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : facets(gv,Dune::Partitions::interiorBorder))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \sa facets(const GV&)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the facets.
   * \param ps   a PartitionSet object that is used to specify the set of Dune::PartitionType to which
   *             the facets must belong.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  /**
   * \relates GridView
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> facets(const GV& gv, PartitionSet<partitions> ps);


  //! Iterates over all edges (entities with dimension 1) of a GridView that belong to the given PartitionSet.
  /**
   * This functions returns an object representing the range of *edges* contained
   * in the GridView gv which belong to the PartitionSet ps. The main purpose of this function
   * is to enable iteration over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all interior edges in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : edges(gv,Dune::Partitions::interior))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \sa edges(const GV&)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the edges.
   * \param ps   a PartitionSet object that is used to specify the set of Dune::PartitionType to which
   *             the edges must belong.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> edges(const GV& gv, PartitionSet<partitions> ps);


  //! Iterates over all vertices (entities with dimension 0) of a GridView that belong to the given PartitionSet.
  /**
   * This functions returns an object representing the range of *vertices* contained
   * in the GridView gv which belong to the PartitionSet ps. The main purpose of this function
   * is to enable iteration over those entities by means of a range-based for loop.
   *
   * Example:
   *
   * \code
   * // iterate over all interior vertices in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : vertices(gv,Dune::Partitions::interior))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \sa vertices(const GV&)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the vertices.
   * \param ps   a PartitionSet object that is used to specify the set of Dune::PartitionType to which
   *             the vertices must belong.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, unsigned int partitions>
  IteratorRange<...> vertices(const GV& gv, PartitionSet<partitions> ps);

  //! \}



  //! \name Generic entity ranges for non-standard parallel partitions
  //! \brief These Entity ranges allow for the maximum flexibility; they are parameterized on both the co(cimension) and the parallel PartitionSet.
  //! \{


  //! Iterates over all entities of a GridView with the given codimension that belong to the given PartitionSet.
  /**
   * This functions returns an object representing the range of entities of codimension cd contained
   * in the GridView gv which belong to the PartitionSet ps. The main purpose of this function is to
   * enable iteration over those entities by means of a range-based for loop:
   *
   * \code
   * // iterate over all interior and border cells in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : entities(gv,Dune::Codim<0>(),Dune::Partitions::interiorBorder))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This function allows you to write loops that are parameterized on the codimension of
   *         the entities. While this allows for extra flexibility (for e.g. some codimension-agnostic
   *         algorithms), it reduces the code readibility. If you don't need this flexibility, consider
   *         using elements(), facets(), edges() or vertices() instead.
   *
   * \remark If you have to iterate over entities with a specific *dimension*, consider using
   *         entities(const GV&,Dim<dim>,PartitionSet<partitions>) instead to improve the readability
   *         of your code.
   *
   * \sa entities(const GV&,Codim<codim>)
   * \sa entities(const GV&,Dim<dim>,PartitionSet<partitions>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the entities.
   * \param cd   a Codim object that is used to specify the codimension of the entities by means
   *             of its template parameter.
   * \param ps   a PartitionSet object that is used to specify the set of Dune::PartitionType to which
   *             the entities must belong.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, int codim, unsigned int partitions>
  IteratorRange<...> entities(const GV& gv gv, Codim<codim> cd, PartitionSet<partitions> ps);


  //! Iterates over all entities of a GridView with the given dimension that belong to the given PartitionSet.
  /**
   * This functions returns an object representing the range of entities of dimension d contained
   * in the GridView gv which belong to the PartitionSet ps. The main purpose of this function is to
   * enable iteration over those entities by means of a range-based for loop:
   *
   * \code
   * // iterate over all interior and border edges in the LeafGridView
   * auto gv = grid.leafGridView();
   * for (auto&& e : entities(gv,Dune::Dim<1>(),Dune::Partitions::interiorBorder))
   * {
   *   std::cout << e.geometry().center() << std::endl;
   * }
   * \endcode
   *
   * \remark This function allows you to write loops that are parameterized on the dimension of
   *         the entities. While this allows for extra flexibility (for e.g. some codimension-agnostic
   *         algorithms), it reduces the code readibility. If you don't need this flexibility, consider
   *         using elements(), facets(), edges() or vertices() instead.
   *
   * \remark If you have to iterate over entities with a specific *codimension*, consider using
   *         entities(const GV&,Codim<codim>,PartitionSet<partitions>) instead to improve the readability
   *         of your code.
   *
   * \sa entities(const GV&,Dim<dim>)
   * \sa entities(const GV&,Codim<codim>,PartitionSet<partitions>)
   *
   * \relates    GridView
   * \param gv   a GridView object that contains the entities.
   * \param cd   a Codim object that is used to specify the dimension of the entities by means
   *             of its template parameter.
   * \param ps   a PartitionSet object that is used to specify the set of Dune::PartitionType to which
   *             the entities must belong.
   * \returns    an unspecified object that is guaranteed to fulfil the interface
   *             of IteratorRange and that can be iterated over using a range-based
   *             for loop.
   */
  template<typename GV, int dim, unsigned int partitions>
  IteratorRange<...> entities(const GV& gv, Dim<dim> d, PartitionSet<partitions> ps);


  //! \}


#else // DOXYGEN


  // ******************************************************************************************
  // Implementations
  // ******************************************************************************************


  /**
   * Master entity range implementation - all other functions that return an entity range eventually
   * end up calling this function. The return type of this function is forwarded by those other functions
   * with the help of decltype.
   *
   * It is thus possible for a grid to use a different type of iterator range object by overloading
   * or partially specializing this single function. The specialization has to be placed either in the
   * Dune namespace or in the namespace of the GridView object.
   */
  template<typename GV, int codim, unsigned int partitions>
  IteratorRange<
    typename GV::template Codim<codim>::template Partition<
      derive_partition_iterator_type<partitions>::value
      >::Iterator
    >
  entities(const GV& gv, Codim<codim>, PartitionSet<partitions>)
  {
    static_assert(0 <= codim && codim <= GV::dimension, "invalid codimension for given GridView");
    const PartitionIteratorType pit = derive_partition_iterator_type<partitions>::value;
    typedef IteratorRange<
      typename GV::template Codim<codim>::template Partition<pit>::Iterator
      > return_type;
    return return_type(gv.template begin<codim,pit>(),gv.template end<codim,pit>());
  }

  /**
   * Entity range implementation without PartitionSet parameter. The default implementation forwards
   * this by adding Dune::Partitions::all, but if you have a GridView-like object that cannot handle
   * parallel partitions correctly, you can provide an overload for this function.
   *
   * All other functions without PartitionSet parameter forward to this function.
   */
  template<typename GV, int codim>
  auto entities(const GV& gv, Codim<codim> cd)
    -> decltype(entities(gv,cd,Partitions::all))
  {
    static_assert(0 <= codim && codim <= GV::dimension, "invalid codimension for given GridView");
    return entities(gv,cd,Partitions::all);
  }


  /**
   * Intersection range implementation.
   */
  template<typename GV, typename Entity>
  IteratorRange<typename GV::IntersectionIterator> intersections(const GV& gv, const Entity& e)
  {
    return IteratorRange<typename GV::IntersectionIterator>(gv.ibegin(e),gv.iend(e));
  }


  /**
   * Remaining implementations - these are mostly copy and paste and making sure to forward to the
   * correct function.
   */

  template<typename GV, int dim, unsigned int partitions>
  auto entities(const GV& gv, Dim<dim>, PartitionSet<partitions>)
    -> decltype(entities(gv,Codim<GV::dimension - dim>(),PartitionSet<partitions>()))
  {
    static_assert(0 <= dim && dim <= GV::dimension, "invalid dimension for given GridView");
    return entities(gv,Codim<GV::dimension - dim>(),PartitionSet<partitions>());
  }

  template<typename GV, int dim>
  auto entities(const GV& gv, Dim<dim>)
    -> decltype(entities(gv,Codim<GV::dimension - dim>()))
  {
    static_assert(0 <= dim && dim <= GV::dimension, "invalid dimension for given GridView");
    return entities(gv,Codim<GV::dimension - dim>());
  }

  template<typename GV, unsigned int partitions>
  auto elements(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Codim<0>(),PartitionSet<partitions>()))
  {
    return entities(gv,Codim<0>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto elements(const GV& gv)
    -> decltype(entities(gv,Codim<0>()))
  {
    return entities(gv,Codim<0>());
  }

  template<typename GV, unsigned int partitions>
  auto facets(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Codim<1>(),PartitionSet<partitions>()))
  {
    return entities(gv,Codim<1>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto facets(const GV& gv)
    -> decltype(entities(gv,Codim<1>()))
  {
    return entities(gv,Codim<1>());
  }

  template<typename GV, unsigned int partitions>
  auto edges(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Dim<1>(),PartitionSet<partitions>()))
  {
    return entities(gv,Dim<1>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto edges(const GV& gv)
    -> decltype(entities(gv,Dim<1>()))
  {
    return entities(gv,Dim<1>());
  }

  template<typename GV, unsigned int partitions>
  auto vertices(const GV& gv, PartitionSet<partitions>)
    -> decltype(entities(gv,Dim<0>(),PartitionSet<partitions>()))
  {
    return entities(gv,Dim<0>(),PartitionSet<partitions>());
  }

  template<typename GV>
  auto vertices(const GV& gv)
    -> decltype(entities(gv,Dim<0>()))
  {
    return entities(gv,Dim<0>());
  }


#endif // DOXYGEN

  /**
   * \} // GIIteration
   */

} // namespace Dune

#endif // DUNE_GRID_COMMON_ENTITYRANGEGENERATORS_HH
