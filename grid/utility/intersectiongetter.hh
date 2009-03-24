// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: p1groundwater.hh 336 2006-05-03 13:09:05Z oliver $
#ifndef DUNE_INTERSECTIONGETTER_HH
#define DUNE_INTERSECTIONGETTER_HH

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/deprecated.hh>

/**
 * @file
 * @author Markus Blatt
 * @brief Utilities to get the intersection iterator depending
 *
 * Depending on whether you work on a leaf or a level you you have to call
 * the methods ilevel{begin,end}() or ileaf{begin,end}() to get the correct iterator.
 *
 * @deprecated This helper class is obsolete when using the new GridView interface
 */
namespace Dune
{
  /**
   * @ingroup GIIntersectionIterator{
   */

  /**
   * @brief A tag to identify that we work on a level of a grid.
   */
  struct LevelTag {} DUNE_DEPRECATED;

  /**
   * @brief A tag to identify that we work on the leaf of a grid.
   */
  struct LeafTag {} DUNE_DEPRECATED;

  /**
   * @brief Utility class to get the Intersection Iterator the right way.
   *
   * Depending on whether one works on the leaf or a level of a grid,
   * the methods ilevel{begin,end}() or ileaf{begin,end}()have to called
   * to get get the correct iterator.
   * This can be done in a generic way using the class.
   *
   * The template parameters are:
   * <dt>Grid</dt><dd>The grid implementation.</dd>
   * <dt>Tag</dt><dd>The tag identifying whether we work on the leaf or a level
   * of the grid. Either LeafTag, or LevelTag.</dd>
   */
  template<typename Grid, typename Tag>
  struct IntersectionIteratorGetter
  {
    /**
     * @brief Tag telling whether we work on the level or leaf.
     *
     * This should be either LeafTag or LevelTag. If it is leaf
     * ::begin() and ::end() call ileafbegin() and ileafend() on
     * the entity or Codim<0>::Iterator.
     */
    typedef Tag TypeTag;

    /**
     * @brief The type of the Intersection Iterator
     */
    typedef typename Grid::template Codim<0>::IntersectionIterator IntersectionIterator;

    /**
     * @brief Get the correct begin Iterator depending on the TypeTag.
     * @param iter the entity or Codim<0>::Iterator
     */
    template<typename T>
    inline static IntersectionIterator begin(T& iter)
    {
      // Trigger a compile error
      IsTrue<is_same<TypeTag,LevelTag>::value||is_same<TypeTag,LevelTag>::value>::yes();
      return iter.ileafbegin();
    }
    /**
     * @brief Get the correct end Iterator depending on the TypeTag.
     * @param iter the entity or Codim<0>::Iterator
     */
    template<typename T>
    inline static IntersectionIterator end(T& iter)
    {
      // Trigger a compile error
      IsTrue<is_same<TypeTag,LevelTag>::value||is_same<TypeTag,LevelTag>::value>::yes();
      return iter.ileafend();
    }
  };

  // Specialization for the leaf
  template<typename Grid>
  struct IntersectionIteratorGetter<Grid,LeafTag>
  {
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

    template<typename T>
    inline static IntersectionIterator begin(T& iter)
    {
      return iter.ileafbegin();
    }

    template<typename T>
    inline static IntersectionIterator end(T& iter)
    {
      return iter.ileafend();
    }
  };

  // Specialization of the level
  template<typename Grid>
  struct IntersectionIteratorGetter<Grid,LevelTag>
  {
    typedef typename Grid::template Codim<0>::LevelIntersectionIterator IntersectionIterator;

    template<typename T>
    inline static IntersectionIterator begin(T& iter)
    {
      return iter.ilevelbegin();
    }

    template<typename T>
    inline static IntersectionIterator end(T& iter)
    {
      return iter.ilevelend();
    }
  };

  /**
   * @}
   */

} // end namespace Dune
#endif
