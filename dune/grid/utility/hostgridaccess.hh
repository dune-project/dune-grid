// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_HOSTGRIDACCESS_HH
#define DUNE_GRID_HOSTGRIDACCESS_HH

#include <string>

#include <dune/grid/geometrygrid/intersection.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  class GeometryGrid;

  template< class >
  class IdentityGrid;


  // HostGridAccess
  // --------------

  template< class Grid >
  struct HostGridAccess;



  /** \class HostGridAccess
   *  \brief provides access to host grid objects from GeometryGrid
   *
   *  \tparam GeometryGrid
   *
   *  \nosubgrouping
   */
  template< class HG, class CoordFunction, class Allocator >
  struct HostGridAccess< GeometryGrid< HG, CoordFunction, Allocator > >
  {
    /** \name Exported Types
     * \{ */

    typedef GeometryGrid< HG, CoordFunction, Allocator > Grid;

    //! type of HostGrid
    typedef typename Grid::HostGrid HostGrid;

    /** \} */

    /** \brief A Traits struct that collects return types of class member methods.
     *
     *  \tparam codim codimension
     */
    template< int codim >
    struct Codim
    {
      //! type of the GeometryGrid entity
      typedef typename Grid::template Codim< codim >::Entity Entity;

      //! type of the host entity
      typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
    };

    //! type of the GeometryGrid leaf intersection
    typedef typename Grid::Traits::LeafIntersection LeafIntersection;
    //! type of the GeometryGrid level intersection
    typedef typename Grid::Traits::LevelIntersection LevelIntersection;

    //! type of the host leaf intersection
    typedef typename HostGrid::Traits::LeafIntersection HostLeafIntersection;
    //! type of the host level intersection
    typedef typename HostGrid::Traits::LevelIntersection HostLevelIntersection;

    /** \brief Get underlying HostGrid.
     *  \param[in] grid  GeometryGrid
     *  \returns HostGrid
     */
    static const HostGrid &hostGrid ( const Grid &grid )
    {
      return grid.hostGrid();
    }

    template< class Entity >
    static const typename Codim< Entity::codimension >::HostEntity &
    hostEntity ( const Entity &entity )
    {
      return hostEntity< Entity::codimension >( entity );
    }

    template< int codim >
    static const typename Codim< codim >::HostEntity &
    hostEntity ( const typename Codim< codim >::Entity &entity )
    {
      return entity.impl().hostEntity();
    }

    template< class HostIntersection >
    static const HostIntersection &
    hostIntersection ( const Intersection< const Grid, GeoGrid::Intersection< const Grid, HostIntersection > > &intersection )
    {
      return intersection.impl().hostIntersection();
    }
  };


  /** \class HostGridAccess
   *  \brief provides access to host grid objects from IdentityGrid
   *
   *  \tparam IdentityGrid  meta grid, whose host grid shall be accessed
   *
   *  \nosubgrouping
   */
  template< class HG >
  struct HostGridAccess< IdentityGrid< HG > >
  {
    /** \name Exported Types
     * \{ */

    typedef IdentityGrid< HG > Grid;

    //! type of HostGrid
    typedef HG HostGrid;

    /** \} */

    /** \brief A Traits struct that collects return types of class member methods.
     *
     *  \tparam codim codimension
     */
    template< int codim >
    struct Codim
    {
      //! type of the IdGrid entity
      typedef typename Grid::template Codim< codim >::Entity Entity;

      //! type of the host entity
      typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
    };

    //! type of the IdGrid leaf intersection
    typedef typename Grid::Traits::LeafIntersection LeafIntersection;
    //! type of the IdGrid level intersection
    typedef typename Grid::Traits::LevelIntersection LevelIntersection;

    //! type of the host leaf intersection
    typedef typename HostGrid::Traits::LeafIntersection HostLeafIntersection;
    //! type of the host level intersection
    typedef typename HostGrid::Traits::LevelIntersection HostLevelIntersection;

    /** \brief Get underlying HostGrid.
     *  \param[in]  grid  grid, whose host grid shall be returned
     *  \returns HostGrid
     */
    static const HostGrid &hostGrid ( const Grid &grid )
    {
      return grid.getHostGrid();
    }

    template< class Entity >
    static const typename Codim< Entity::codimension >::HostEntity &
    hostEntity ( const Entity &entity )
    {
      return hostEntity< Entity::codimension >( entity );
    }

    template< int codim >
    static const typename Codim< codim >::HostEntity &
    hostEntity ( const typename Codim< codim >::Entity &entity )
    {
      return *entity.impl().hostEntity_;
    }

    static const HostLeafIntersection &
    hostIntersection ( const LeafIntersection &intersection )
    {
      return *intersection.impl().hostIterator_;
    }

    static const HostLevelIntersection &
    hostIntersection ( const LevelIntersection &intersection )
    {
      return *intersection.impl().hostIterator_;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_HOSTGRIDACCESS_HH
