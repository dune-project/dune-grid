// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_HOSTGRIDACCESS_HH
#define DUNE_GEOGRID_HOSTGRIDACCESS_HH

#include <string>

namespace Dune
{

  // HostGridAccess
  // --------------

  /** \class HostGridAccess
   *  \brief provides access to host grid objects
   *
   *  \tparam GeometryGrid
   *
   *  \nosubgrouping
   */
  template< class GeometryGrid >
  struct HostGridAccess
  {
    /** \name Exported Types
     * \{ */

    //! type of HostGrid
    typedef typename GeometryGrid :: HostGrid HostGrid;

    /** \} */

    /** \brief A Traits struct that collects return types of class member methods.
     *
     *  \tparam codim codimension
     */
    template< int codim >
    struct Codim
    {
      //! type of the GeometryGrid entity
      typedef typename GeometryGrid :: template Codim< codim > :: Entity Entity;
      //! type of the GeometryGrid entity pointer
      typedef typename GeometryGrid :: template Codim< codim > :: EntityPointer
      EntityPointer;

      //! type of the host entity
      typedef typename HostGrid :: template Codim< codim > :: Entity HostEntity;
      //! type of the host entity pointer
      typedef typename HostGrid :: template Codim< codim > :: EntityPointer
      HostEntityPointer;
    };

    //! type of the GeometryGrid leaf intersection
    typedef typename GeometryGrid :: Traits :: LeafIntersection LeafIntersection;
    //! type of the GeometryGrid level intersection
    typedef typename GeometryGrid :: Traits :: LevelIntersection LevelIntersection;

    //! type of the host leaf intersection
    typedef typename HostGrid :: Traits :: LeafIntersection HostLeafIntersection;
    //! type of the host level intersection
    typedef typename HostGrid :: Traits :: LevelIntersection HostLevelIntersection;

    /** \brief Get underlying HostGrid.
     *  \param[in] grid  GeometryGrid
     *  \returns HostGrid
     */
    static const HostGrid &hostGrid ( const GeometryGrid &grid )
    {
      return grid.hostGrid();
    }

    /** \brief Get underlying HostEntity of given GeometryGridEntity.
     *
     *  Return type is exported by the Codim struct.
     *  \param[in] grid  GeometryGrid
     *  \param[in] entity  GeometryGridEntity
     *  \returns HostGridEntity
     */
    template< int codim >
    static const typename Codim< codim > :: HostEntity &
    getHostEntity ( const GeometryGrid &grid,
                    const typename Codim< codim > :: Entity &entity )
    {
      return grid.template getHostEntity< codim > ( entity );
    }

    /** \brief Get underlying HostEntityPointer of given GeometryGridEntityPointer.
     *
     *  Return type is exported by the Codim struct.
     *  \param[in] grid  GeometryGrid
     *  \param[in] entityPointer  GeometryGridEntityPointer
     *  \returns HostEntityPointer
     */
    template< int codim >
    static const typename Codim< codim > :: HostEntityPointer &
    getHostEntityPointer ( const GeometryGrid &grid,
                           const typename Codim< codim > :: EntityPointer &entityPointer )
    {
      return grid.template getHostEntityPointer< codim > ( entityPointer );
    }

    static const HostLeafIntersection &
    getIntersection ( const LeafIntersection &intersection )
    {
      return GeometryGrid :: getRealImplementation( intersection ).hostIntersection();
    }

    static const HostLevelIntersection &
    getIntersection ( const LevelIntersection &intersection )
    {
      return GeometryGrid :: getRealImplementation( intersection ).hostIntersection();
    }
  };

}

#endif
