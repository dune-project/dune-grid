// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_HOST_GRID_ACCESS_HH
#define DUNE_GEOGRID_HOST_GRID_ACCESS_HH

#include <string>

namespace Dune
{
  // HostGridAccess
  // ------------

  /** \class HostGridAccess
   *  \brief provides access to HostGrid Entities and EntityPointers
   *
   *  \tparam GeometryGrid
   *
   *  \nosubgrouping
   */

  template< class GeometryGrid >
  class HostGridAccess
  {

    typedef HostGridAccess< GeometryGrid > ThisType;

  public:

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
      //! type of HostEntity
      typedef typename GeometryGrid :: HostGrid :: template Codim< codim > :: Entity HostEntity;
      //! type of HostEntityPointer
      typedef typename GeometryGrid :: HostGrid :: template Codim< codim > :: EntityPointer HostEntityPointer;
    };

    /** \brief Get underlying HostGrid.
     *  \param[in] grid  GeometryGrid
     *  \returns HostGrid
     */
    static const HostGrid & hostGrid (const GeometryGrid & grid)
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
    template < int codim >
    static const typename ThisType :: template Codim< codim > :: HostEntity &
    getHostEntity ( const GeometryGrid & grid,
                    const typename GeometryGrid :: template Codim< codim > :: Entity & entity
                    )
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
    template < int codim >
    static const typename ThisType :: template Codim< codim > :: HostEntityPointer &
    getHostEntityPointer ( const GeometryGrid & grid,
                           const typename GeometryGrid :: template Codim< codim > :: EntityPointer & entityPointer
                           )
    {
      return grid.template getHostEntityPointer< codim > ( entityPointer );
    }

  };

}

#endif
