// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_HIERARCHICSEARCH_HH
#define DUNE_GRID_HIERARCHICSEARCH_HH

/**
   @file
   @brief Utility class for hierarchically searching for an Entity
   containing a given point.
 */

#include <cstddef>
#include <sstream>
#include <string>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /**
     @brief Search an IndexSet for an Entity containing a given point.
   */
  template<class Grid, class IS>
  class HierarchicSearch
  {
    //! get dimension from the grid
    constexpr static int dim = Grid::dimension;

    //! get world dimension from the grid
    constexpr static int dimw = Grid::dimensionworld;

    //! get coord type from the grid
    typedef typename Grid::ctype ct;

    //! get entity from the grid
    typedef typename Grid::template Codim<0>::Entity Entity;

    //! type of HierarchicIterator
    typedef typename Grid::HierarchicIterator HierarchicIterator;

    static std::string formatEntityInformation ( const Entity &e ) {
      const typename Entity::Geometry &geo = e.geometry();
      std::ostringstream info;
      info << "level=" << e.level() << " "
           << "partition=" << e.partitionType() << " "
           << "center=(" << geo.center() << ") "
           << "corners=[(" << geo.corner(0) << ")";
      for(int i = 1; i < geo.corners(); ++i)
        info << " (" << e.geometry().corner(i) << ")";
      info << "]";
      return info.str();
    }

    /**
       internal helper method

       @param[in] entity Entity whose children should be searched
       @param[in] global Point you are searching for

       Search the child entity containing point global. Recursively
       recursively continue until we found an entity that is part of
       the IndexSet.
     */
    Entity hFindEntity ( const Entity &entity,
                         const FieldVector<ct,dimw>& global) const
    {
      // type of element geometry
      typedef typename Entity::Geometry Geometry;
      // type of local coordinate
      typedef typename Geometry::LocalCoordinate LocalCoordinate;

      const int childLevel = entity.level()+1 ;
      // loop over all child Entities
      const HierarchicIterator end = entity.hend( childLevel );
      for( HierarchicIterator it = entity.hbegin( childLevel ); it != end; ++it )
      {
        Entity child = *it;
        Geometry geo = child.geometry();

        LocalCoordinate local = geo.local(global);
        if (referenceElement( geo ).checkInside(local))
        {
          // return if we found the leaf, else search through the child entites
          if( indexSet_.contains( child ) )
            return child;
          else
            return hFindEntity( child, global );
        }
      }
      std::ostringstream children;
      HierarchicIterator it = entity.hbegin( childLevel );
      if(it != end) {
        children << "{" << formatEntityInformation(*it) << "}";
        for( ++it; it != end; ++it )
          children << " {" << formatEntityInformation(*it) << "}";
      }
      DUNE_THROW(Exception, "{" << className(*this) << "} Unexpected "
                 "internal Error: none of the children of the entity "
                 "{" << formatEntityInformation(entity) << "} contains "
                 "coordinate (" << global << ").  Children are: "
                 "[" << children.str() << "].");
    }

  public:
    /**
       @brief Construct a HierarchicSearch object from a Grid and an IndexSet
     */
    HierarchicSearch(const Grid & g, const IS & is) : grid_(g), indexSet_(is) {}

    /**
       @brief Search the IndexSet of this HierarchicSearch for an Entity
       containing point global.

       \exception GridError No element of the coarse grid contains the given
                            coordinate.
     */
    Entity findEntity(const FieldVector<ct,dimw>& global) const
    { return findEntity<All_Partition>(global); }

    /**
       @brief Search the IndexSet of this HierarchicSearch for an Entity
       containing point global.

       \exception GridError No element of the coarse grid contains the given
                            coordinate.
     */
    template<PartitionIteratorType partition>
    Entity findEntity(const FieldVector<ct,dimw>& global) const
    {
      typedef typename Grid::LevelGridView LevelGV;
      const LevelGV &gv = grid_.levelGridView(0);

      //! type of LevelIterator
      typedef typename LevelGV::template Codim<0>::template Partition<partition>::Iterator LevelIterator;

      // type of element geometry
      typedef typename Entity::Geometry Geometry;
      // type of local coordinate
      typedef typename Geometry::LocalCoordinate LocalCoordinate;

      // loop over macro level
      const LevelIterator end = gv.template end<0, partition>();
      for (LevelIterator it = gv.template begin<0, partition>(); it != end; ++it)
      {
        Entity entity = *it;
        Geometry geo = entity.geometry();

        LocalCoordinate local = geo.local( global );
        if( !referenceElement( geo ).checkInside( local ) )
          continue;

        if( (int(dim) != int(dimw)) && ((geo.global( local ) - global).two_norm() > 1e-8) )
          continue;

        // return if we found the leaf, else search through the child entites
        if( indexSet_.contains( entity ) )
          return entity;
        else
          return hFindEntity( entity, global );
      }
      DUNE_THROW( GridError, "Coordinate " << global << " is outside the grid." );
    }

  private:
    const Grid& grid_;
    const IS&   indexSet_;
  };

} // end namespace Dune

#endif // DUNE_GRID_HIERARCHICSEARCH_HH
