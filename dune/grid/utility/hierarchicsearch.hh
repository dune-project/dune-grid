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
    enum {dim=Grid::dimension};

    //! get world dimension from the grid
    enum {dimw=Grid::dimensionworld};

    //! get coord type from the grid
    typedef typename Grid::ctype ct;

    //! get entity from the grid
    typedef typename Grid::template Codim<0>::Entity Entity;

    //! type of EntityPointer
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

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

       @param[in] entity Entity whos children should be searched
       @param[in] global Point you are searching for

       Search the child entity containing point global. Recursively
       recursively continue until we found an entity that is part of
       the IndexSet.
     */
    EntityPointer hFindEntity ( const Entity &entity,
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
        const Entity &child = *it;
        const Geometry &geo = child.geometry();

        LocalCoordinate local = geo.local(global);
        if (ReferenceElements<double, dim>::general( child.type() ).checkInside(local))
        {
          // return if we found the leaf, else search through the child entites
          if( indexSet_.contains( child ) )
            return EntityPointer( child );
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
    EntityPointer findEntity(const FieldVector<ct,dimw>& global) const
    { return findEntity<All_Partition>(global); }

    /**
       @brief Search the IndexSet of this HierarchicSearch for an Entity
       containing point global.

       \exception GridError No element of the coarse grid contains the given
                            coordinate.
     */
    template<PartitionIteratorType partition>
    EntityPointer findEntity(const FieldVector<ct,dimw>& global) const
    {
      typedef typename Grid::template Partition<partition>::LevelGridView
      LevelGV;
      const LevelGV &gv = grid_.template levelView<partition>(0);

      //! type of LevelIterator
      typedef typename LevelGV::template Codim<0>::Iterator LevelIterator;

      // type of element geometry
      typedef typename Entity::Geometry Geometry;
      // type of local coordinate
      typedef typename Geometry::LocalCoordinate LocalCoordinate;

      // loop over macro level
      LevelIterator it = gv.template begin<0>();
      LevelIterator end = gv.template end<0>();
      for (; it != end; ++it)
      {
        const Entity &entity = *it;
        const Geometry &geo = entity.geometry();

        LocalCoordinate local = geo.local( global );
        if( !ReferenceElements< double, dim >::general( geo.type() ).checkInside( local ) )
          continue;

        if( (int(dim) != int(dimw)) && ((geo.global( local ) - global).two_norm() > 1e-8) )
          continue;

        // return if we found the leaf, else search through the child entites
        if( indexSet_.contains( entity ) )
          return EntityPointer( entity );
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
