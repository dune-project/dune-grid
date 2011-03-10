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

       @param[in] e      EntityPointer whos children should be searched
       @param[in] global Point you are searching for

       Search the child entity containing point global. Recursively
       recursively continue until we found an entity that is part of
       the IndexSet.
     */
    EntityPointer hFindEntity ( const Entity &e,
                                const FieldVector<ct,dimw>& global) const
    {
      // loop over all child Entities
      const HierarchicIterator end = e.hend(e.level()+1);
      for( HierarchicIterator it = e.hbegin( e.level()+1 ); it != end; ++it )
      {
        FieldVector<ct,dim> local = it->geometry().local(global);
        if (GenericReferenceElements<double, dim>::general(it->type()).checkInside(local))
        {
          // return if we found the leaf, else search through the child entites
          if( is.contains( *it ) )
            return EntityPointer( it );
          else
            return hFindEntity( *it, global );
        }
      }
      std::ostringstream children;
      HierarchicIterator it = e.hbegin( e.level()+1 );
      if(it != end) {
        children << "{" << formatEntityInformation(*it) << "}";
        for( ++it; it != end; ++it )
          children << " {" << formatEntityInformation(*it) << "}";
      }
      DUNE_THROW(Exception, "{" << className(*this) << "} Unexpected "
                 "internal Error: none of the children of the entity "
                 "{" << formatEntityInformation(e) << "} contains "
                 "coordinate (" << global << ").  Children are: "
                 "[" << children.str() << "].");
    }

  public:
    /**
       @brief Construct a HierarchicSearch object from a Grid and an IndexSet
     */
    HierarchicSearch(const Grid & _g, const IS & _is) : g(_g), is(_is) {};

    /**
       @brief Search the IndexSet of this HierarchicSearch for an Entity
       containing point global.

       \exception GridError No element of the coarse grid contains the given
                            coordinate.
     */
    template<PartitionIteratorType partition = All_Partition>
    EntityPointer findEntity(const FieldVector<ct,dimw>& global) const
    {
      typedef typename Grid::template Partition<partition>::LevelGridView
      LevelGV;
      const LevelGV &gv = g.template levelView<partition>(0);

      //! type of LevelIterator
      typedef typename LevelGV::template Codim<0>::Iterator LevelIterator;

      // loop over macro level
      LevelIterator it = gv.template begin<0>();
      LevelIterator end = gv.template end<0>();
      for (; it != end; ++it)
      {
        const Entity &e = *it;
        const typename Entity::Geometry &geo = e.geometry();

        FieldVector< ct, dim > local = geo.local( global );
        if( !GenericReferenceElements< double, dim >::general( geo.type() ).checkInside( local ) )
          continue;

        if( (int(dim) != int(dimw)) && ((geo.global( local ) - global).two_norm() > 1e-8) )
          continue;

        // return if we found the leaf, else search through the child entites
        if( is.contains( *it ) )
          return EntityPointer( it );
        else
          return hFindEntity( *it, global );
      }
      DUNE_THROW( GridError, "Coordinate " << global << " is outside the grid." );
    }

  private:
    const Grid& g;
    const IS& is;
  };

} // end namespace Dune

#endif // DUNE_GRID_HIERARCHICSEARCH_HH
