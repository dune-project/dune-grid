// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CHECK_GEOMETRYINFATHER_CC
#define DUNE_CHECK_GEOMETRYINFATHER_CC

#include <dune/common/typetraits.hh>

#include "checkgeometry.cc"

/** \file
    \brief A test for the Method Geometry::geometryInFather()
 */

/** \brief Test the Method Geometry::geometryInFather()

   This test works by comparing the output of geometryInFather with the vertex positions
   obtained by directly expressing the son vertex positions in local coordinates of the
   father.  That should work for all grid implementations that realize truly nested
   grids.  One exception is UGGrid with parametrized boundaries.
 */
template <class GridType>
void checkGeometryInFather(const GridType& grid)
{
  using namespace Dune;

  // count the number of different vertices
  unsigned int differentVertexCoords = 0;

  typedef typename GridType::Traits::LocalIdSet IdSet;
  typedef typename GridType::ctype ctype;

  const IdSet &idSet = grid.localIdSet();

  // We need at least two levels to do any checking
  if (grid.maxLevel()==0)
  {
    dwarn << "SKIPPING check geometryInFather(), because grid has only one level! \n";
    return;
  }

  // Loop over all levels except the lowest one
  for (int i=1; i<=grid.maxLevel(); i++) {

    typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
    ElementIterator eIt    = grid.template lbegin<0>(i);
    ElementIterator eEndIt = grid.template lend<0>(i);

    for (; eIt!=eEndIt; ++eIt)
    {
      // check geometry
      checkGeometry( eIt->geometry() );

      // check the father method
      if( eIt->hasFather() )
      {
        typedef typename GridType::template Codim<0>::EntityPointer EntityPointer;
        EntityPointer father = eIt->father();

        // check geometry
        checkGeometry( father->geometry() );

        while ( father->hasFather() )
        {
          if( father->level() == 0 )
            DUNE_THROW( GridError, "A level zero entity returns hasFather()=true." );
          EntityPointer grandPa = father->father();
          typedef typename GridType :: Traits :: HierarchicIterator HierarchicIterator;

          const int mxl = grandPa->level() + 1;

          bool foundChild = false;
          const HierarchicIterator end = grandPa->hend( mxl );
          for( HierarchicIterator sons = grandPa->hbegin( mxl ); sons != end; ++sons )
          {
            // check geometry
            checkGeometry( sons->geometry() );

            if( father != sons )
            {
              if( idSet.id( *father ) == idSet.id( *sons ) )
                DUNE_THROW( GridError, "Two different entities have the same id." );
            }
            else
              foundChild = true;
          }

          if( !foundChild )
            DUNE_THROW( GridError, "Cannot find child in its own father." );
          father = grandPa;
        }
      }

      // hierarchy check
      {
        typedef typename GridType :: Traits :: HierarchicIterator HierarchicIterator;
        const int mxl = grid.maxLevel();
        int countChildren = 0;

        HierarchicIterator end = eIt->hend(mxl);
        for(HierarchicIterator sons = eIt->hbegin(mxl);
            sons != end; ++sons)
        {
          ++countChildren;
          int count = sons->level();
          if( sons->hasFather() )
          {
            typedef typename GridType::template Codim<0>::EntityPointer EntityPointer;
            EntityPointer father = sons->father();
            --count;
            while ( father->hasFather() )
            {
              father = father->father();
              --count;
            }
          }
          assert( count == 0 );
        }

        if( eIt->isLeaf () && countChildren > 0 )
          DUNE_THROW(GridError, "leaf entity has children ==> entity is not leaf");
      }

      if ( eIt->hasFather() )
      {
        // Get geometry in father
        typedef typename GridType::template Codim<0>::Entity::Geometry Geometry;
        typedef typename GridType::template Codim<0>::Entity::LocalGeometry LocalGeometry;

        const LocalGeometry& geometryInFather = eIt->geometryInFather();
        checkLocalGeometry( geometryInFather, eIt->father()->type(), "geometryInFather" );

        // //////////////////////////////////////////////////////
        //   Check for types and constants
        // //////////////////////////////////////////////////////

        dune_static_assert((is_same<
                                typename Geometry::ctype,
                                typename GridType::ctype>::value == true),"Geometry has wrong ctype");

        dune_static_assert((static_cast<int>(Geometry::dimension)
                            == static_cast<int>(GridType::dimension)),"Geometry has wrong dimension");

        dune_static_assert((static_cast<int>(Geometry::mydimension)
                            == static_cast<int>(GridType::dimension)),"Geometry has wrong mydimension");

        dune_static_assert((static_cast<int>(Geometry::coorddimension)
                            == static_cast<int>(GridType::dimensionworld)),"Geometry has wrong coorddimension");

        dune_static_assert((static_cast<int>(Geometry::dimensionworld)
                            == static_cast<int>(GridType::dimensionworld)),"Geometry has wrong dimensionworld");

        // ///////////////////////////////////////////////////////
        //   Check the different methods
        // ///////////////////////////////////////////////////////
        if (geometryInFather.type() != eIt->type())
          DUNE_THROW(GridError, "Type of geometry and geometryInFather differ!");

        if (geometryInFather.corners() != eIt->geometry().corners())
          DUNE_THROW(GridError, "entity and geometryInFather have different number of corners!");

        // Compute the element center just to have an argument for the following methods
        typename LocalGeometry::GlobalCoordinate center(0);

        for (int j=0; j<geometryInFather.corners(); j++)
          center += geometryInFather.corner( j );

        if (geometryInFather.integrationElement(center) <=0)
          DUNE_THROW(GridError, "nonpositive integration element found!");

        /** \todo Missing local() */
        /** \todo Missing global() */
        /** \todo Missing jacobianInverse() */
        /** \todo Missing checkInside() */

        // /////////////////////////////////////////////////////////////////////////////////////
        // Check whether the positions of the vertices of geometryInFather coincide
        // with the ones computed 'by hand'.  This only works if the grids really are nested!
        // /////////////////////////////////////////////////////////////////////////////////////
        for( int j=0; j < geometryInFather.corners(); ++j )
        {
          const typename Geometry::GlobalCoordinate cornerInFather
            = eIt->father()->geometry().global( geometryInFather.corner( j ) );
          const typename Geometry::GlobalCoordinate &cornerInSon = eIt->geometry().corner( j );

          if( (cornerInFather - cornerInSon).infinity_norm() > 1e-7 )
          {
            ++differentVertexCoords;
            std :: cout << "geometryInFather yields different vertex position "
                        << "(son: " << cornerInSon
                        << ", father: " << cornerInFather << ")." << std :: endl;
          }
        }
      }
    }

  }

  if( differentVertexCoords > 0 )
  {
    std :: cerr << "Warning: geometryInFather yields different vertex positions." << std :: endl;
    std :: cerr << "         This behaviour may be correct if the grid is not"
                << " nested geometrically." << std :: endl;
  }

}

#endif
