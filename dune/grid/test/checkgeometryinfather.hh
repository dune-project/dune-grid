// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKGEOMETRYINFATHER_HH
#define DUNE_GRID_TEST_CHECKGEOMETRYINFATHER_HH

#include <dune/common/typetraits.hh>

#include "checkgeometry.hh"

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
  using std::sqrt;

  using ctype = typename GridType::ctype;

  // define a tolerance for floating-point checks
  const ctype tolerance = sqrt(std::numeric_limits< ctype >::epsilon());

  // count the number of different vertices
  unsigned int differentVertexCoords = 0;

  typedef typename GridType::Traits::LocalIdSet IdSet;

  const IdSet &idSet = grid.localIdSet();

  // We need at least two levels to do any checking
  if (grid.maxLevel()==0)
  {
    dwarn << "SKIPPING check geometryInFather(), because grid has only one level! \n";
    return;
  }

  // Loop over all levels except the lowest one
  for (int level = 1; level <= grid.maxLevel(); ++level)
  {
    typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
    ElementIterator eIt    = grid.levelGridView(level).template begin<0>();
    ElementIterator eEndIt = grid.levelGridView(level).template end<0>();

    for (; eIt!=eEndIt; ++eIt)
    {
      // check geometry
      checkGeometry( eIt->geometry() );

      // check the father method
      if( eIt->hasFather() )
      {
        typedef typename GridType::template Codim<0>::Entity Entity;
        Entity father = eIt->father();

        // check geometry
        checkGeometry( father.geometry() );

        while ( father.hasFather() )
        {
          if( father.level() == 0 )
            DUNE_THROW( GridError, "A level zero entity returns hasFather()=true." );
          Entity grandPa = father.father();
          typedef typename GridType :: Traits :: HierarchicIterator HierarchicIterator;

          const int mxl = grandPa.level() + 1;

          bool foundChild = false;
          const HierarchicIterator end = grandPa.hend( mxl );
          for( HierarchicIterator sons = grandPa.hbegin( mxl ); sons != end; ++sons )
          {
            // check geometry
            checkGeometry( sons->geometry() );

            if( father != *sons )
            {
              if( idSet.id( father ) == idSet.id( *sons ) )
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
            typedef typename GridType::template Codim<0>::Entity Entity;
            Entity father = sons->father();
            --count;
            while ( father.hasFather() )
            {
              father = father.father();
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
        checkLocalGeometry( geometryInFather, eIt->father().type(), "geometryInFather" );

        // //////////////////////////////////////////////////////
        //   Check for types and constants
        // //////////////////////////////////////////////////////

        static_assert((std::is_same<
                         typename Geometry::ctype,
                         typename GridType::ctype>::value == true),
                      "Geometry has wrong ctype");

        static_assert((static_cast<int>(Geometry::mydimension)
                       == static_cast<int>(GridType::dimension)),
                      "Geometry has wrong mydimension");

        static_assert((static_cast<int>(Geometry::coorddimension)
                       == static_cast<int>(GridType::dimensionworld)),
                      "Geometry has wrong coorddimension");

        // ///////////////////////////////////////////////////////
        //   Check the different methods
        // ///////////////////////////////////////////////////////
        if (geometryInFather.type() != eIt->type())
          DUNE_THROW(GridError, "Type of geometry and geometryInFather differ!");

        if (geometryInFather.corners() != eIt->geometry().corners())
          DUNE_THROW(GridError, "entity and geometryInFather have different number of corners!");

        // Compute the element center just to have an argument for the following methods
        auto refElement = referenceElement( geometryInFather );

        typename LocalGeometry::GlobalCoordinate center = refElement.position(0,0);

        if (geometryInFather.integrationElement(center) <=0)
          DUNE_THROW(GridError, "nonpositive integration element found!");

        // /////////////////////////////////////////////////////////////////////////////////////
        // check global transformation of geometryInFather
        // /////////////////////////////////////////////////////////////////////////////////////
        for( int j=0; j < eIt->geometry().corners(); ++j )
        {
          // create a global vertex
          const typename Geometry::GlobalCoordinate cornerViaSon =
            eIt->geometry().corner(j);
          // map to child local
          const typename Geometry::LocalCoordinate cornerInSon =
            eIt->geometry().local(cornerViaSon);
          // map to father
          const typename Geometry::LocalCoordinate cornerInFather =
            geometryInFather.global(cornerInSon);
          // map father to global
          const typename Geometry::GlobalCoordinate cornerViaFather =
            eIt->father().geometry().global(cornerInFather);

          if( (cornerViaFather - cornerViaSon).infinity_norm() > tolerance )
          {
            ++differentVertexCoords;
            std :: cout << "global transformation of geometryInFather yields different vertex position "
                        << "(son: " << cornerViaSon
                        << ", father: " << cornerViaFather << ")." << std :: endl;
          }
        }

        const typename Geometry::LocalCoordinate X(0.2);
        typename Geometry::LocalCoordinate x = X;

        typename GridType::template Codim<0>::Entity e(*eIt);
        while (e.level() != 0)
        {
            x = e.geometryInFather().global(x);
            e = e.father();
        }

        if ((e.geometry().global(x)-eIt->geometry().global(X)).two_norm() > tolerance)
        {
          std::cerr << "Warning: mapping broken! " << e.geometry().global(x)
                    << " vs. "  << eIt->geometry().global(X)
                    << "\tchild " << eIt->geometry().center()
                    << "\tfather " << e.geometry().center()
                    << "\tat "  << x
                    << "\tmaps to " << X << std::endl;
        }

        for( int j=0; j < eIt->geometry().corners(); ++j )
        {
          // create a global vertex
          const typename Geometry::GlobalCoordinate cornerViaSon =
            eIt->geometry().corner(j);
          // map to child local
          const typename Geometry::LocalCoordinate cornerInSon =
            eIt->geometry().local(cornerViaSon);
          // map to father local
          const typename Geometry::LocalCoordinate cornerInFather =
            geometryInFather.global(cornerInSon);
          // map father to
          const typename Geometry::LocalCoordinate cornerFromGlobal =
            eIt->father().geometry().local(cornerViaSon);

          if( (cornerInFather - cornerFromGlobal).infinity_norm() > tolerance )
          {
            ++differentVertexCoords;
            std :: cout << "geometryInFather().global() and yields different vertex position than father().geometry().local()"
                        << "(from son: " << cornerInFather
                        << ", from global: " << cornerFromGlobal << ")." << std :: endl;
            break;
          }
        }

        // /////////////////////////////////////////////////////////////////////////////////////
        // check local transformation of geometryInFather
        // /////////////////////////////////////////////////////////////////////////////////////
        for( int j=0; j < eIt->geometry().corners(); ++j )
        {
          // create a global vertex
          const typename Geometry::GlobalCoordinate global =
            eIt->geometry().corner(j);
          // map to child local
          const typename Geometry::LocalCoordinate cornerInSon =
            eIt->geometry().local(global);
          // map global to father
          const typename Geometry::LocalCoordinate cornerInFather =
            eIt->father().geometry().local(global);

          // map from father to son
          const typename Geometry::LocalCoordinate cornerViaFather =
            geometryInFather.local(cornerInFather);

          if( (cornerViaFather - cornerInSon).infinity_norm() > tolerance )
          {
            ++differentVertexCoords;
            std :: cout << "local transformation of geometryInFather yields different vertex position "
                        << "(global: " << global
                        << ", son: " << cornerInSon
                        << ", father: " << cornerViaFather << ")." << std :: endl;
          }
        }

        /** \todo Missing jacobianInverse() */
        /** \todo Missing checkInside() */

        // /////////////////////////////////////////////////////////////////////////////////////
        // Check whether the positions of the vertices of geometryInFather coincide
        // with the ones computed 'by hand'.  This only works if the grids really are nested!
        // /////////////////////////////////////////////////////////////////////////////////////
        for( int j=0; j < geometryInFather.corners(); ++j )
        {
          const typename Geometry::GlobalCoordinate cornerViaFather
            = eIt->father().geometry().global( geometryInFather.corner( j ) );

          const typename Geometry::GlobalCoordinate &cornerViaSon = eIt->geometry().corner( j );

          if( (cornerViaFather - cornerViaSon).infinity_norm() > tolerance )
          {
            ++differentVertexCoords;
            std :: cout << "geometryInFather yields different vertex position "
                        << "(son: " << cornerViaSon
                        << ", father: " << cornerViaFather << ")." << std :: endl;
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

#endif // DUNE_GRID_TEST_CHECKGEOMETRYINFATHER_HH
