// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_GRIDFAMILTY_HH
#define DUNE_ALBERTAGRID_GRIDFAMILTY_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides the GridFamily for AlbertaGrid
 */

#include <dune/common/parallel/communication.hh>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/gridview.hh>
#include <dune/grid/albertagrid/geometryreference.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;

  template< int codim, int dim, class GridImp >
  class AlbertaGridEntity;

  template< int codim, class GridImp >
  class AlbertaGridEntityPointer;

  template< int codim, class Grid >
  class AlbertaGridEntitySeed;

  template< int mydim, int cdim, class GridImp >
  class AlbertaGridGeometry;

  template< int mydim, int cdim, class GridImp >
  class AlbertaGridGlobalGeometry;

  template< int dim, int dimworld >
  class AlbertaGridHierarchicIndexSet;

  template< class GridImp >
  class AlbertaGridHierarchicIterator;

  template< class GridImp >
  class AlbertaGridLeafIntersection;

  template< class GridImp >
  class AlbertaGridLeafIntersectionIterator;

  template< int dim, int dimworld >
  class AlbertaGridIdSet;

  template< int dim, int dimworld >
  class AlbertaGridIndexSet;

  template< int codim, PartitionIteratorType pitype, class GridImp >
  class AlbertaGridLeafIterator;

  template< int codim, PartitionIteratorType pitype, class GridImp >
  class AlbertaGridLevelIterator;




  // AlbertaGridFamily
  // -----------------

  template <int dim, int dimworld>
  struct AlbertaGridFamily
  {
    typedef AlbertaGrid< dim, dimworld > GridImp;

    typedef Alberta::Real ctype;

    static const int dimension = dim;
    static const int dimensionworld = dimworld;

    typedef AlbertaGridIndexSet< dim, dimworld > LevelIndexSetImp;
    typedef AlbertaGridIndexSet< dim, dimworld > LeafIndexSetImp;

    typedef AlbertaGridIdSet< dim, dimworld > IdSetImp;
    typedef unsigned int IdType;

    struct Traits
    {
      typedef GridImp Grid;

      typedef Dune::Intersection< const GridImp, AlbertaGridLeafIntersection< const GridImp > > LeafIntersection;
      typedef Dune::Intersection< const GridImp, AlbertaGridLeafIntersection< const GridImp > > LevelIntersection;
      typedef Dune::IntersectionIterator
      < const GridImp, AlbertaGridLeafIntersectionIterator< const GridImp >, AlbertaGridLeafIntersection< const GridImp > >
      LeafIntersectionIterator;
      typedef Dune::IntersectionIterator
      < const GridImp, AlbertaGridLeafIntersectionIterator< const GridImp >, AlbertaGridLeafIntersection< const GridImp > >
      LevelIntersectionIterator;

      typedef Dune::EntityIterator< 0, const GridImp, AlbertaGridHierarchicIterator< const GridImp > > HierarchicIterator;

      typedef IdType GlobalIdType;
      typedef IdType LocalIdType;

      template< int cd >
      struct Codim
      {
        typedef AlbertaGridGlobalGeometry< dim-cd, dimworld, const Grid > GeometryImpl;
        typedef AlbertaGridGeometry< dim-cd, dim, const Grid > LocalGeometryImpl;
        typedef Dune::Geometry< dim-cd, dimworld, const Grid, AlbertaGridGlobalGeometry > Geometry;
        typedef Dune::Geometry< dim-cd, dim, const Grid, LocalGeometryReference > LocalGeometry;

        typedef AlbertaGridEntity< cd, dim, const GridImp > EntityImpl;
        typedef Dune::Entity< cd, dim, const GridImp, AlbertaGridEntity > Entity;

        typedef Dune::EntitySeed< const GridImp, AlbertaGridEntitySeed< cd, const GridImp > > EntitySeed;

        template <PartitionIteratorType pitype>
        struct Partition
        {
          typedef Dune::EntityIterator< cd, const GridImp, AlbertaGridLevelIterator< cd, pitype, const GridImp > > LevelIterator;
          typedef Dune::EntityIterator< cd, const GridImp, AlbertaGridLeafIterator< cd, pitype, const GridImp > > LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
      };

      typedef Dune::GridView< AlbertaLevelGridViewTraits< const GridImp > > LevelGridView;
      typedef Dune::GridView< AlbertaLeafGridViewTraits< const GridImp > > LeafGridView;

      typedef IndexSet< GridImp, LevelIndexSetImp, int, std::array< GeometryType, 1 > > LevelIndexSet;
      typedef IndexSet< GridImp, LeafIndexSetImp, int, std::array< GeometryType, 1 > > LeafIndexSet;
      typedef AlbertaGridHierarchicIndexSet< dim, dimworld > HierarchicIndexSet;
      typedef IdSet<GridImp,IdSetImp,IdType> GlobalIdSet;
      typedef IdSet<GridImp,IdSetImp,IdType> LocalIdSet;

      typedef Dune::Communication< No_Comm > Communication;

      /**
       * \deprecated Use Communication instead! Will be removed after Dune 2.9.
       */
      [[deprecated("Use Communication instead!")]]
      typedef Communication CollectiveCommunication;
    };
  };

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTAGRID_GRIDFAMILTY_HH
