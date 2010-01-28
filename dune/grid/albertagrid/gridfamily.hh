// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_GRIDFAMILTY_HH
#define DUNE_ALBERTAGRID_GRIDFAMILTY_HH

/** \file
 *  \author Martin Nolte
 *  \brief  provides the GridFamily for AlbertaGrid
 */

#include <dune/common/collectivecommunication.hh>

#include <dune/grid/common/defaultgridview.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/entitypointer.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/hierarchiciterator.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/common/leafiterator.hh>
#include <dune/grid/common/leveliterator.hh>

#include <dune/grid/albertagrid/misc.hh>

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

      typedef Dune::Intersection< const GridImp, AlbertaGridLeafIntersection > LeafIntersection;
      typedef Dune::Intersection< const GridImp, AlbertaGridLeafIntersection > LevelIntersection;
      typedef Dune::IntersectionIterator
      < const GridImp, AlbertaGridLeafIntersectionIterator, AlbertaGridLeafIntersection >
      LeafIntersectionIterator;
      typedef Dune::IntersectionIterator
      < const GridImp, AlbertaGridLeafIntersectionIterator, AlbertaGridLeafIntersection >
      LevelIntersectionIterator;

      typedef Dune::HierarchicIterator<const GridImp, AlbertaGridHierarchicIterator> HierarchicIterator;

      typedef IdType GlobalIdType;
      typedef IdType LocalIdType;

      template< int cd >
      struct Codim
      {
        // IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
        typedef Dune::Geometry<dim-cd, dimworld, const GridImp, AlbertaGridGlobalGeometry> Geometry;
        typedef Dune::Geometry<dim-cd, dim, const GridImp, AlbertaGridGeometry> LocalGeometry;

        typedef Dune::Entity< cd, dim, const GridImp, AlbertaGridEntity > Entity;

        typedef AlbertaGridEntityPointer< cd, const GridImp > EntityPointerImpl;
        typedef Dune::EntityPointer< const GridImp, EntityPointerImpl > EntityPointer;

        template <PartitionIteratorType pitype>
        struct Partition
        {
          typedef Dune::LevelIterator<cd,pitype,const GridImp,AlbertaGridLevelIterator> LevelIterator;
          typedef Dune::LeafIterator<cd,pitype,const GridImp,AlbertaGridLeafIterator> LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
      };

      template <PartitionIteratorType pitype>
      struct Partition
      {
        typedef Dune::GridView<DefaultLevelGridViewTraits<const GridImp,pitype> >
        LevelGridView;
        typedef Dune::GridView<DefaultLeafGridViewTraits<const GridImp,pitype> >
        LeafGridView;
      };

      typedef IndexSet< GridImp, LevelIndexSetImp, int > LevelIndexSet;
      typedef IndexSet< GridImp, LeafIndexSetImp, int > LeafIndexSet;
      typedef AlbertaGridHierarchicIndexSet< dim, dimworld > HierarchicIndexSet;
      typedef IdSet<GridImp,IdSetImp,IdType> GlobalIdSet;
      typedef IdSet<GridImp,IdSetImp,IdType> LocalIdSet;

      typedef Dune::CollectiveCommunication< AlbertaGrid<dim,dimworld> > CollectiveCommunication;
    };
  };

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTAGRID_GRIDFAMILTY_HH
