// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ENTITYSEED_HH
#define DUNE_ALBERTA_ENTITYSEED_HH

#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/meshpointer.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int codim, class Grid >
  class AlbertaGridEntitySeed;



  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



#if HAVE_ALBERTA

  // AlbertaGridEntitySeed (for higher codimension)
  // ----------------------------------------------

  template< int codim, int dim, int dimworld >
  class AlbertaGridEntitySeed< codim, const AlbertaGrid< dim, dimworld > >
  {
  public:
    typedef AlbertaGrid< dim, dimworld > Grid;

    static const int codimension = codim;
    static const int dimension = dim;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = dimworld;

    typedef Alberta::MeshPointer< dimension > MeshPointer;
    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef typename ElementInfo::Seed Seed;

    typedef typename Grid::template Codim< codimension >::Entity Entity;

    AlbertaGridEntitySeed ( )
    {}

    AlbertaGridEntitySeed ( const ElementInfo &elementInfo, int subEntity )
      : seed_( elementInfo.seed() ),
        subEntity_( subEntity )
    {}

    bool isValid () const
    {
      return seed_.isValid();
    }

    ElementInfo elementInfo ( const MeshPointer &mesh ) const { return ElementInfo( mesh, seed_ ); }
    int subEntity () const { return subEntity_; }

  private:
    Seed seed_;
    int subEntity_;
  };



  // AlbertaGridEntitySeed (for codimension 0)
  // -----------------------------------------

  template< int dim, int dimworld >
  class AlbertaGridEntitySeed< 0, const AlbertaGrid< dim, dimworld > >
  {
  public:
    typedef AlbertaGrid< dim, dimworld > Grid;

    static const int codimension = 0;
    static const int dimension = dim;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = dimworld;

    typedef Alberta::MeshPointer< dimension > MeshPointer;
    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef typename ElementInfo::Seed Seed;

    typedef typename Grid::template Codim< codimension >::Entity Entity;

    AlbertaGridEntitySeed ( )
    {}

    explicit AlbertaGridEntitySeed ( const ElementInfo &elementInfo )
      : seed_( elementInfo.seed() )
    {}

    bool isValid () const
    {
      return seed_.isValid();
    }

    ElementInfo elementInfo ( const MeshPointer &mesh ) const { return ElementInfo( mesh, seed_ ); }
    int subEntity () const { return 0; }

  private:
    Seed seed_;
  };

#endif // #if HAVE_ALBERTA

} // end namespace Dune

#endif // #ifndef DUNE_ALBERTA_ENTITYSEED_HH
