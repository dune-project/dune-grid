// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_LEVEL_HH
#define DUNE_ALBERTA_LEVEL_HH

#include <cassert>
#include <cstdlib>

#include <dune/grid/albertagrid/meshpointer.hh>
#include <dune/grid/albertagrid/dofadmin.hh>
#include <dune/grid/albertagrid/dofvector.hh>

#if HAVE_ALBERTA

namespace Dune
{

  // AlbertaGridLevelProvider
  // ------------------------

  template< int dim >
  class AlbertaGridLevelProvider
  {
    typedef AlbertaGridLevelProvider< dim > This;

    typedef unsigned char Level;

    typedef Alberta::DofVectorPointer< Level > DofVectorPointer;
    typedef Alberta::DofAccess< dim, 0 > DofAccess;

    typedef Alberta::FillFlags< dim > FillFlags;

    static const Level isNewFlag = (1 << 7);
    static const Level levelMask = (1 << 7) - 1;

    class SetLocal;
    class CalcMaxLevel;

    template< Level flags >
    struct ClearFlags;

    struct Interpolation;

  public:
    typedef Alberta::ElementInfo< dim > ElementInfo;
    typedef Alberta::MeshPointer< dim > MeshPointer;
    typedef Alberta::HierarchyDofNumbering< dim > DofNumbering;

    Level operator() ( const Alberta::Element *element ) const
    {
      const Level *array = (Level *)level_;
      return array[ dofAccess_( element, 0 ) ] & levelMask;
    }

    Level operator() ( const ElementInfo &elementInfo ) const
    {
      return (*this)( elementInfo.el() );
    }

    bool isNew ( const Alberta::Element *element ) const
    {
      const Level *array = (Level *)level_;
      return ((array[ dofAccess_( element, 0 ) ] & isNewFlag) != 0);
    }

    bool isNew ( const ElementInfo &elementInfo ) const
    {
      return isNew( elementInfo.el() );
    }

    Level maxLevel () const
    {
      CalcMaxLevel calcFromCache;
      level_.forEach( calcFromCache );
#ifndef NDEBUG
      CalcMaxLevel calcFromGrid;
      mesh().leafTraverse( calcFromGrid, FillFlags::nothing );
      assert( calcFromCache.maxLevel() == calcFromGrid.maxLevel() );
#endif
      return calcFromCache.maxLevel();;
    }

    MeshPointer mesh () const
    {
      return MeshPointer( level_.dofSpace()->mesh );
    }

    void markAllOld ()
    {
      ClearFlags< isNewFlag > clearIsNew;
      level_.forEach( clearIsNew );
    }

    void create ( const DofNumbering &dofNumbering )
    {
      const Alberta::DofSpace *const dofSpace = dofNumbering.dofSpace( 0 );
      dofAccess_ = DofAccess( dofSpace );

      level_.create( dofSpace, "Element level" );
      assert( level_ );
      level_.template setupInterpolation< Interpolation >();

      SetLocal setLocal( level_ );
      mesh().hierarchicTraverse( setLocal, FillFlags::nothing );
    }

    void release ()
    {
      level_.release();
      dofAccess_ = DofAccess();
    }

  private:
    DofVectorPointer level_;
    DofAccess dofAccess_;
  };



  // AlbertaGridLevelProvider::SetLocal
  // ----------------------------------

  template< int dim >
  class AlbertaGridLevelProvider< dim >::SetLocal
  {
    DofVectorPointer level_;
    DofAccess dofAccess_;

  public:
    explicit SetLocal ( const DofVectorPointer &level )
      : level_( level ),
        dofAccess_( level.dofSpace() )
    {}

    void operator() ( const Alberta::ElementInfo< dim > &elementInfo ) const
    {
      Level *const array = (Level *)level_;
      array[ dofAccess_( elementInfo, 0 ) ] = elementInfo.level();
    }
  };



  // AlbertaGridLevelProvider::CalcMaxLevel
  // --------------------------------------

  template< int dim >
  class AlbertaGridLevelProvider< dim >::CalcMaxLevel
  {
    Level maxLevel_;

  public:
    CalcMaxLevel ()
      : maxLevel_( 0 )
    {}

    void operator() ( const Level &dof )
    {
      maxLevel_ = std::max( maxLevel_, Level( dof & levelMask ) );
    }

    void operator() ( const Alberta::ElementInfo< dim > &elementInfo )
    {
      maxLevel_ = std::max( maxLevel_, Level( elementInfo.level() ) );
    }

    Level maxLevel () const
    {
      return maxLevel_;
    }
  };



  // AlbertaGridLevelProvider::ClearFlags
  // ------------------------------------

  template< int dim >
  template< typename AlbertaGridLevelProvider< dim >::Level flags >
  struct AlbertaGridLevelProvider< dim >::ClearFlags
  {
    void operator() ( Level &dof ) const
    {
      dof &= ~flags;
    }
  };



  // AlbertaGridLevelProvider::Interpolation
  // ---------------------------------------

  template< int dim >
  struct AlbertaGridLevelProvider< dim >::Interpolation
  {
    static const int dimension = dim;

    typedef Alberta::Patch< dimension > Patch;

    static void interpolateVector ( const DofVectorPointer &dofVector,
                                    const Patch &patch )
    {
      const DofAccess dofAccess( dofVector.dofSpace() );
      Level *array = (Level *)dofVector;

      for( int i = 0; i < patch.count(); ++i )
      {
        const Alberta::Element *const father = patch[ i ];
        assert( (array[ dofAccess( father, 0 ) ] & levelMask) < levelMask );
        const Level childLevel = (array[ dofAccess( father, 0 ) ] + 1) | isNewFlag;
        for( int i = 0; i < 2; ++i )
        {
          const Alberta::Element *child = father->child[ i ];
          array[ dofAccess( child, 0 ) ] = childLevel;
        }
      }
    }
  };

}

#endif // #if HAVE_ALBERTA

#endif
