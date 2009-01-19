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

    typedef Alberta::DofVectorPointer< int > DofVectorPointer;
    typedef Alberta::DofAccess< dim, 0 > DofAccess;

    typedef Alberta::FillFlags< dim > FillFlags;

    struct SetLocal;
    struct CalcMaxLevel;
    struct Interpolation;

  public:
    typedef Alberta::ElementInfo< dim > ElementInfo;
    typedef Alberta::MeshPointer< dim > MeshPointer;
    typedef Alberta::HierarchyDofNumbering< dim > DofNumbering;

  private:
    DofVectorPointer level_;
    DofAccess dofAccess_;

  public:
    int operator() ( const Alberta::Element *element ) const
    {
      const int *array = (int *)level_;
      return std::abs( array[ dofAccess_( element, 0 ) ] );
    }

    int operator() ( const ElementInfo &elementInfo ) const
    {
      return (*this)( elementInfo.el() );
    }

    bool isNew ( const Alberta::Element *element ) const
    {
      const int *array = (int *)level_;
      return (array[ dofAccess_( element, 0 ) ] < 0);
    }

    bool isNew ( const ElementInfo &elementInfo ) const
    {
      return isNew( elementInfo.el() );
    }

    int maxLevel () const
    {
      const int maxLevel = Alberta::maxAbs( level_ );
#ifndef NDEBUG
      CalcMaxLevel calcMaxLevel;
      mesh().leafTraverse( calcMaxLevel, FillFlags::nothing );
      assert( maxLevel == calcMaxLevel.maxLevel() );
#endif
      return maxLevel;
    }

    MeshPointer mesh () const
    {
      return MeshPointer( level_.dofSpace()->mesh );
    }

    void markAllOld ()
    {
      Alberta::abs( level_ );
    }

    void create ( const DofNumbering &dofNumbering )
    {
      const Alberta::DofSpace *const dofSpace = dofNumbering.dofSpace( 0 );
      dofAccess_ = DofAccess( dofSpace );

      level_.create( dofSpace, "Element level" );
      assert( !(!level_) );
      level_.template setupInterpolation< Interpolation >();

      SetLocal setLocal( level_ );
      mesh().hierarchicTraverse( setLocal, FillFlags::nothing );
    }

    void release ()
    {
      level_.release();
      dofAccess_ = DofAccess();
    }
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
      int *const array = (int *)level_;
      array[ dofAccess_( elementInfo, 0 ) ] = elementInfo.level();
    }
  };



  // AlbertaGridLevelProvider::CalcMaxLevel
  // --------------------------------------

  template< int dim >
  class AlbertaGridLevelProvider< dim >::CalcMaxLevel
  {
    int maxLevel_;

  public:
    CalcMaxLevel ()
      : maxLevel_( 0 )
    {}

    void operator() ( const Alberta::ElementInfo< dim > &elementInfo )
    {
      maxLevel_ = std::max( maxLevel_, elementInfo.level() );
    }

    int maxLevel () const
    {
      return maxLevel_;
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
      int *array = (int *)dofVector;

      for( int i = 0; i < patch.count(); ++i )
      {
        const Alberta::Element *const father = patch[ i ];
        const int fatherLevel = std::abs( array[ dofAccess( father, 0 ) ] );
        for( int i = 0; i < 2; ++i )
        {
          const Alberta::Element *child = father->child[ i ];
          array[ dofAccess( child, 0 ) ] = -(fatherLevel+1);
        }
      }
    }
  };

}

#endif // #if HAVE_ALBERTA

#endif
