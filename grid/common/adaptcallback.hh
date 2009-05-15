// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ADAPTCALLBACK_HH
#define DUNE_ADAPTCALLBACK_HH

/** \file
 *  \author Martin Nolte
 *  \brief  interfaces and wrappers needed for the callback adaptation provided
 *          by AlbertaGrid and ALUGrid
 */

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Grid, class Impl >
  class AdaptDataHandle;



  // AdaptDataHandleInterface
  // ------------------------

  template< class Grid, class Impl >
  class AdaptDataHandleInterface
  {
    typedef AdaptDataHandleInterface< Grid, Impl > This;

    friend class AdaptDataHandle< Grid, Impl >;

  public:
    typedef typename Grid::template Codim< 0 >::Entity Entity;

  private:
    AdaptDataHandleInterface ()
    {}

    AdaptDataHandleInterface ( const This & );
    This &operator= ( const This & );

  public:
    void preAdapt ( const unsigned int estimateAdditionalElements )
    {
      asImp().preAdapt( estimateAdditionalElements );
    }

    void postAdapt ()
    {
      asImp().postAdapt();
    }

    void preCoarsening ( const Entity &father ) const
    {
      asImp().preCoarsening( father );
    }

    void postRefinement ( const Entity &father ) const
    {
      asImp().postRefinement( father );
    }

  protected:
    const Impl &asImp () const
    {
      return static_cast< const Impl & >( *this );
    }

    Impl &asImp ()
    {
      return static_cast< Impl & >( *this );
    }
  };



  // AdaptDataHandle
  // ---------------

  template< class Grid, class Impl >
  class AdaptDataHandle
    : public AdaptDataHandleInterface< Grid, Impl >
  {
    typedef AdaptDataHandle< Grid, Impl > This;
    typedef AdaptDataHandleInterface< Grid, Impl > Base;

  public:
    typedef typename Base::Entity Entity;

  protected:
    AdaptDataHandle ()
    {}

  private:
    AdaptDataHandle ( const This & );
    This &operator= ( const This & );

    void preAdapt ( const unsigned int estimateAdditionalElements );
    void postAdapt ();
    void preCoarsening ( const Entity &father ) const;
    void postRefinement ( const Entity &father ) const;
  };



  // RestrictProlongWrapper
  // ----------------------

  template< class Grid, class DofManager, class RestrictProlongOperator >
  class RestrictProlongWrapper
    : public AdaptDataHandle
      < Grid, RestrictProlongWrapper< Grid, DofManager, RestrictProlongOperator > >
  {
    typedef RestrictProlongWrapper< Grid, DofManager, RestrictProlongOperator > This;
    typedef AdaptDataHandle< Grid, This > Base;

    DofManager &dofManager_;
    RestrictProlongOperator &rpOp_;

  public:
    typedef typename Base::Entity Entity;

    RestrictProlongWrapper ( DofManager &dofManager, RestrictProlongOperator &rpOp )
      : dofManager_( dofManager ),
        rpOp_( rpOp )
    {}

    void preAdapt ( const unsigned int estimatedAdditionalElements )
    {
      dofManager_.reserveMemory( estimatedAdditionalElements );
    }

    void postAdapt ()
    {
      dofManager_.compress();
    }

    void preCoarsening ( const Entity &father ) const
    {
      typedef typename Entity::HierarchicIterator HIterator;

      bool initialize = true;
      const int childLevel = father.level()+1;
      const HIterator end = father.hend( childLevel );
      for( HIterator it = father.hbegin( childLevel ); it != end; ++it )
      {
        restrictLocal( father, *it, initialize );
        initialize = false;
      }
    }

    void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      dofManager_.indexSetRestrictProlong().restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      rpOp_.restrictLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
    }

    void postRefinement ( const Entity &father ) const
    {
      typedef typename Entity::HierarchicIterator HIterator;

      bool initialize = true;
      const int childLevel = father.level()+1;
      const HIterator end = father.hend( childLevel );
      for( HIterator it = father.hbegin( childLevel ); it != end; ++it )
      {
        prolongLocal( father, *it, initialize );
        initialize = false;
      }
    }

    void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      dofManager_.indexSetRestrictProlong().prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
      rpOp_.prolongLocal( const_cast< Entity & >( father ), const_cast< Entity & >( son ), initialize );
    }
  };



  // CombinedAdaptProlongRestrict
  // ----------------------------

  //! class for combining 2 index sets together for adaptation process
  template <class A, class B >
  class CombinedAdaptProlongRestrict
  {
    //! space A and B
    const A & _a;
    const B & _b;
  public:
    //! constructor storing the two references
    CombinedAdaptProlongRestrict ( const A & a, const B & b ) : _a ( a ) , _b ( b )
    {}

    //! restrict data to father
    template <class EntityType>
    void restrictLocal ( EntityType &father, EntityType &son, bool initialize ) const
    {
      _a.restrictLocal(father,son,initialize);
      _b.restrictLocal(father,son,initialize);
    }

    //! prolong data to children
    template <class EntityType>
    void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
    {
      _a.prolongLocal(father,son,initialize);
      _b.prolongLocal(father,son,initialize);
    }
  };

} // end namespace Dune

#endif
