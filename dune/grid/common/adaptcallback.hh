// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_ADAPTCALLBACK_HH
#define DUNE_GRID_COMMON_ADAPTCALLBACK_HH

/** \file
 *  \author Martin Nolte
 *  \brief  interfaces and wrappers needed for the callback adaptation provided
 *          by AlbertaGrid and dune-ALUGrid
 */

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Grid, class Impl >
  class AdaptDataHandle;



  // AdaptDataHandleInterface
  // ------------------------

  /** \brief Interface class for the Grid's adapt method where the
            parameter is a AdaptDataHandleInterface.
  */
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
    /** \brief call back for activity to take place on father
              and all descendants before the descendants are removed

       \param father  entity where all descendants are going to be removed
    */
    void preCoarsening ( const Entity &father )
    {
      asImp().preCoarsening( father );
    }

    /** \brief call back for activity to take place on newly created
              elements below the father element.

       \param father  entity where all descendants were newly created
    */
    void postRefinement ( const Entity &father )
    {
      asImp().postRefinement( father );
    }

    void restrictLocal( const Entity &father, const Entity& son, bool initialize )
    {
      asImp().restrictLocal( father, son, initialize );
    }

    void prolongLocal( const Entity &father, const Entity& son, bool initialize )
    {
      asImp().prolongLocal( father, son, initialize );
    }

  protected:
    const Impl &asImp () const { return static_cast< const Impl & >( *this ); }
    Impl &asImp () { return static_cast< Impl & >( *this ); }
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

    void preCoarsening ( const Entity &father );
    void postRefinement ( const Entity &father );
  };


  // CombinedAdaptProlongRestrict
  // ----------------------------

  //! class for combining 2 index sets together for adaptation process
  template <class A, class B >
  class CombinedAdaptProlongRestrict
  {
    //! space A and B
    A& _a;
    B& _b;
  public:
    //! constructor storing the two references
    CombinedAdaptProlongRestrict ( A& a, B& b ) : _a ( a ) , _b ( b )
    {}

    //! restrict data to father
    template <class Entity>
    void restrictLocal ( const Entity &father, const Entity &son, bool initialize )
    {
      _a.restrictLocal(father,son,initialize);
      _b.restrictLocal(father,son,initialize);
    }

    //! prolong data to children
    template <class Entity>
    void prolongLocal ( const Entity &father, const Entity &son, bool initialize )
    {
      _a.prolongLocal(father,son,initialize);
      _b.prolongLocal(father,son,initialize);
    }
  };

} // end namespace Dune

#endif
