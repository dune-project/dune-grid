// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_YASPGRIDIDSET_HH
#define DUNE_GRID_YASPGRIDIDSET_HH


namespace Dune {

  //========================================================================
  /*!
     \brief persistent, globally unique Ids

   */
  //========================================================================

  template<class GridImp>
  class YaspGlobalIdSet : public IdSet<GridImp,YaspGlobalIdSet<GridImp>,
                              typename std::remove_const<GridImp>::type::PersistentIndexType >
                          /*
                             We used the remove_const to extract the Type from the mutable class,
                             because the const class is not instantiated yet.
                           */
  {
    typedef YaspGlobalIdSet< GridImp > This;

  public:
    //! define the type used for persistent indices
    typedef typename std::remove_const<GridImp>::type::PersistentIndexType IdType;

    using IdSet<GridImp, This, IdType>::subId;

    //! Only default-constructible
    YaspGlobalIdSet()
    {}

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return e.impl().persistentIndex();
    }

    //! get id of subentity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    IdType subId (const typename std::remove_const<GridImp>::type::Traits::template Codim< 0 >::Entity &e,
                  int i, unsigned int codim ) const
    {
      return e.impl().subPersistentIndex(i,codim);
    }

  };

} // namespace Dune

#endif //  DUNE_GRID_YASPGRIDIDSET_HH
