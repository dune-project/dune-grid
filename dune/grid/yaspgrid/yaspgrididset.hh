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
                              typename remove_const<GridImp>::type::PersistentIndexType >
                          /*
                             We used the remove_const to extract the Type from the mutable class,
                             because the const class is not instantiated yet.
                           */
  {
    typedef YaspGlobalIdSet< GridImp > This;

  public:
    //! define the type used for persisitent indices
    typedef typename remove_const<GridImp>::type::PersistentIndexType IdType;

    using IdSet<GridImp, This, IdType>::subId;

    //! constructor stores reference to a grid
    explicit YaspGlobalIdSet ( const GridImp &g )
      : grid( g )
    {}

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return grid.getRealImplementation(e).persistentIndex();
    }

    //! get id of subentity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    IdType subId (const typename remove_const<GridImp>::type::Traits::template Codim< 0 >::Entity &e,
                  int i, unsigned int codim ) const
    {
      return grid.getRealImplementation(e).subPersistentIndex(i,codim);
    }

  private:
    const GridImp& grid;
  };

} // namespace Dune

#endif //  DUNE_GRID_YASPGRIDIDSET_HH
