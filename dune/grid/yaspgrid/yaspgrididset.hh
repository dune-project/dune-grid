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

  template< class GridImp >
  class YaspGlobalIdSet
    : public IdSet< GridImp, YaspGlobalIdSet< GridImp >, typename remove_const< GridImp >::type::PersistentIndexType >
  {
    typedef YaspGlobalIdSet< GridImp > This;
    typedef IdSet< GridImp, This, typename remove_const< GridImp >::type::PersistentIndexType > Base;

  public:
    //! define the type used for persisitent indices
    typedef typename Base::IdType IdType;

    using Base::subId;

    //! get id of an entity
    template< int cd >
    IdType id ( const typename remove_const< GridImp >::type::Traits::template Codim< cd >::Entity &e ) const
    {
      return GridImp::getRealImplementation( e ).persistentIndex();
    }

    //! get id of subentity
    template< int cc >
    IdType subId ( const typename remove_const< GridImp >::type::Traits::template Codim< cc >::Entity &e,
                   int i, unsigned int codim ) const
    {
      return GridImp::getRealImplementation( e ).subPersistentIndex( i, codim );
    }
  };

} // namespace Dune

#endif //  DUNE_GRID_YASPGRIDIDSET_HH
