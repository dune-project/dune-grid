// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_LBDATAHANDLE_HH
#define DUNE_ALUGRID_LBDATAHANDLE_HH

#include <dune/grid/alugrid/3d/datahandle.hh>

#include <dune/grid/common/datahandleif.hh>

namespace Dune
{

  template< class Grid, class DataHandleImpl, class Data >
  class ALUGridLoadBalanceDataHandle
  {
    typedef typename Grid :: Traits :: HierarchicIterator HierarchicIterator;

  public:
    typedef typename Grid :: ObjectStreamType ObjectStream;

    typedef CommDataHandleIF< DataHandleImpl, Data > DataHandle;

    static const int dimension = Grid :: dimension;

    template< int codim >
    struct Codim
    {
      typedef typename Grid :: Traits :: template Codim< codim > :: Entity Entity;
      typedef typename Grid :: Traits :: template Codim< codim > :: EntityPointer
      EntityPointer;
    };

    typedef typename Codim< 0 > :: Entity Element;

  private:
    const Grid &grid_;
    DataHandle &dataHandle_;

  public:
    ALUGridLoadBalanceDataHandle ( const Grid &grid, DataHandle &dataHandle )
      : grid_( grid ),
        dataHandle_( dataHandle )
    {}

    void inlineData ( ObjectStream &stream, const Element &element ) const
    {
      inlineElementData( stream, element );

      const int maxLevel = grid_.maxLevel();
      const HierarchicIterator end = element.hend( maxLevel );
      for( HierarchicIterator it = element.hbegin( maxLevel ); it != end; ++it )
        inlineElementData( stream, *it );
    }

    void xtractData ( ObjectStream &stream, const Element &element, size_t newElements )
    {
      xtractElementData( stream, element );

      const int maxLevel = grid_.maxLevel();
      const HierarchicIterator end = element.hend( maxLevel );
      for( HierarchicIterator it = element.hbegin( maxLevel ); it != end; ++it )
        xtractElementData( stream, *it );
    }

    void compress ()
    {}

  private:
    void inlineElementData ( ObjectStream &stream, const Element &element ) const
    {
      // call element data direct without creating entity pointer
      if( dataHandle_.contains( dimension, 0 ) )
      {
        inlineEntityData( stream, element );
      }

      // now call all higher codims
      inlineCodimData< 1 >( stream, element );
      inlineCodimData< 2 >( stream, element );
      inlineCodimData< 3 >( stream, element );
    }

    void xtractElementData ( ObjectStream &stream, const Element &element )
    {
      // call element data direct without creating entity pointer
      if( dataHandle_.contains( dimension, 0 ) )
      {
        xtractEntityData( stream, element );
      }

      // now call all higher codims
      xtractCodimData< 1 >( stream, element );
      xtractCodimData< 2 >( stream, element );
      xtractCodimData< 3 >( stream, element );
    }

    template< int codim >
    void inlineCodimData ( ObjectStream &stream, const Element &element ) const
    {
      typedef typename Codim< codim > :: EntityPointer EntityPointer;

      if( dataHandle_.contains( dimension, codim ) )
      {
        const int numSubEntities = element.template count< codim >();
        for( int i = 0; i < numSubEntities; ++i )
        {
          const EntityPointer pEntity = element.template entity< codim >( i );
          inlineEntityData< codim >( stream, *pEntity );
        }
      }
    }

    template< int codim >
    void xtractCodimData ( ObjectStream &stream, const Element &element )
    {
      typedef typename Codim< codim > :: EntityPointer EntityPointer;

      if( dataHandle_.contains( dimension, codim ) )
      {
        const int numSubEntities = element.template count< codim >();
        for( int i = 0; i < numSubEntities; ++i )
        {
          const EntityPointer pEntity = element.template entity< codim >( i );
          xtractEntityData< codim >( stream, *pEntity );
        }
      }
    }

    template< int codim >
    void inlineEntityData ( ObjectStream &stream,
                            const typename Codim< codim > :: Entity &entity ) const
    {
      const size_t size = dataHandle_.size( entity );
      stream.write( size );
      dataHandle_.gather( stream, entity );
    }

    template< int codim >
    void xtractEntityData ( ObjectStream &stream,
                            const typename Codim< codim > :: Entity &entity )
    {
      size_t size = 0;
      stream.read( size );
      dataHandle_.scatter( stream, entity, size );
    }
  };

}

#endif
