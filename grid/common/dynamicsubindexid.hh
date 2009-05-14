// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DYNAMICCODIMSUBINDEXID_HH
#define DUNE_DYNAMICCODIMSUBINDEXID_HH

#include <dune/grid/genericgeometry/misc.hh>

namespace Dune
{

  // DynamicSubIndex
  // ---------------

  template< class Grid, class IndexSet >
  class DynamicSubIndex
  {
    typedef DynamicSubIndex< Grid, IndexSet > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    static const unsigned int dimension = remove_const< Grid >::type::dimension;

    typedef typename Traits::template Codim< 0 >::Entity Element;

  public:
    typedef typename IndexSet::IndexType IndexType;

  private:
    struct Caller
    {
      virtual ~Caller ()
      {}

      virtual IndexType
      subIndex ( const IndexSet &indexSet, const Element &e, int i ) const = 0;
    };

    template< int codim >
    struct CallerImpl
      : public Caller
    {
      virtual IndexType
      subIndex ( const IndexSet &indexSet, const Element &e, int i ) const
      {
        typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
        const unsigned int tid = GenericGeometry::topologyId( e.type() );
        const int j = Numbering::template generic2dune< codim >( tid, i );
        return indexSet.template subIndex< codim >( e, j );
      }

      static void apply ( const Caller *(&caller)[ dimension+1 ] )
      {
        caller[ codim ] = new CallerImpl< codim >;
      }
    };

    // prohibit copying and assignment
    DynamicSubIndex ( const This & );
    This &operator= ( const This & );

  public:
    explicit DynamicSubIndex ( const IndexSet &indexSet )
      : indexSet_( indexSet )
    {
      GenericGeometry::ForLoop< CallerImpl, 0, dimension >::apply( caller_ );
    }

    ~DynamicSubIndex ()
    {
      for( unsigned int codim = 0; codim <= dimension; ++codim )
        delete caller_[ codim ];
    }

    IndexType operator() ( const Element &e, int i, unsigned int codim ) const
    {
      assert( codim <= dimension );
      return caller_[ codim ]->subIndex( indexSet_, e, i );
    }

  private:
    const IndexSet &indexSet_;
    const Caller *caller_[ dimension+1 ];
  };



  // DynamicSubId
  // ------------

  template< class Grid, class IdSet >
  class DynamicSubId
  {
    typedef DynamicSubId< Grid, IdSet > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    static const unsigned int dimension = remove_const< Grid >::type::dimension;

    typedef typename Traits::template Codim< 0 >::Entity Element;

  public:
    typedef typename IdSet::IdType IdType;

  private:
    struct Caller
    {
      virtual ~Caller ()
      {}

      virtual IdType
      subId ( const IdSet &idSet, const Element &e, int i ) const = 0;
    };

    template< int codim >
    struct CallerImpl
      : public Caller
    {
      virtual IdType
      subId ( const IdSet &idSet, const Element &e, int i ) const
      {
        typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
        const unsigned int tid = GenericGeometry::topologyId( e.type() );
        const int j = Numbering::template generic2dune< codim >( tid, i );
        return idSet.template subId< codim >( e, j );
      }

      static void apply ( const Caller *(&caller)[ dimension+1 ] )
      {
        caller[ codim ] = new CallerImpl< codim >;
      }
    };

    // prohibit copying and assignment
    DynamicSubId ( const This & );
    This &operator= ( const This & );

  public:
    explicit DynamicSubId ( const IdSet &idSet )
      : idSet_( idSet )
    {
      GenericGeometry::ForLoop< CallerImpl, 0, dimension >::apply( caller_ );
    }

    ~DynamicSubId ()
    {
      for( int codim = 0; codim <= dimension; ++codim )
        delete caller_[ codim ];
    }

    IdType operator() ( const Element &e, int i, unsigned int codim ) const
    {
      assert( codim <= dimension );
      return caller_[ codim ]->subId( idSet_, e, i );
    }

  private:
    const IdSet &idSet_;
    const Caller *caller_[ dimension+1 ];
  };

}

#endif // #ifndef DUNE_DYNAMICCODIMSUBINDEXID_HH
