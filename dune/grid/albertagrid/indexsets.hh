// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRIDINDEXSETS_HH
#define DUNE_ALBERTAGRIDINDEXSETS_HH

#include <array>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/stdstreams.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/albertagrid/indexstack.hh>
#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/dofadmin.hh>
#include <dune/grid/albertagrid/dofvector.hh>
#include <dune/grid/albertagrid/elementinfo.hh>
#include <dune/grid/albertagrid/gridfamily.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {
    typedef Dune::IndexStack< int, 100000 > IndexStack;
  }



  // AlbertaGridHierarchicIndexSet
  // -----------------------------

  template< int dim, int dimworld >
  class AlbertaGridHierarchicIndexSet
    : public IndexSet< AlbertaGridFamily< dim, dimworld >, AlbertaGridHierarchicIndexSet< dim,dimworld >, int, std::array< GeometryType, 1 > >
  {
    typedef AlbertaGridHierarchicIndexSet< dim, dimworld > This;
    typedef IndexSet< AlbertaGridFamily< dim, dimworld >, AlbertaGridHierarchicIndexSet< dim,dimworld >, int, std::array< GeometryType, 1 > > Base;

    friend class AlbertaGrid< dim, dimworld >;

  public:
    typedef AlbertaGrid< dim, dimworld > Grid;
    typedef AlbertaGridFamily< dim, dimworld > GridFamily;

    typedef typename Base::IndexType IndexType;

    typedef typename Base::Types Types;

    static const int dimension = GridFamily::dimension;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef Alberta::HierarchyDofNumbering< dimension > DofNumbering;

  private:
    typedef typename GridFamily::Traits Traits;

    typedef Alberta::DofVectorPointer< IndexType > IndexVectorPointer;

    class InitEntityNumber;

    template< int codim >
    struct CreateEntityNumbers;

    template< int codim >
    struct RefineNumbering;

    template< int codim >
    struct CoarsenNumbering;

    explicit AlbertaGridHierarchicIndexSet ( const DofNumbering &dofNumbering );

    static Alberta::IndexStack *currentIndexStack;

  public:
    typedef Alberta::IndexStack IndexStack;

    //! return true if entity is contained in set
    template< class Entity >
    bool contains ( const Entity & ) const
    {
      return true;
    }

    using Base::index;
    using Base::subIndex;

    //! return hierarchic index of given entity
    template< int cc >
    IndexType index ( const typename Traits::template Codim< cc >::Entity &entity ) const
    {
      typedef AlbertaGridEntity< cc, dim, const Grid > EntityImp;
      const EntityImp &entityImp = entity.impl();
      return subIndex( entityImp.elementInfo(), entityImp.subEntity(), cc );
    }

    //! return subIndex of given enitiy's sub entity
    template< int cc >
    IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &entity, int i, unsigned int codim ) const
    {
      typedef AlbertaGridEntity< cc, dim, const Grid > EntityImp;
      const EntityImp &entityImp = entity.impl();

      int k = i;
      if( cc > 0 )
      {
        auto refElement = ReferenceElements< Alberta::Real, dimension >::simplex();
        k = refElement.subEntity( entityImp.subEntity(), cc, i, codim );
      }

      const int j = entityImp.grid().generic2alberta( codim, k );
      return subIndex( entityImp.elementInfo(), j, codim );
    }

    //! return size of set for given GeometryType
    std::size_t size ( const GeometryType &type ) const
    {
      return (type.isSimplex() ? size( dimension - type.dim() ) : 0);
    }

    //! return size of set
    std::size_t size ( int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return indexStack_[ codim ].size();
    }

    Types types ( int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return {{ GeometryTypes::simplex( dimension - codim ) }};
    }

    //! return geometry types this set has indices for
    const std::vector< GeometryType > &geomTypes( int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return geomTypes_[ codim ];
    }

    IndexType subIndex ( const ElementInfo &elementInfo, int i, unsigned int codim ) const
    {
      assert( !elementInfo == false );
      return subIndex( elementInfo.element(), i, codim );
    }

    /** \brief obtain hierarchic subindex
     *
     *  \param[in]  element  pointer to ALBERTA element
     *  \param[in]  i        number of the subelement (in ALBERTA numbering)
     *  \param[in]  codim    codimension
     */
    IndexType subIndex ( const Alberta::Element *element, int i, unsigned int codim ) const
    {
      IndexType *array = (IndexType *)entityNumbers_[ codim ];
      const IndexType subIndex = array[ dofNumbering_( element, codim, i ) ];
      assert( (subIndex >= 0) && (subIndex < IndexType(size( codim ))) );
      return subIndex;
    }

    void preAdapt ()
    {
      // set global pointer to index stack
      if( !IndexVectorPointer::supportsAdaptationData )
      {
        assert( currentIndexStack == nullptr );
        currentIndexStack = indexStack_;
      }
    }

    void postAdapt ()
    {
      // remove global pointer to index stack
      if( !IndexVectorPointer::supportsAdaptationData )
        currentIndexStack = nullptr;
    }

    void create ();
    void read ( const std::string &filename );
    bool write ( const std::string &filename ) const;

    void release ()
    {
      for( int i = 0; i <= dimension; ++i )
        entityNumbers_[ i ].release();
    }

  private:
    template< int codim >
    static IndexStack &getIndexStack ( const IndexVectorPointer &dofVector )
    {
      IndexStack *indexStack;
      if( IndexVectorPointer::supportsAdaptationData )
        indexStack = dofVector.template getAdaptationData< IndexStack >();
      else
        indexStack = &currentIndexStack[ codim ];
      assert( indexStack != 0 );
      return *indexStack;
    }

    // access to the dof vectors
    const DofNumbering &dofNumbering_;

    // index stacks providing new numbers during adaptation
    IndexStack indexStack_[ dimension+1 ];

    // dof vectors storing the (persistent) numbering
    IndexVectorPointer entityNumbers_[ dimension+1 ];

    // all geometry types contained in the grid
    std::vector< GeometryType > geomTypes_[ dimension+1 ];
  };



  // AlbertaGridHierarchicIndexSet::InitEntityNumber
  // -----------------------------------------------

  template< int dim, int dimworld >
  class AlbertaGridHierarchicIndexSet< dim, dimworld >::InitEntityNumber
  {
    IndexStack &indexStack_;

  public:
    InitEntityNumber ( IndexStack &indexStack )
      : indexStack_( indexStack )
    {}

    void operator() ( int &dof )
    {
      dof = indexStack_.getIndex();
    }
  };



  // AlbertaGridHierarchicIndexSet::CreateEntityNumbers
  // --------------------------------------------------

  template< int dim, int dimworld >
  template< int codim >
  struct AlbertaGridHierarchicIndexSet< dim, dimworld >::CreateEntityNumbers
  {
    static void setup ( AlbertaGridHierarchicIndexSet< dim, dimworld > &indexSet );

    static void apply ( const Alberta::HierarchyDofNumbering< dimension > &dofNumbering,
                        AlbertaGridHierarchicIndexSet< dim, dimworld > &indexSet );

    static void apply ( const std::string &filename,
                        const Alberta::MeshPointer< dimension > &mesh,
                        AlbertaGridHierarchicIndexSet< dim, dimworld > &indexSet );
  };



  // AlbertaGridHierarchicIndexSet::RefineNumbering
  // ----------------------------------------------

  template< int dim, int dimworld >
  template< int codim >
  struct AlbertaGridHierarchicIndexSet< dim, dimworld >::RefineNumbering
  {
    static const int dimension = dim;
    static const int codimension = codim;

  private:
    typedef Alberta::DofAccess< dimension, codimension > DofAccess;

    explicit RefineNumbering ( const IndexVectorPointer &dofVector )
      : indexStack_( getIndexStack< codimension >( dofVector ) ),
        dofVector_( dofVector ),
        dofAccess_( dofVector.dofSpace() )
    {}

  public:
    void operator() ( const Alberta::Element *child, int subEntity );

    typedef Alberta::Patch< dimension > Patch;
    static void interpolateVector ( const IndexVectorPointer &dofVector,
                                    const Patch &patch );

  private:
    IndexStack &indexStack_;
    IndexVectorPointer dofVector_;
    DofAccess dofAccess_;
  };



  // AlbertaGridHierarchicIndexSet::CoarsenNumbering
  // -----------------------------------------------

  template< int dim, int dimworld >
  template< int codim >
  struct AlbertaGridHierarchicIndexSet< dim, dimworld >::CoarsenNumbering
  {
    static const int dimension = dim;
    static const int codimension = codim;

  private:
    typedef Alberta::DofAccess< dimension, codimension > DofAccess;

    explicit CoarsenNumbering ( const IndexVectorPointer &dofVector )
      : indexStack_( getIndexStack< codimension >( dofVector ) ),
        dofVector_( dofVector ),
        dofAccess_( dofVector.dofSpace() )
    {}

  public:
    void operator() ( const Alberta::Element *child, int subEntity );

    typedef Alberta::Patch< dimension > Patch;
    static void restrictVector ( const IndexVectorPointer &dofVector,
                                 const Patch &patch );
  private:
    IndexStack &indexStack_;
    IndexVectorPointer dofVector_;
    DofAccess dofAccess_;
  };



  // AlbertaGridIndexSet
  // -------------------

  template< int dim, int dimworld >
  class AlbertaGridIndexSet
    : public IndexSet< AlbertaGrid< dim, dimworld >, AlbertaGridIndexSet< dim, dimworld >, int, std::array< GeometryType, 1 > >
  {
    typedef AlbertaGridIndexSet< dim, dimworld > This;
    typedef IndexSet< AlbertaGrid< dim, dimworld >, AlbertaGridIndexSet< dim, dimworld >, int, std::array< GeometryType, 1 > > Base;

  public:
    typedef AlbertaGrid< dim, dimworld > Grid;

    typedef typename Base::IndexType IndexType;

    typedef typename Base::Types Types;

    static const int dimension = Grid::dimension;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef Alberta::HierarchyDofNumbering< dimension > DofNumbering;

  private:
    typedef typename Grid::Traits Traits;

    template< int codim >
    struct Insert;

  public:
    explicit AlbertaGridIndexSet ( const DofNumbering &dofNumbering )
      : dofNumbering_( dofNumbering )
    {
      for( int codim = 0; codim <= dimension; ++codim )
      {
        indices_[ codim ] = 0;
        geomTypes_[ codim ].push_back( GeometryTypes::simplex( dimension - codim ) );
      }
    }

    ~AlbertaGridIndexSet ()
    {
      for( int codim = 0; codim <= dimension; ++codim )
        delete[] indices_[ codim ];
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      const int codim = Entity::codimension;

      const AlbertaGridEntity< codim, dim, const Grid > &entityImp
        = entity.impl();
      const Alberta::Element *element = entityImp.elementInfo().el();

      const IndexType *const array = indices_[ codim ];
      const IndexType subIndex = array[ dofNumbering_( element, codim, entityImp.subEntity() ) ];

      return (subIndex >= 0);
    }

    using Base::index;
    using Base::subIndex;

    //! return hierarchic index of given entity
    template< int cc >
    IndexType index ( const typename Traits::template Codim< cc >::Entity &entity ) const
    {
      typedef AlbertaGridEntity< cc, dim, const Grid > EntityImp;
      const EntityImp &entityImp = entity.impl();
      return subIndex( entityImp.elementInfo(), entityImp.subEntity(), cc );
    }

    //! return subIndex of given enitiy's sub entity
    template< int cc >
    IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &entity, int i, unsigned int codim ) const
    {
      typedef AlbertaGridEntity< cc, dim, const Grid > EntityImp;
      const EntityImp &entityImp = entity.impl();

      int k = i;
      if( cc > 0 )
      {
        auto refElement = ReferenceElements< Alberta::Real, dimension >::simplex();
        k = refElement.subEntity( entityImp.subEntity(), cc, i, codim );
      }

      const int j = entityImp.grid().generic2alberta( codim, k );
      return subIndex( entityImp.elementInfo(), j, codim );
    }

    std::size_t size ( const GeometryType &type ) const
    {
      return (type.isSimplex() ? size( dimension - type.dim() ) : 0);
    }

    std::size_t size ( int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return size_[ codim ];
    }

    Types types ( int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return {{ GeometryTypes::simplex( dimension - codim ) }};
    }

    const std::vector< GeometryType > &geomTypes( int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return geomTypes_[ codim ];
    }

    template< class Iterator >
    void update ( const Iterator &begin, const Iterator &end )
    {
      for( int codim = 0; codim <= dimension; ++codim )
      {
        delete[] indices_[ codim ];

        const unsigned int dofSize = dofNumbering_.size( codim );
        indices_[ codim ] = new IndexType[ dofSize ];
        for( unsigned int i = 0; i < dofSize; ++i )
          indices_[ codim ][ i ] = -1;

        size_[ codim ] = 0;
      }

      for( Iterator it = begin; it != end; ++it )
      {
        const AlbertaGridEntity< 0, dim, const Grid > &entityImp
          = it->impl();
        const Alberta::Element *element = entityImp.elementInfo().el();
        Hybrid::forEach( std::make_index_sequence< dimension+1 >{},
          [ & ]( auto i ){ Insert< i >::apply( element, *this ); } );
      }
    }

  private:
    IndexType subIndex ( const ElementInfo &elementInfo, int i, unsigned int codim ) const
    {
      assert( !elementInfo == false );
      return subIndex( elementInfo.element(), i, codim );
    }

    /** \brief obtain subindex
     *
     *  \param[in]  element  pointer to ALBERTA element
     *  \param[in]  i        number of the subelement (in ALBERTA numbering)
     *  \param[in]  codim    codimension
     */
    IndexType subIndex ( const Alberta::Element *element, int i, unsigned int codim ) const
    {
      const IndexType *const array = indices_[ codim ];
      const IndexType subIndex = array[ dofNumbering_( element, codim, i ) ];
      assert( (subIndex >= 0) && (static_cast<unsigned int>(subIndex) < size( codim )) );
      return subIndex;
    }

    // access to the dof vectors
    const DofNumbering &dofNumbering_;

    // an array of indices for each codimension
    IndexType *indices_[ dimension+1 ];

    // the size of each codimension
    IndexType size_[ dimension+1 ];

    // all geometry types contained in the grid
    std::vector< GeometryType > geomTypes_[ dimension+1 ];
  };



  // AlbertaGridIndexSet::Insert
  // ---------------------------

  template< int dim, int dimworld >
  template< int codim >
  struct AlbertaGridIndexSet< dim, dimworld >::Insert
  {
    static void apply ( const Alberta::Element *const element,
                        AlbertaGridIndexSet< dim, dimworld > &indexSet )
    {
      int *const array = indexSet.indices_[ codim ];
      IndexType &size = indexSet.size_[ codim ];

      for( int i = 0; i < Alberta::NumSubEntities< dim, codim >::value; ++i )
      {
        int &index = array[ indexSet.dofNumbering_( element, codim, i ) ];
        if( index < 0 )
          index = size++;
      }
    }
  };



  // AlbertaGridIdSet
  // ----------------

  //! hierarchic index set of AlbertaGrid
  template< int dim, int dimworld >
  class AlbertaGridIdSet
    : public IdSet< AlbertaGrid< dim, dimworld >, AlbertaGridIdSet< dim, dimworld >, unsigned int >
  {
    typedef AlbertaGridIdSet< dim, dimworld > This;
    typedef IdSet< AlbertaGrid< dim, dimworld >, This, unsigned int > Base;

    friend class AlbertaGrid< dim, dimworld >;

  public:
    //! export type of id
    typedef typename Base::IdType IdType;

  private:
    typedef AlbertaGrid< dim, dimworld > Grid;

    static const int dimension = Grid::dimension;

    typedef typename Grid::HierarchicIndexSet HierarchicIndexSet;

    // create id set, only allowed for AlbertaGrid
    AlbertaGridIdSet ( const HierarchicIndexSet &hIndexSet )
      : hIndexSet_( hIndexSet )
    {}

  public:
    /** \copydoc IdSet::id(const EntityType &e) const */
    template< class Entity >
    IdType id ( const Entity &e ) const
    {
      const int codim = Entity::codimension;
      return id< codim >( e );
    }

    /** \copydoc IdSet::id(const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity &e) const */
    template< int codim >
    IdType id ( const typename Grid::template Codim< codim >::Entity &e ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      const IdType index = hIndexSet_.index( e );
      return (index << 2) | IdType( codim );
    }

    /** \copydoc IdSet::subId(const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity &e,int i,unsigned int codim) const */
    IdType subId ( const typename Grid::template Codim< 0 >::Entity &e, int i, unsigned int subcodim ) const
    {
      assert( int( subcodim ) <= dimension );
      const IdType index = hIndexSet_.subIndex( e, i, subcodim );
      return (index << 2) | IdType( subcodim );
    }

    template< int codim >
    IdType subId ( const typename Grid::template Codim< codim >::Entity &e, int i, unsigned int subcodim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) && (int( codim + subcodim ) <= dimension) );
      const IdType index = hIndexSet_.subIndex( e, i, subcodim );
      return (index << 2) | IdType( codim + subcodim );
    }

    template< class Entity >
    IdType subId ( const Entity &e, int i, unsigned int subcodim ) const
    {
      return subId< Entity::codimension >( e, i, subcodim );
    }

  private:
    // prohibit copying
    AlbertaGridIdSet ( const This & );

    const HierarchicIndexSet &hIndexSet_;
  };

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTAGRIDINDEXSETS_HH
