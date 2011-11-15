// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/geometry/referenceelements.hh>

namespace Dune
{
  template <int d, int w>
  class AlbertaGrid;
  template <class G>
  struct HasLevelIntersections
  {
    static const bool value = true;
  };
  template <int d, int w>
  struct HasLevelIntersections<AlbertaGrid<d,w> >
  {
    static const bool value = false;
  };
  //****************************************************************
  //
  // --GrapeGridDisplay, GrapeGridDisplay for given grid
  //
  //****************************************************************

  template<class GridType>
  inline GrapeGridDisplay<GridType>::
  GrapeGridDisplay(const GridType &grid, const int myrank )
    :
#if HAVE_GRAPE
      setGridPartIter_(0)
      , entityIndex(GrapeGridDisplay<GridType>::template getEntityIndex<LeafIndexSetType>)
      , vertexIndex(GrapeGridDisplay<GridType>::template getVertexIndex<LeafIndexSetType>),
#endif
      grid_(grid)
      , hasLevelIntersections_(HasLevelIntersections<GridType>::value)
      , gridPart_(0)
      , indexSet_( (void *)(&grid.leafIndexSet()) )
      , lid_(grid.localIdSet())
      , myRank_(myrank)
      , hmesh_ (0)
  {
#if HAVE_GRAPE
    GrapeInterface<dim,dimworld>::init();
    if(!hmesh_) hmesh_ = setupHmesh();
#endif
  }

  template<class GridType>
  template<class GridPartType>
  inline GrapeGridDisplay<GridType>::
  GrapeGridDisplay(const GridPartType &gridPart, const int myrank )
    :
#if HAVE_GRAPE
      setGridPartIter_(&SetIter<GridPartType>::setGPIterator)
      , entityIndex(GrapeGridDisplay<GridType>::
                    template getEntityIndex<typename GridPartType::IndexSetType>)
      , vertexIndex(GrapeGridDisplay<GridType>::
                    template getVertexIndex<typename GridPartType::IndexSetType>),
#endif
      grid_(gridPart.grid())
      , hasLevelIntersections_(HasLevelIntersections<GridType>::value)
      , gridPart_((void *) &gridPart)
      , indexSet_( (void *)(&gridPart.indexSet()) )
      , lid_(grid_.localIdSet())
      , myRank_(myrank)
      , hmesh_ (0)
  {
#if HAVE_GRAPE
    GrapeInterface<dim,dimworld>::init();
    if(!hmesh_) hmesh_ = setupHmesh();
#endif
  }

  template< class GridType >
  template< class VT >
  inline GrapeGridDisplay< GridType >
  ::GrapeGridDisplay ( const GridView< VT > &gridView, const int myrank )
    :
#if HAVE_GRAPE
      setGridPartIter_( &GridViewIterators< VT >::set ),
      entityIndex( getEntityIndex< typename GridView< VT >::IndexSet > ),
      vertexIndex( getVertexIndex< typename GridView< VT >::IndexSet > ),
#endif
      grid_( gridView.grid() ),
      hasLevelIntersections_(HasLevelIntersections<GridType>::value),
      gridPart_( (void *)&gridView ),
      indexSet_( (void *)&gridView.indexSet() ),
      lid_( grid_.localIdSet() ),
      myRank_( myrank ),
      hmesh_( 0 )
  {
#if HAVE_GRAPE
    GrapeInterface< dim, dimworld >::init();
    if( !hmesh_ )
      hmesh_ = setupHmesh();
#endif
  }

  template<class GridType>
  inline GrapeGridDisplay<GridType>::
  ~GrapeGridDisplay()
  {
#if HAVE_GRAPE
    dune_.delete_iter(&hel_);
    ThisType::deleteStackEntry(stackEntry_);
    deleteHmesh();
#endif
  }

#if HAVE_GRAPE
  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  deleteStackEntry(StackEntryType & stackEntry)
  {
    while( !stackEntry.empty() )
    {
      STACKENTRY * entry = stackEntry.top();
      stackEntry.pop();

      DUNE_ELEM * elem = (DUNE_ELEM *) entry->hel.user_data;
      delete elem;
      delete entry;
    }
  }


  //****************************************************************
  //
  // --GridDisplay, Some Subroutines needed in display
  //
  //****************************************************************
  /** hmesh functionen **/

  template<class GridType>
  template <class IntersectionIteratorType>
  inline void GrapeGridDisplay<GridType>::
  checkNeighbors(IntersectionIteratorType & nit,
                 const IntersectionIteratorType & endnit,
                 DUNE_ELEM * he)
  {
    typedef typename GridType::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorType :: Intersection IntersectionType;
    int lastElNum = -1;

    // check all faces for boundary or not
    for( ; nit != endnit; ++nit )
    {
      const IntersectionType &intersection = *nit;
      const int number = mapDune2GrapeFace(he->type, intersection.indexInInside());
      assert( (number >= 0) && (number < MAX_EL_FACE) );

      if( number != lastElNum )
      {
        he->bnd[ number ]
          = (intersection.boundary() ? intersection.boundaryId() : 0);
        if( intersection.neighbor() )
        {
          if( intersection.outside()->partitionType() != InteriorEntity )
            he->bnd[ number ] = 2*(Entity :: dimensionworld) + (number+1);
        }
        lastElNum = number;
      }
    }
  }

  template<class GridType>
  template <class Entity>
  inline void GrapeGridDisplay<GridType>::
  el_update_base (Entity& en, DUNE_ELEM * he)
  {
    typedef typename Entity::Geometry DuneGeometryType;
    typedef typename DuneGeometryType :: ctype ctype;

    enum { dim      = Entity::dimension };
    enum { dimworld = Entity::dimensionworld };

    typedef FieldVector<ctype, dimworld> CoordinateType;

    const DuneGeometryType &geometry = en.geometry();

    he->eindex = this->entityIndex(indexSet_,en);
    he->level  = en.level();

    // if not true, only the macro level is drawn
    he->has_children = 1;

    // know the type
    int geomType = convertToGrapeType( geometry.type(), dim );
    he->type = geomType;

    // get pointer to coordinates and copy from geometry
    double ** vpointer = he->vpointer;

    // number of corners and number of vertices schould be the same
    // grape visual does not work for other situations
    assert( en.template count<dim>() == geometry.corners() );
    assert( geometry.corners() <= MAX_EL_DOF );

    for( int i = 0; i < geometry.corners(); ++i )
    {
      const int grapeVx = mapDune2GrapeVertex( geomType, i );
      he->vindex[ i ] = this->vertexIndex( indexSet_, en, grapeVx );

      assert( Entity::dimensionworld <= 3 );
      const CoordinateType coord = geometry.corner( grapeVx );
      for( int j = 0; j < Entity::dimensionworld ; ++j )
      {
        // here the mapping from dune to grape elements is done
        // it's only different for quads and hexas
        vpointer[ i ][ j ] = coord[ j ];
      }
    }
  }

  template<class GridType>
  template <class EntityPointerType>
  inline int GrapeGridDisplay<GridType>::
  el_update (EntityPointerType * it, DUNE_ELEM * he)
  {
    typedef typename GridType::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry DuneGeometryType;
    typedef typename DuneGeometryType :: ctype ctype;

    enum { dim      = Entity::dimension };
    enum { dimworld = Entity::dimensionworld };

    typedef FieldVector<ctype, dimworld> CoordinateType;

    Entity &en = (*it[0]);

    // only for debuging, becsaue normaly references are != NULL
    if(&en)
    {
      el_update_base( en, he );

      {
        if( en.hasBoundaryIntersections() )
        {
          // reset the boundary information
          for(int i=0; i < MAX_EL_FACE; ++i) he->bnd[i] = -1;

          if( en.isLeaf() )
          {
            typedef typename Entity::LeafIntersectionIterator IntersectionIterator;
            IntersectionIterator endnit = en.ileafend();
            IntersectionIterator nit    = en.ileafbegin();

            checkNeighbors(nit,endnit,he);
          }
          else if(hasLevelIntersections_)
          {
            typedef typename Entity::LevelIntersectionIterator IntersectionIterator;
            IntersectionIterator endnit = en.ilevelend();
            IntersectionIterator nit    = en.ilevelbegin();

            checkNeighbors(nit,endnit,he);
          }
        }
        else
        {
          // if no boundary intersections, then all faces are interior
          for(int i=0; i < MAX_EL_FACE; ++i) he->bnd[i] = 0;
        }
      }

      // for data displaying
      he->actElement = it;
      return 1;

    } // end if(&en)
    else
    {
      he->actElement = 0;
      return 0;
    }
  }

  template<class GridType>
  template <class EntityPointerType, class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  el_update (EntityPointerType * it, DUNE_ELEM * he, GridPartType& gridPart)
  {
    typedef typename GridType::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry DuneGeometryType;
    typedef typename DuneGeometryType :: ctype ctype;

    enum { dim      = Entity::dimension };
    enum { dimworld = Entity::dimensionworld };

    typedef FieldVector<ctype, dimworld> CoordinateType;

    const Entity &en = (*it[0]);

    // only for debuging, becsaue normaly references are != NULL
    if(&en)
    {
      el_update_base ( en , he );

      {
        typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
        // reset the boundary information
        for(int i=0; i < MAX_EL_FACE; ++i) he->bnd[i] = -1;

        IntersectionIteratorType endnit = gridPart.iend(en);
        IntersectionIteratorType nit = gridPart.ibegin(en);

        checkNeighbors(nit,endnit,he);
      }

      // for data displaying
      he->actElement = it;
      return 1;

    } // end if(&en)
    else
    {
      he->actElement = 0;
      return 0;
    }
  }


  template< class GridType >
  template< class EntityPointer, class VT >
  inline int GrapeGridDisplay< GridType >
  ::el_update ( EntityPointer *it, DUNE_ELEM *he, const GridView< VT > &gridView )
  {
    typedef typename GridView< VT >::template Codim< 0 >::Entity Entity;
    typedef typename Entity::Geometry DuneGeometryType;
    typedef typename DuneGeometryType::ctype ctype;

    //const int dim      = Entity::dimension;
    const int dimworld = Entity::dimensionworld;

    typedef FieldVector< ctype, dimworld > CoordinateType;

    Entity &entity = **it;
    assert( &entity != 0 );

    el_update_base ( entity, he );

    typedef typename GridView< VT >::IntersectionIterator IntersectionIterator;
    // reset the boundary information
    for( int i = 0; i < MAX_EL_FACE; ++i )
      he->bnd[ i ] = -1;

    IntersectionIterator iend = gridView.iend( entity );
    IntersectionIterator iit  = gridView.ibegin( entity );

    checkNeighbors( iit, iend, he );

    // for data displaying
    he->actElement = it;
    return 1;
  }


  template<class GridType>
  template<PartitionIteratorType pitype>
  inline int GrapeGridDisplay<GridType>::
  first_leaf (DUNE_ELEM * he)
  {
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LeafIterator LeafIteratorType;

    if(he->liter) dune_.delete_iter(he);

    he->liter   = 0;
    he->enditer = 0;

    LeafIteratorType * it    = new LeafIteratorType ( grid_.template leafbegin<0, pitype> () );
    LeafIteratorType * endit = new LeafIteratorType ( grid_.template leafend  <0, pitype> () );

    he->liter   = (void *) it;
    he->enditer = (void *) endit;

    if(it[0] == endit[0])
    {
      this->template delete_leaf<pitype>(he);
      return 0;
    }

    return el_update(it,he);
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline int GrapeGridDisplay<GridType>::
  next_leaf (DUNE_ELEM * he)
  {
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LeafIterator LeafIteratorType;

    LeafIteratorType * it    = (LeafIteratorType *) he->liter;
    LeafIteratorType * endit = (LeafIteratorType *) he->enditer;
    assert( it );
    assert( endit );

    if( ++it[0] != endit[0] )
    {
      return el_update(it,he);
    }
    else
    {
      this->template delete_leaf<pitype> (he);
    }
    return 0;
  }

  template<class GridType>
  template<class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  first_item (DUNE_ELEM * he)
  {
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;

    if(he->liter) dune_.delete_iter(he);

    assert( he->gridPart );
    GridPartType & gridPart = *((GridPartType *) he->gridPart);

    assert( he->liter   == 0 );
    assert( he->enditer == 0 );

    IteratorType * it    = new IteratorType ( gridPart.template begin<0> () );
    IteratorType * endit = new IteratorType ( gridPart.template end  <0> () );

    he->liter   = (void *) it;
    he->enditer = (void *) endit;

    if(it[0] == endit[0])
    {
      this->template delete_iterators<IteratorType> (he);
      return 0;
    }

    return el_update(it,he,gridPart);
  }

  template<class GridType>
  template<class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  next_item (DUNE_ELEM * he)
  {
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;

    assert( he->gridPart );
    GridPartType & gridPart = *((GridPartType *) he->gridPart);

    IteratorType * it    = (IteratorType *) he->liter;
    IteratorType * endit = (IteratorType *) he->enditer;
    assert( it );
    assert( endit );

    if( ++it[0] != endit[0] )
    {
      return el_update(it,he,gridPart);
    }
    else
    {
      this->template delete_iterators<IteratorType> (he);
    }
    return 0;
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline int GrapeGridDisplay<GridType>::
  first_level (DUNE_ELEM * he, int lvl)
  {
    if(he->liter) dune_.delete_iter(he);

    assert( he->liter   == 0 );
    assert( he->enditer == 0 );

    // for leaf level, lvl has the value -1
    int level = (lvl < 0) ? grid_.maxLevel() : lvl;

    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LevelIterator LevelIteratorType;

    // class copy constructor
    LevelIteratorType * it    = new LevelIteratorType( grid_.template lbegin<0,pitype> (level) );
    LevelIteratorType * endit = new LevelIteratorType( grid_.template lend<0,pitype>   (level) );

    he->liter   = (void *) it;
    he->enditer = (void *) endit;

    if(it[0] == endit[0])
    {
      this->template delete_level<pitype>(he);
      return 0;
    }

    return el_update(it,he);
  }


  template<class GridType>
  template<PartitionIteratorType pitype>
  inline int GrapeGridDisplay<GridType>::
  next_level (DUNE_ELEM * he)
  {
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LevelIterator LevelIteratorType;

    LevelIteratorType * it    = ((LevelIteratorType *) he->liter);
    LevelIteratorType * endit = ((LevelIteratorType *) he->enditer);

    assert( it );
    assert( endit );
    if( ++it[0] != endit[0] )
    {
      return el_update(it,he);
    }
    else
    {
      this->template delete_level<pitype>(he);

      // clear all hierachic iterators
      while(!hierList_.empty())
      {
        HierarchicIteratorType * hit = hierList_.back();
        hierList_.pop_back();
        delete hit;
      }
      assert( hierList_.size () == 0 );
    }
    return 0;
  }


  template<class GridType>
  template<class EntityPointerType>
  inline int GrapeGridDisplay<GridType>::
  child_update(EntityPointerType * it, DUNE_ELEM * he)
  {
    typedef typename  GridType :: template Codim<0> :: Entity EntityType;

    const EntityType & en = (*it[0]);

    HierarchicIteratorType * hit = (HierarchicIteratorType *) he->hiter;

    const EntityType *newEn = (!hit) ? (&en) : (hit[0].operator -> ()) ;

    assert( newEn );

    // if entity is leaf, then no first child
    if( newEn->isLeaf() ) return 0;

    int childLevel = newEn->level() + 1;

    // store former pointer for removal later
    if(hit) hierList_.push_back( hit );

    // create HierarchicIterator with default constructor
    hit = new HierarchicIteratorType ( newEn->hbegin ( childLevel ) );

    assert( hit != 0 );
    if( hit[0] != newEn->hend( childLevel ) )
    {
      he->hiter = (void *) hit;
      return el_update( hit, he);
    }
    else
    {
      hierList_.pop_back();
      delete hit;
      return 0;
    }
  }

  template<class GridType>
  template<class EntityPointerType>
  inline int GrapeGridDisplay<GridType>::
  child_n_update(EntityPointerType *it, DUNE_ELEM * he)
  {
    typedef typename  GridType::Traits::template Codim<0>::Entity EntityType;

    const EntityType &en = (*it[0]);

    int childLevel = en.level();
    HierarchicIteratorType * hit = (HierarchicIteratorType *) he->hiter;
    assert( hit );

    //HierarchicIteratorType ehit = hit[0]->hend(childLevel);
    HierarchicIteratorType ehit = en.hend(childLevel);
    if( ++hit[0] != ehit )
    {
      return el_update(hit,he);
    }

    hierList_.remove( hit );
    delete hit;
    he->hiter = 0;

    return 0;
  }



  template<class GridType>
  template<class IteratorType>
  inline void GrapeGridDisplay<GridType>::
  delete_iterators(DUNE_ELEM * he)
  {
    IteratorType * it  = ((IteratorType *) he->liter);
    if(it)
    {
      IteratorType * endit = ((IteratorType *) he->enditer);
      assert( endit );

      delete it;
      delete endit;

      he->actElement = 0;
      he->liter      = 0;
      he->enditer    = 0;
    }
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline void GrapeGridDisplay<GridType>::
  delete_leaf (DUNE_ELEM * he)
  {
    assert( he );
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LeafIterator IteratorType;

    this->template delete_iterators<IteratorType> (he);
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline void GrapeGridDisplay<GridType>::
  delete_level (DUNE_ELEM * he)
  {
    assert( he );
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LevelIterator IteratorType;

    this->template delete_iterators<IteratorType> (he);
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline void GrapeGridDisplay<GridType>::
  delete_hier (DUNE_ELEM * he)
  {
    assert( he );
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LevelIterator IteratorType;

    this->template delete_iterators<IteratorType> (he);

    // clear all hierachic iterators
    while(!hierList_.empty())
    {
      HierarchicIteratorType * hit = hierList_.back();
      hierList_.pop_back();
      delete hit;
    }
    assert( hierList_.size () == 0 );
  }


  template<class GridType>
  inline int GrapeGridDisplay<GridType>::
  first_child(DUNE_ELEM * he)
  {
    typedef typename GridType :: template Codim<0> ::
    EntityPointer EntityPointerType;
    return child_update( ((EntityPointerType *) he->liter) , he);
  }


  template<class GridType>
  inline int GrapeGridDisplay<GridType>::
  next_child(DUNE_ELEM * he)
  {
    typedef typename GridType :: template Codim<0> ::
    EntityPointer EntityPointerType;
    return child_n_update( ((EntityPointerType *) he->liter), he);
  }


  template<class GridType>
  inline void * GrapeGridDisplay<GridType>::
  copy_iterator (const void * i)
  {
    std::cerr << "ERROR: copt_iterator not implemented! file = " << __FILE__ << ", line = " << __LINE__ << "\n";
    DUNE_THROW(NotImplemented,"method copy_iterator not implemented!");
    return 0 ;
  }

  // checkInside
  template< class GridType >
  template< class Entity >
  inline int GrapeGridDisplay< GridType >
  ::checkInside( const Entity &entity, const double *c )
  {
    const int dim = Entity::dimension;

    FieldVector< double, dim > local;
    for( int i = 0; i < dim; ++i )
      local[ i ] = c[ i ];

    // see hmesh doc page 32, if point is inside, -1 has to be returned
    // otherwise local face , grrrr
    const GenericReferenceElement< double, dim > &refElement
      = GenericReferenceElements< double, dim >::general( entity.type() );
    return (refElement.checkInside( local ) ? -1 : 0);
  }

  // check inside
  template<class GridType>
  inline int GrapeGridDisplay<GridType>::
  checkWhetherInside(DUNE_ELEM * he, const double * w)
  {
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    EntityPointerType * ep = (EntityPointerType *) he->actElement;
    assert( ep );
    return checkInside(*(ep[0]),w);
  }

  // world to local
  template<class GridType> template <class EntityType>
  inline void GrapeGridDisplay<GridType>::
  local_to_world(const EntityType &en, const double * c, double * w)
  {
    const int dim = EntityType::dimension;
    const int dimworld = EntityType::dimensionworld;

    FieldVector< double, dim > local;
    for( int i = 0; i < dim; ++i )
      local[ i ] = c[ i ];

    FieldVector< double, dimworld > global = en.geometry().global( local );
    for( int i = 0; i < dimworld; ++i )
      w[ i ] = global[ i ];
  }


  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  local2world (DUNE_ELEM * he, const double * c, double * w)
  {
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    EntityPointerType * ep = (EntityPointerType *) he->actElement;
    assert( ep );
    local_to_world(*(ep[0]),c,w);
    return;
  }

  // world to local
  template< class GridType >
  template< class Entity >
  inline int GrapeGridDisplay< GridType >
  ::world_to_local( const Entity &entity, const double *w, double *c )
  {
    const int dim = Entity::dimension;
    const int dimworld = Entity::dimensionworld;

    const typename Entity::Geometry &geometry = entity.geometry();

    FieldVector< double, dimworld > global;
    for( int i = 0; i < dimworld; ++i )
      global[ i ] = w[ i ];
    FieldVector< double, dim > local = geometry.local( global );
    for( int i = 0; i < dim; ++i )
      c[ i ] = local[ i ];

    const GenericReferenceElement< double, dim > &refElement
      = GenericReferenceElements< double, dim >::general( geometry.type() );
    return (refElement.checkInside( local ) ? -1 : 0);
  }

  // world to local
  template<class GridType>
  inline int GrapeGridDisplay<GridType>::
  world2local(DUNE_ELEM * he, const double * w, double * c)
  {
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    EntityPointerType * ep = (EntityPointerType *) he->actElement;
    assert( ep );
    return world_to_local(*(ep[0]),w,c);
  }


  // world to local
  template<class GridType>
  template<PartitionIteratorType pitype>
  inline void GrapeGridDisplay<GridType>::
  selectIterators(DUNE_DAT * dune, void * gridPart, setGridPartIterators_t * func) const
  {
    // if pointer are 0, then no evaluation is done
    dune->first_child = 0;
    dune->next_child  = 0;

    if(dune->iteratorType == g_LeafIterator)
    {
      dune->first_macro = &IterationMethods<pitype>::fst_leaf;
      dune->next_macro  = &IterationMethods<pitype>::nxt_leaf;
      dune->delete_iter = &IterationMethods<pitype>::del_leaf;

      return ;
    }

    if(dune->iteratorType == g_LevelIterator)
    {
      dune->first_macro = &IterationMethods<pitype>::first_lev;
      dune->next_macro  = &IterationMethods<pitype>::next_lev;
      dune->delete_iter = &IterationMethods<pitype>::del_level;

      return ;
    }

    if(dune->iteratorType == g_HierarchicIterator)
    {
      dune->first_macro = &IterationMethods<pitype>::first_mac;
      dune->next_macro  = &IterationMethods<pitype>::next_lev;
      dune->delete_iter = &IterationMethods<pitype>::del_hier;

      dune->first_child = &IterationMethods<pitype>::fst_child;
      dune->next_child  = &IterationMethods<pitype>::nxt_child;

      return ;
    }

    if(dune->iteratorType == g_GridPart)
    {
      static bool called = false;
      if(!func)
      {
        std::string name("Null");
        if(!called)
          std::cerr << "No function for data '" <<name<<"' and therefore no GridPart! Defaulting Iterator to LeafIterator! \n";

        dune->first_macro = &IterationMethods<pitype>::fst_leaf;
        dune->next_macro  = &IterationMethods<pitype>::nxt_leaf;

        dune->gridPart = 0;
        called = true;
        return ;
      }

      assert( func );
      assert( gridPart );
      // set first and next methods due to grid part
      func(dune,gridPart);
      called = false;

      return ;
    }

    // wrong iteratorType here
    assert(false);
    abort();
  }


  // setIterationsMethods
  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  setIterationMethods(DUNE_DAT * dune, DUNE_FDATA * data) const
  {
    if(dune->delete_iter) dune->delete_iter(dune->all);

    void * gridPart = 0;
    setGridPartIterators_t * func = 0;

    if(data)
    {
      gridPart = data->gridPart;
      func = data->setGridPartIterators;
    }
    else if(gridPart_ && setGridPartIter_)
    {
      gridPart = gridPart_;
      func = setGridPartIter_;
    }

    switch(dune->partitionIteratorType)
    {
    case g_All_Partition :            selectIterators<All_Partition> (dune,gridPart,func) ;
      return;
    case g_Interior_Partition :       selectIterators<Interior_Partition> (dune,gridPart,func) ;
      return;
    case g_InteriorBorder_Partition : selectIterators<InteriorBorder_Partition> (dune,gridPart,func) ;
      return;
    case g_Overlap_Partition :        selectIterators<Overlap_Partition> (dune,gridPart,func) ;
      return;
    case g_OverlapFront_Partition :   selectIterators<OverlapFront_Partition> (dune,gridPart,func) ;
      return;
    case g_Ghost_Partition :          selectIterators<Ghost_Partition> (dune,gridPart,func) ;
      return;
    default : assert(false);
      abort();
      return ;
    }
  }

  // setIterationsMethods
  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  changeIterationMethods(int iterType, int partType, DUNE_FDATA * data )
  {
    dune_.iteratorType = iterType;
    dune_.partitionIteratorType = partType;
    setIterationMethods(&dune_,data);
  }


  // check inside
  template<class GridType>
  inline int GrapeGridDisplay<GridType>::
  check_inside(DUNE_ELEM * he, const double * w)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    return disp[0].checkWhetherInside(he,w);
  }
  // local to world
  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  ctow (DUNE_ELEM * he, const double * c, double * w)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    disp[0].local2world(he,c,w);
    return ;
  }

  // world to local
  template<class GridType>
  inline int GrapeGridDisplay<GridType>::
  wtoc(DUNE_ELEM * he, const double * w, double * c)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    return disp[0].world2local(he,w,c);
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  setIterationModus(DUNE_DAT * dat, DUNE_FDATA * func)
  {
    MyDisplayType * disp = (MyDisplayType *) dat->all->display;
    disp[0].setIterationMethods(dat,func);
  }

  template<class GridType>
  inline void * GrapeGridDisplay<GridType>::getHmesh()
  {
    if(!hmesh_) hmesh_ = setupHmesh();
    return (void *) hmesh_;
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  addMyMeshToTimeScene(void * timescene, double time, int proc)
  {
    GrapeInterface<dim,dimworld>::addHmeshToTimeScene(timescene,time,this->getHmesh(),proc);
  }

  template<class GridType>
  inline void * GrapeGridDisplay<GridType>::setupHmesh()
  {
    int maxlevel = grid_.maxLevel();

    int noe = grid_.size(0);
    int nov = grid_.size(dim);

    // set pointer to me
    hel_.display = (void *) this;

    // set dune pointers
    DUNE_DAT * dune = &dune_;

    dune->wtoc         = wtoc;
    dune->ctow         = ctow;
    dune->check_inside = check_inside;

    // set method to select iterators
    dune->setIterationModus = &setIterationModus;

    dune->get_stackentry = &getStackEn;
    dune->free_stackentry = &freeStackEn;

    dune->all          = &hel_;
    dune->partition    = myRank_;

    dune->iteratorType          = g_LeafIterator;
    dune->partitionIteratorType = g_All_Partition;

    setIterationMethods(dune,0);

    std::string gridName( "Grid" );
    /* return hmesh with no data */
    return GrapeInterface<dim,dimworld>::
           setupHmesh(noe,nov,maxlevel,dune, gridName.c_str());
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::deleteHmesh()
  {
    if( hmesh_ )
    {
      GrapeInterface<dim,dimworld>::deleteHmesh(hmesh_);
    }
  }


  template<class GridType>
  inline void * GrapeGridDisplay<GridType>::
  getStackEntry(StackEntryType & stackEntry)
  {
    STACKENTRY * entry = 0;

    if( stackEntry.empty() )
    {
      entry = new STACKENTRY ();
      DUNE_ELEM * elem = new DUNE_ELEM ();
      assert( elem );
      entry->hel.user_data = (void *)elem;
    }
    else
    {
      entry = stackEntry.top();
      stackEntry.pop();
    }
    assert( entry );
    return( (void *) entry);
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  freeStackEntry(StackEntryType & stackEntry, void * entry)
  {
    assert( entry );
    stackEntry.push( ((STACKENTRY *) entry) );
  }

  template<class GridType>
  inline void * GrapeGridDisplay<GridType>::
  getStackEn(DUNE_DAT * dune)
  {
    MyDisplayType * disp = (MyDisplayType *) dune->all->display;
    return MyDisplayType::getStackEntry(disp->stackEntry_);
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  freeStackEn(DUNE_DAT * dune, void * entry)
  {
    MyDisplayType * disp = (MyDisplayType *) dune->all->display;
    MyDisplayType::freeStackEntry(disp->stackEntry_,entry);
  }

#endif

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::display()
  {
#if HAVE_GRAPE
    /* call handle mesh in g_hmesh.c */
    GrapeInterface<dim,dimworld>::handleMesh ( hmesh_ , true );
#endif
    return ;
  }

  template<class GridType>
  inline const GridType & GrapeGridDisplay<GridType>::getGrid() const
  {
    return grid_;
  }

} // end namespace Dune
