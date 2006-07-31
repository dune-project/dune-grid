// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{

  //****************************************************************
  //
  // --GrapeGridDisplay, GrapeGridDisplay for given grid
  //
  //****************************************************************

  template<class GridType>
  inline GrapeGridDisplay<GridType>::
  GrapeGridDisplay(const GridType &grid, const int myrank )
    : grid_(grid)
      , indexSet_( (void *)(&grid.leafIndexSet()) )
      , lid_(grid.localIdSet())
      , myRank_(myrank)
      , hmesh_ (0)
      , entityIndex(GrapeGridDisplay<GridType>::template getEntityIndex<LeafIndexSetType>)
      , vertexIndex(GrapeGridDisplay<GridType>::template getVertexIndex<LeafIndexSetType>)
  {
    GrapeInterface<dim,dimworld>::init();
    if(!hmesh_) hmesh_ = setupHmesh();
  }

  template<class GridType>
  template<class GridPartType>
  inline GrapeGridDisplay<GridType>::
  GrapeGridDisplay(const GridPartType &gridPart, const int myrank )
    : grid_(gridPart.grid())
      , indexSet_( (void *)(&gridPart.indexSet()) )
      , lid_(grid_.localIdSet())
      , myRank_(myrank)
      , hmesh_ (0)
      , entityIndex(GrapeGridDisplay<GridType>::
                    template getEntityIndex<typename GridPartType::IndexSetType>)
      , vertexIndex(GrapeGridDisplay<GridType>::
                    template getVertexIndex<typename GridPartType::IndexSetType>)
  {
    GrapeInterface<dim,dimworld>::init();
    if(!hmesh_) hmesh_ = setupHmesh();
  }

  template<class GridType>
  inline GrapeGridDisplay<GridType>::
  ~GrapeGridDisplay()
  {
    //deleteHmesh();
  }

  //****************************************************************
  //
  // --GridDisplay, Some Subroutines needed in display
  //
  //****************************************************************
  /** hmesh functionen **/

  template<class GridType>
  template <class EntityPointerType>
  inline int GrapeGridDisplay<GridType>::
  el_update (EntityPointerType * it, DUNE_ELEM * he)
  {
    typedef typename GridType::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::IntersectionIterator IntersectionIterator;
    typedef typename Entity::Geometry DuneGeometryType;

    enum { dim      = Entity::dimension };
    enum { dimworld = Entity::dimensionworld };

    Entity &en = (*it[0]);

    // only for debuging, becsaue normaly references are != NULL
    if(&en)
    {
      const DuneGeometryType &geometry = en.geometry();

      he->eindex = this->entityIndex(indexSet_,en);
      he->level  = en.level();

      // if not true, only the macro level is drawn
      he->has_children = 1;

      // know the type
      int geomType = convertToGrapeType ( geometry.type() , dim );
      he->type = geomType;

      {
        // set the vertex coordinates
        double (* vpointer)[3] = he->vpointer;

        // number of corners and number of vertices schould be the same
        // grape visual does not work for other situations
        assert( en.template count<dim>() == geometry.corners() );

        for(int i= 0 ; i<geometry.corners(); ++i)
        {
          const int grapeVx = mapDune2GrapeVertex(geomType,i);
          he->vindex[i] = this->vertexIndex(indexSet_, en, grapeVx);

          for(int j = 0; j < Entity::dimensionworld ; ++j)
          {
            // here the mapping from dune to grape elements is done
            // it's only different for quads and hexas
            vpointer[i][j] = geometry[ grapeVx ][j] ;
          }
        }
      } // end set all vertex coordinates

      {
        // reset the boundary information
        for(int i=0; i < MAX_EL_FACE; i++) he->bnd[i] = -1;

        IntersectionIterator endnit = en.iend();
        IntersectionIterator nit    = en.ibegin();

        // value < zero otherwise first test fails
        int lastElNum = -1;

        // check all faces for boundary or not
        while ( nit != endnit )
        {
          int num = nit.numberInSelf();
          assert( num >= 0 );
          assert( num < MAX_EL_FACE );

          if(num != lastElNum)
          {
            he->bnd[num] = ( nit.boundary() ) ? nit.boundaryId() : 0;
            //if(nit.levelNeighbor() || nit.leafNeighbor() )
            if( nit.neighbor() )
              if(nit.outside()->partitionType() != InteriorEntity )
                he->bnd[num] = 2*(Entity::dimensionworld) + (nit.numberInSelf()+1);
            lastElNum = num;
          }
          ++nit;
        }
      }

      {
        // for this type of element we have to swap the faces
        if(he->type == g_hexahedron)
        {
          int help_bnd [MAX_EL_FACE];
          for(int i=0; i < MAX_EL_FACE; i++) help_bnd[i] = he->bnd[i] ;

          assert( MAX_EL_FACE == 6 );
          // do the mapping from dune to grape hexa
          he->bnd[0] = help_bnd[4];
          he->bnd[1] = help_bnd[5];
          he->bnd[3] = help_bnd[1];
          he->bnd[4] = help_bnd[3];
          he->bnd[5] = help_bnd[0];
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
  template<PartitionIteratorType pitype, class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  first_entity (DUNE_ELEM * he)
  {
    GridPartType & gridPart = *((GridPartType*) he->gridPart);
    typedef typename GridPartType :: Traits:: template Codim<0> ::
    IteratorType IteratorType;

    he->liter   = 0;
    he->enditer = 0;

    IteratorType * it    = new IteratorType ( gridPart.template begin<0>() );
    IteratorType * endit = new IteratorType ( gridPart.template end  <0>() );

    if(it[0] == endit[0])
    {
      he->actElement = 0;
      delete it;
      delete endit;
      return 0;
    }

    he->liter   = (void *) it;
    he->enditer = (void *) endit;
    return el_update(it,he);
  }

  template<class GridType>
  template<PartitionIteratorType pitype, class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  next_entity (DUNE_ELEM * he)
  {
    typedef typename GridPartType :: Traits:: template Codim<0> ::
    IteratorType IteratorType;

    IteratorType * it    = (IteratorType *) he->liter;
    IteratorType * endit = (IteratorType *) he->enditer;
    assert( it );
    assert( endit );

    if( ++it[0] != endit[0] )
    {
      return el_update(it,he);
    }
    else
    {
      delete it;
      delete endit;
      he->liter      = 0;
      he->enditer    = 0;
      he->actElement = 0;
    }
    return 0;
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline int GrapeGridDisplay<GridType>::
  first_leaf (DUNE_ELEM * he)
  {
    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LeafIterator LeafIteratorType;

    he->liter   = 0;
    he->enditer = 0;

    LeafIteratorType * it    = new LeafIteratorType ( grid_.template leafbegin<0, pitype> () );
    LeafIteratorType * endit = new LeafIteratorType ( grid_.template leafend  <0, pitype> () );

    if(it[0] == endit[0])
    {
      he->actElement = 0;
      delete it;
      delete endit;
      return 0;
    }

    he->liter   = (void *) it;
    he->enditer = (void *) endit;
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
      delete it;
      delete endit;
      he->liter      = 0;
      he->enditer    = 0;
      he->actElement = 0;
    }
    return 0;
  }

  template<class GridType>
  template<class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  first_item (DUNE_ELEM * he)
  {
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;

    GridPartType & gridPart = *((GridPartType *) he->gridPart);

    he->liter   = 0;
    he->enditer = 0;

    IteratorType * it    = new IteratorType ( gridPart.template begin<0> () );
    IteratorType * endit = new IteratorType ( gridPart.template end  <0> () );

    if(it[0] == endit[0])
    {
      he->actElement = 0;
      delete it;
      delete endit;
      return 0;
    }

    he->liter   = (void *) it;
    he->enditer = (void *) endit;
    return el_update(it,he);
  }

  template<class GridType>
  template<class GridPartType>
  inline int GrapeGridDisplay<GridType>::
  next_item (DUNE_ELEM * he)
  {
    typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;

    IteratorType * it    = (IteratorType *) he->liter;
    IteratorType * endit = (IteratorType *) he->enditer;
    assert( it );
    assert( endit );

    if( ++it[0] != endit[0] )
    {
      return el_update(it,he);
    }
    else
    {
      delete it;
      delete endit;
      he->liter      = 0;
      he->enditer    = 0;
      he->actElement = 0;
    }
    return 0;
  }

  template<class GridType>
  template<PartitionIteratorType pitype>
  inline int GrapeGridDisplay<GridType>::
  first_level (DUNE_ELEM * he, int lvl)
  {
    he->liter   = 0;
    he->enditer = 0;

    // for leaf level, lvl has the value -1
    int level = (lvl < 0) ? grid_.maxLevel() : lvl;

    typedef typename GridType :: template Codim<0> ::
    template Partition<pitype> :: LevelIterator LevelIteratorType;

    // class copy constructor
    LevelIteratorType * it    = new LevelIteratorType( grid_.template lbegin<0,pitype> (level) );
    LevelIteratorType * endit = new LevelIteratorType( grid_.template lend<0,pitype>   (level) );

    if(it[0] == endit[0])
    {
      he->actElement = 0;
      delete it;
      delete endit;
      return 0;
    }

    he->liter   = (void *) it;
    he->enditer = (void *) endit;
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
      delete it;
      delete endit;
      he->liter   = 0;
      he->enditer = 0;
      he->actElement = 0;
      he->hiter = 0;

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

    EntityType & en = (*it[0]);

    HierarchicIteratorType * hit = (HierarchicIteratorType *) he->hiter;

    EntityType *newEn = (!hit) ? (&en) : (hit[0].operator -> ()) ;

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

    EntityType &en = (*it[0]);

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
  template<class GridType> template <class EntityType>
  inline int GrapeGridDisplay<GridType>::
  checkInside(EntityType &en, const double * c)
  {
    enum { dim = EntityType::dimension };

    for(int i=0; i<dim; i++) localVec_[i] = c[i];

    // see hmesh doc page 32, if point is inside, -1 has to be returned
    // otherwise local face , grrrr
    int isInside = (en.geometry().checkInside(localVec_) == true) ? -1 : 0;

    return isInside;
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
  local_to_world(EntityType &en, const double * c, double * w)
  {
    enum { dim      = EntityType::dimension };
    enum { dimworld = EntityType::dimensionworld };

    for(int i=0; i<dim; i++) localVec_[i] = c[i];

    globalVec_ = en.geometry().global(localVec_);

    for(int i=0; i<dimworld; i++) w[i] = globalVec_[i];
    return;
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
  template<class GridType> template <class EntityType>
  inline int GrapeGridDisplay<GridType>::
  world_to_local(EntityType &en, const double * w, double * c)
  {
    enum { dim      = EntityType::dimension };
    enum { dimworld = EntityType::dimensionworld };

    for(int i=0; i<dimworld; i++) globalVec_[i] = w[i];

    localVec_ = en.geometry().local(globalVec_);

    for(int i=0; i<dim; ++i) c[i] = localVec_[i];

    return (en.geometry().checkInside(localVec_) == true) ? -1 : 0;
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
  selectIterators(DUNE_DAT * dune, DUNE_FUNC * func) const
  {
    if(dune->iteratorType == g_LeafIterator)
    {
      dune->first_macro = &IterationMethods<pitype>::fst_leaf;
      dune->next_macro  = &IterationMethods<pitype>::nxt_leaf;

      // if pointer are 0, then nor evaluation is done
      dune->first_child = 0;
      dune->next_child  = 0;

      return ;
    }

    if(dune->iteratorType == g_LevelIterator)
    {
      dune->first_macro = &IterationMethods<pitype>::first_lev;
      dune->next_macro  = &IterationMethods<pitype>::next_lev;

      dune->first_child = 0;
      dune->next_child  = 0;

      return ;
    }

    if(dune->iteratorType == g_HierarchicIterator)
    {
      dune->first_macro = &IterationMethods<pitype>::first_mac;
      dune->next_macro  = &IterationMethods<pitype>::next_lev;

      dune->first_child = &IterationMethods<pitype>::fst_child;
      dune->next_child  = &IterationMethods<pitype>::nxt_child;

      return ;
    }

    if(dune->iteratorType == g_GridPart)
    {
      bool validFunction = (func) ? (func->all) ? true : false : false;
      if(validFunction)
      {
        validFunction = (func->all->setGridPartIterators) ? true : false;
      }

      if(!validFunction)
      {
        std::cerr << "No function or function data, therefore no GridPart! Defaulting Iterator to LeafIterator! \n";
        dune->first_macro = &IterationMethods<pitype>::fst_leaf;
        dune->next_macro  = &IterationMethods<pitype>::nxt_leaf;

        // if pointer are 0, then nor evaluation is done
        dune->first_child = 0;
        dune->next_child  = 0;

        return ;
      }

      assert( func );
      assert( func->all );
      assert( func->all->setGridPartIterators );

      // set first and next methods due to grid part
      func->all->setGridPartIterators(dune,func->all->gridPart);

      return ;
    }

    // wrong iteratorType here
    assert(false);
    abort();
  }


  // setIterationsMethods
  template<class GridType>
  inline void GrapeGridDisplay<GridType>::
  setIterationMethods(DUNE_DAT * dune, DUNE_FUNC * func ) const
  {
    switch(dune->partitionIteratorType)
    {
    case g_All_Partition :            selectIterators<All_Partition> (dune,func) ;
      return;
    case g_Interior_Partition :       selectIterators<Interior_Partition> (dune,func) ;
      return;
    case g_InteriorBorder_Partition : selectIterators<InteriorBorder_Partition> (dune,func) ;
      return;
    case g_Overlap_Partition :        selectIterators<Overlap_Partition> (dune,func) ;
      return;
    case g_OverlapFront_Partition :   selectIterators<OverlapFront_Partition> (dune,func) ;
      return;
    case g_Ghost_Partition :          selectIterators<Ghost_Partition> (dune,func) ;
      return;
    default : assert(false);
      abort();
      return ;
    }
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
  setIterationModus(DUNE_DAT * dat, DUNE_FUNC * func)
  {
    MyDisplayType * disp = (MyDisplayType *) dat->all->display;
    disp[0].setIterationMethods(dat,func);
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::display()
  {
    /* call handle mesh in g_hmesh.c */
    GrapeInterface<dim,dimworld>::handleMesh ( hmesh_ , true );
    return ;
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
  inline const GridType & GrapeGridDisplay<GridType>::getGrid() const
  {
    return grid_;
  }

  template<class GridType>
  inline void * GrapeGridDisplay<GridType>::setupHmesh()
  {
    // default set all coordinates to zero
    for(int i=0; i<MAX_EL_DOF; i++)
      for(int j=0; j<3; j++)
      {
        hel_.vpointer[i][j] = 0.0;
      }

    int maxlevel = grid_.maxLevel();

    int noe = grid_.size(0);
    int nov = grid_.size(dim);

    hel_.display = (void *) this;
    hel_.liter   = 0;
    hel_.enditer = 0;

    hel_.hiter    = 0;

    hel_.actElement = 0;

    DUNE_DAT * dune = &dune_;

    dune->copy         = 0; // no copy at the moment
    dune->wtoc         = wtoc;
    dune->ctow         = ctow;
    dune->check_inside = check_inside;

    // set method to select iterators
    dune->setIterationModus = &setIterationModus;

    dune->all          = &hel_;
    dune->partition    = myRank_;

    dune->iteratorType          = g_LeafIterator;
    dune->partitionIteratorType = g_All_Partition;

    setIterationMethods(dune,0);

    /* return hmesh with no data */
    return GrapeInterface<dim,dimworld>::setupHmesh(0,noe,nov,maxlevel,0,dune);
  }

  template<class GridType>
  inline void GrapeGridDisplay<GridType>::deleteHmesh()
  {
    if( hmesh_ )
    {
      GrapeInterface<dim,dimworld>::deleteHmesh(hmesh_);
    }
  }

} // end namespace Dune
