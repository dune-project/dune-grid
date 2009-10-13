// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{

  //****************************************************************
  //
  // --CombinedGrapeDisplay, CombinedGrapeDisplay for given grid
  //
  //****************************************************************

  template<class DisplayType>
  inline CombinedGrapeDisplay<DisplayType>::
  CombinedGrapeDisplay() : disp_(0) , hmesh_ (0)
  {
    GrapeInterface<dim,dimworld>::init();
  }

  template<class DisplayType>
  inline CombinedGrapeDisplay<DisplayType>::
  ~CombinedGrapeDisplay()
  {
    if( hmesh_ )
    {
      GrapeInterface<dim,dimworld>::deleteHmesh(hmesh_);
    }

    // still work to do
    /*
       for(size_t i=0 ;i<vecFdata_.size(); i++)
       {
       if( vecFdata_[i] ) DisplayType::deleteDuneFunc(vecFdata_[i]);
       vecFdata_[i] = 0;
       }
     */
    DisplayType::deleteStackEntry(stackEntry_);
  }


  //****************************************************************
  //
  // --GridDisplay, Some Subroutines needed in display
  //
  //****************************************************************
  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  first_macro (DUNE_ELEM * he)
  {
    grditer_ = dispList_.begin();
    enditer_ = dispList_.end();

    partEnd_  = gridPartList_.end();
    partIter_ = gridPartList_.begin();

    disp_ = 0;

    return callFirstMacro(he);
  }

  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  callFirstMacro(DUNE_ELEM * he)
  {
    if(grditer_ != enditer_)
    {
      disp_ = *grditer_;
      GrapeInterface<dim,dimworld>::setThread( disp_->myRank() );
      he->display = (void *) disp_;
      void * gridPart = he->gridPart;

      // set appropriate grid part
      if( partIter_ != partEnd_) he->gridPart = *partIter_;

      // call first macro of current display
      int ret = disp_->firstMacro(he);

      // set value from before
      he->gridPart = gridPart;
      he->display = (void *) this;
      return ret;
    }
    return 0;
  }

  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  next_macro (DUNE_ELEM * he)
  {
    if( disp_ )
    {
      void * gridPart = he->gridPart;
      he->display = (void *) disp_;

      // set appropriate grid part
      if( partIter_ != partEnd_) he->gridPart = *partIter_;

      int ret = disp_->nextMacro(he);

      if(!ret)
      {
        ++grditer_;
        if( partIter_ != partEnd_) ++partIter_;
        disp_ = 0;

        return callFirstMacro(he);
      }
      else
      {
        he->display  = (void *)this;
        he->gridPart = gridPart;
        return ret;
      }
    }
    else
    {
      return 0;
    }
  }

  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  first_child(DUNE_ELEM * he)
  {
    if(disp_)
    {
      he->display = (void *) disp_;
      int ret = disp_->firstChild(he);
      he->display = (void *)this;
      return ret;
    }
    else
      return 0;
  }


  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  next_child(DUNE_ELEM * he)
  {
    if( disp_ )
    {
      he->display = (void *) disp_;
      int ret = disp_->nextChild(he);
      he->display = (void *)this;
      return ret;
    }
    else
      return 0;
  }


  template<class DisplayType>
  inline void * CombinedGrapeDisplay<DisplayType>::
  copy_iterator (const void * i)
  {
    std::cerr << "ERROR: copt_iterator not implemented! file = " << __FILE__ << ", line = " << __LINE__ << "\n";
    abort () ;
    return 0 ;
  }

  // check inside
  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  checkInside(DUNE_ELEM * he, const double * w)
  {
    assert( disp_ );
    return disp_->checkWhetherInside(he,w);
  }
  // check inside
  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  check_inside(DUNE_ELEM * he, const double * w)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    return disp[0].checkInside(he,w);
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  local_to_world (DUNE_ELEM * he, const double * c, double * w)
  {
    assert( disp_ );
    disp_->local2world(he,c,w);
    return ;
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  ctow (DUNE_ELEM * he, const double * c, double * w)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    disp[0].local_to_world(he,c,w);
    return;
  }

  // world to local
  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  world_to_local(DUNE_ELEM * he, const double * w, double * c)
  {
    assert(disp_);
    return disp_->world2local(he,w,c);
  }

  // world to local
  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  wtoc(DUNE_ELEM * he, const double * w, double * c)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    return disp[0].world_to_local(he,w,c);
  }

  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  first_mac (DUNE_ELEM * he)
  {
    MyDisplayType & disp = *((MyDisplayType *) he->display);
    return disp.first_macro(he);
  }


  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  next_mac (DUNE_ELEM * he)
  {
    MyDisplayType & disp = *((MyDisplayType *) he->display);
    return disp.next_macro(he);
  }

  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  fst_child (DUNE_ELEM * he)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    return disp->first_child(he);
  }


  template<class DisplayType>
  inline int CombinedGrapeDisplay<DisplayType>::
  nxt_child (DUNE_ELEM * he)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    return disp->next_child(he);
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  setIterationModus(DUNE_DAT * dat, DUNE_FDATA * func)
  {
    MyDisplayType * disp = (MyDisplayType *) dat->all->display;
    disp->setIterationMethods(dat,func);
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  setIterationMethods(DUNE_DAT * dat, DUNE_FDATA * func)
  {
    enditer_ = dispList_.end();
    const int iteratorType = dat->iteratorType;
    const int partitionIteratorType = dat->partitionIteratorType;

    gridPartList_.clear();

    for(grditer_ = dispList_.begin(); grditer_ != enditer_; ++grditer_)
    {
      DisplayType & disp = * (*grditer_);

      DUNE_FDATA * data = 0;
      if(func)
      {
        std::vector < DUNE_FDATA * > & vec = disp.getFdataVec();
        data = vec[func->mynum];
        assert( data->gridPart );
        gridPartList_.push_back( data->gridPart );
      }

      disp.changeIterationMethods(iteratorType,partitionIteratorType,data);
    }

    assert( (func) ? (gridPartList_.size() == dispList_.size()) : true );
  }

  template<class DisplayType>
  inline void * CombinedGrapeDisplay<DisplayType>::
  getStackEn(DUNE_DAT * dat)
  {
    MyDisplayType * disp = (MyDisplayType *) dat->all->display;
    assert( disp );
    return DisplayType::getStackEntry(disp->stackEntry_);
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  freeStackEn(DUNE_DAT * dat, void * entry)
  {
    MyDisplayType * disp = (MyDisplayType *) dat->all->display;
    assert( disp );
    DisplayType::freeStackEntry(disp->stackEntry_,entry);
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::display()
  {
    /* call handle mesh in g_hmesh.c */
    GrapeInterface<dim,dimworld>::handleMesh ( hmesh_ );
    return ;
  }

  template<class DisplayType>
  inline void * CombinedGrapeDisplay<DisplayType>::getHmesh()
  {
    if(!hmesh_) hmesh_ = setupHmesh();
    return (void *) hmesh_;
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  evalCoord (DUNE_ELEM *he, DUNE_FDATA *df, const double *coord, double * val)
  {
    assert( disp_ );
    std::vector < DUNE_FDATA * > & vec = disp_->getFdataVec();
    DUNE_FDATA * data = vec[df->mynum];
    data->evalCoord(he,data,coord,val);
    return ;
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  evalDof (DUNE_ELEM *he, DUNE_FDATA *df,int localNum, double * val)
  {
    assert( disp_ );
    std::vector < DUNE_FDATA * > & vec = disp_->getFdataVec();
    DUNE_FDATA * data = vec[df->mynum];
    data->evalDof(he,data,localNum,val);
    return ;
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  evalCoordWrap (DUNE_ELEM *he, DUNE_FDATA *df, const double *coord, double * val)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    assert( disp );
    disp->evalCoord(he,df,coord,val);
    return ;
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  evalDofWrap (DUNE_ELEM *he, DUNE_FDATA *df,int localNum, double * val)
  {
    MyDisplayType * disp = (MyDisplayType *) he->display;
    assert( disp );
    disp->evalDof(he,df,localNum,val);
    return ;
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::addDisplay(DisplayType & disp)
  {
    dispList_.push_back( &disp );
    if(!hmesh_) hmesh_ = setupHmesh();

    if(disp.hasData())
    {
      // get functions data vector of given display
      std::vector < DUNE_FDATA * > & vec = disp.getFdataVec();

      // only copy functions for the first partition, because
      // all functions should be the same on every partition
      if(vec.size() > vecFdata_.size())
      {
        assert( vecFdata_.size() == 0 );
        // mem leak
        vecFdata_.clear();
        vecFdata_.resize( vec.size() );

        for(size_t n = 0; n < vecFdata_.size(); n++)
        {
          assert( n < vec.size());

          vecFdata_[n] = DisplayType::createDuneFunc();
          assert( vecFdata_[n] );

          // save pointer
          void * f_data = vecFdata_[n]->f_data;

          // copy, note that comp and f_data pointer are overwritten
          // we have to reset them
          std::memcpy(vecFdata_[n],vec[n],sizeof(DUNE_FDATA));

          DUNE_FDATA * data = vecFdata_[n];
          // set f_data pointer to our allocated value
          data->f_data = f_data;

          // we only need the two new functions for evaluation
          data->evalCoord = this->evalCoordWrap;
          data->evalDof = this->evalDofWrap;

          // not needed here
          data->comp = 0;

          // add function data to hmesh
          GrapeInterface<dim,dimworld>::addDataToHmesh(hmesh_,vecFdata_[n]);
        }
      }
    }
  }

  template<class DisplayType>
  inline void CombinedGrapeDisplay<DisplayType>::
  addMyMeshToGlobalTimeScene(double time, int proc)
  {
    if(!hmesh_) hmesh_ = setupHmesh();
    GrapeInterface<dim,dimworld>::addHmeshToGlobalTimeScene(time,this->getHmesh(),proc);
  }

  template<class DisplayType>
  inline void * CombinedGrapeDisplay<DisplayType>::setupHmesh()
  {
    int noe = 0, nov = 0;
    int maxlevel = 0;

    enditer_ = dispList_.end();
    for(grditer_ = dispList_.begin(); grditer_ != enditer_; ++grditer_)
    {
      const GridType & grid = (*grditer_)->getGrid();
      maxlevel = std::max( maxlevel, grid.maxLevel());
      noe += grid.size(0);
      nov += grid.size(dim);
    }

    // set display pointer
    hel_.display = (void *) this;

    // set dune data
    DUNE_DAT * dune = &dune_;

    dune->first_macro = &first_mac;
    dune->next_macro  = &next_mac;

    dune->first_child = &fst_child;
    dune->next_child  = &nxt_child;

    dune->wtoc         = wtoc;
    dune->ctow         = ctow;
    dune->check_inside = check_inside;

    // set method to select iterators
    dune->setIterationModus = &setIterationModus;

    dune->get_stackentry  = &getStackEn;
    dune->free_stackentry = &freeStackEn;

    dune->all          = &hel_;
    dune->partition    = __MaxPartition-1;

    dune->iteratorType          = g_LeafIterator;
    dune->partitionIteratorType = g_All_Partition;

    this->setIterationMethods(dune,0);

    /* return hmesh with no data */
    return GrapeInterface<dim,dimworld>::setupHmesh(noe,nov,maxlevel,dune);
  }

} // end namespace Dune
