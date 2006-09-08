// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "geldesc.hh"

namespace Dune
{

  //*******************************************************************
  //
  //  Routines for evaluation of the data
  //
  //*******************************************************************
  template<class EvalImp>
  inline void EvalFunctionData<EvalImp>::
  evalCoord (DUNE_ELEM *he, DUNE_FDATA *df, const double *coord, double * val)
  {
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    EntityPointerType * ep = (EntityPointerType *) he->actElement;
    assert( ep );
    evalCoordNow(*ep[0],df,coord,val);
  }

  template<class EvalImp>
  inline void EvalFunctionData<EvalImp>::
  evalDof (DUNE_ELEM *he, DUNE_FDATA *df,int localNum, double * val)
  {
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    EntityPointerType * ep = (EntityPointerType *) he->actElement;
    assert( ep );
    int geomType = he->type;
    evalDofNow( *ep[0] ,geomType,df,localNum,val);
    return ;
  }

  //*******************************************************************
  //  --EvalDiscreteFunctions
  //*******************************************************************
  // evaluate scalar functions, means val has length 1
  template <class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalScalar (EntityType &en, int geomType,
              DiscreteFunctionType & func, LocalFunctionType &lf,
              const int * comp, int localNum, double * val)
  {
    enum { polynomialOrder = FunctionSpaceType :: polynomialOrder };
    static const GrapeLagrangePoints<ctype,dim,dimworld,polynomialOrder> lagrangePoints;
    const FieldVector<ctype,dim> & localPoint =
      lagrangePoints.getPoint(geomType,polynomialOrder,localNum);

    static RangeType tmp_;
    // evaluate local function on local lagrange point
    lf.evaluateLocal(en,localPoint,tmp_);

    // dimval == 1 here
    // 0 because we only have one value (dimVal == 1)
    val[0] = tmp_[comp[0]];
    return;
  }

  template <class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalVector (EntityType &en, int geomType,
              DiscreteFunctionType & func, LocalFunctionType &lf,
              const int * comp, int vlength , int localNum, double * val)
  {
    enum { dim = EntityType::dimension };
    assert( comp );
    // dimval == dimension here
    // in that case we have only one dof that has to be returned for all
    // corners , kind of hack, but works for the moment
    if(func.getFunctionSpace().polynomOrder() == 0)
    {
      enum { polynomialOrder = FunctionSpaceType :: polynomialOrder };
      static const GrapeLagrangePoints<ctype,dim,dimworld,polynomialOrder> lagrangePoints;
      const FieldVector<ctype,dim> & localPoint =
        lagrangePoints.getPoint(geomType,polynomialOrder,localNum);

      static RangeType tmp_;
      // evaluate local function on local lagrange point
      lf.evaluateLocal(en,localPoint,tmp_);

      for(int i=0; i<vlength; i++)
      {
        val[i] = tmp_[comp[i]];
      }
      return;
    }
    else
    {
      std::cerr << "ERROR: evalVector for polOrd > 0 not implemented! file = " << __FILE__ << ", line = " << __LINE__ << "\n";
      abort();
    }
    return;
  }

  template <class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalDofNow (EntityType &en, int geomType, DUNE_FDATA *df , int localNum, double * val)
  {
    assert( df );
    assert( df->discFunc );

    DiscreteFunctionType & func = *((DiscreteFunctionType *)  ( df->discFunc));

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    enum { dim = EntityType::dimension };
    {
      const int * comp = df->comp;
      assert( comp );

      LocalFuncType lf = func.localFunction( en );

      int dimVal = df->dimVal;
      switch (dimVal)
      {
      case 1 :
      {
        evalScalar(en,geomType, func,lf,comp,localNum,val);
        return;
      }

      case dim :
      {
        evalVector(en,geomType,func,lf,df->comp,dimVal,localNum,val);
        return;
      }
      default :
      {
        assert(false);
        evalVector(en,geomType,func,lf,df->comp,dimVal,localNum,val);
        return;
      }
      }
      return;
    }
  }

  template<class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalCoordNow(EntityType &en, DUNE_FDATA *df , const double *coord, double * val)
  {
    assert( coord );
    enum { dim = GridType::dimension };

    assert( df );
    assert( df->discFunc );
    DiscreteFunctionType & func = *((DiscreteFunctionType *) ( df->discFunc));

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;


    const int * comp = df->comp;
    assert( comp );

    // get local function
    LocalFunctionType lf = func.localFunction( en );

    static DomainType domTmp_;

    // convert double to FieldVector
    for(int i=0; i<dim; i++) domTmp_[i] = coord[i];

    static RangeType tmp_;
    // evaluate local function on local (on reference element)
    // point == domTmp
    lf.evaluateLocal( en , domTmp_ , tmp_);

    const int dimVal = df->dimVal;
    for(int i=0; i<dimVal; ++i) val[i] = tmp_[comp[i]];
    return;
  }


  //*******************************************************************
  //  --EvalVectorData
  //*******************************************************************
  template <class GridType, class VectorType, class IndexSetImp >
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>::
  evalVectorLinear (EntityType &en, int geomType,
                    VectorType & func, const IndexSetImp & set,
                    const int * comp, int vlength , int localNum, double * val)
  {
    if(! set.contains(en) ) return ;

    int idx = vlength * set.template subIndex<dim> (en,localNum) ;
    for(int i=0; i<vlength; ++i) val[i] = func[idx + comp[i]];
    return ;
  }

  template <class GridType, class VectorType, class IndexSetImp >
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>::
  evalVectorConst (EntityType &en, int geomType,
                   VectorType & func, const IndexSetImp & set,
                   const int * comp, int vlength , int localNum, double * val)
  {
    if(! set.contains(en) ) return ;

    int idx = vlength * set.index(en) ;
    val[0] = func[idx + comp[0]];
    return ;
  }

  template <class GridType, class VectorType, class IndexSetImp >
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>::
  evalDofNow (EntityType &en, int geomType, DUNE_FDATA *df, int localNum, double * val)
  {
    assert( df );
    assert( df->discFunc );

    VectorType & func = *((VectorType *)  (df->discFunc));

    assert( df->indexSet );
    const IndexSetImp * set = (const IndexSetImp *) df->indexSet;

    const int * comp = df->comp;
    assert( comp );
    int dimRange = df->dimRange;

    int polOrd = df->polyOrd;

    if( polOrd > 0 )
      evalVectorLinear(en,geomType,func,*set,comp,dimRange,localNum,val);
    else
      evalVectorConst(en,geomType,func,*set,comp,dimRange,localNum,val);
    return ;
  }

  template <class GridType, class VectorType,class IndexSetImp>
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>::
  evalCoordNow(EntityType &en, DUNE_FDATA *df , const double *coord, double * val)
  {
    assert( false );
    abort();
  }

  //****************************************************************
  //
  // --GrapeDataDisplay, GrapeDataDisplay for given grid
  //
  //****************************************************************
  template <class GridType>
  inline GrapeDataDisplay<GridType>::
  GrapeDataDisplay (const GridType &grid, const int myrank ) :
    GrapeGridDisplay < GridType > (grid,myrank) , vecFdata_ (0)
  {}

  template <class GridType>
  template <class GridPartType>
  inline GrapeDataDisplay<GridType>::
  GrapeDataDisplay (const GridPartType &gridPart, const int myrank ) :
    GrapeGridDisplay < GridType > (gridPart,myrank) , vecFdata_ (0)
  {}

  template <class GridType>
  inline GrapeDataDisplay<GridType>::~GrapeDataDisplay()
  {
    GrapeInterface<dim,dimworld>::deleteFunctions(this->hmesh_);

    for(size_t i=0 ; i<vecFdata_.size(); i++)
    {
      if( vecFdata_[i] ) deleteDuneFunc(vecFdata_[i]);
      vecFdata_[i] = 0;
    }
  }

  template<class GridType>
  inline void GrapeDataDisplay<GridType>::
  deleteDuneFunc(DUNE_FDATA * fd)
  {
    if( fd )
    {
      int * comps = fd->comp;
      if( comps ) delete [] comps;

      F_DATA * f_data = (F_DATA *) fd->f_data;
      if( f_data )
      {
        char * name = f_data->name;
        if(name)
        {
          std::string tmp(name);
          // use free here, because mem has be allocated with malloc
          // when name has been overwritten, dont free memory
          if(tmp == fd->name)
          {
            std::free(name);
            f_data->name = 0;
          }
        }

        delete f_data;
        fd->f_data = 0;
      }
      delete fd;
    }
  }

  template<class GridType>
  inline typename GrapeDataDisplay<GridType>::DUNE_FDATA * GrapeDataDisplay<GridType>::
  createDuneFunc()
  {
    DUNE_FDATA * func = new DUNE_FDATA();
    assert (func );
    F_DATA * f_data = new F_DATA ();
    f_data->name = 0;

    func->f_data = (void *) f_data;
    return func;
  }


  //****************************************************************
  //
  // --GrapeDataDisplay, Some Subroutines needed in display
  //
  //****************************************************************
  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  dataDisplay(DiscFuncType &func, bool vector)
  {
    /* add function data */
    this->addData(func,func.name(),0.0,vector);

    /* display mesh */
    GrapeInterface<dim,dimworld>::handleMesh ( this->hmesh_ );
    return ;
  }

  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  addData(DiscFuncType &func, double time)
  {
    typedef typename DiscFuncType::FunctionSpaceType FunctionSpaceType;
    enum { dimR = FunctionSpaceType::DimRange };
    int comp[dimR];
    for(int i=0; i<dimR; ++i) comp[i] = i;
    std::string name = func.name();
    // name, base_name, next, dimVal, comp
    DATAINFO dinf = { name.c_str() , name.c_str() , 0 , dimR , (int *) &comp };
    addData(func,&dinf,time);
  }

  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  addData(DiscFuncType &func , std::string name , double time , bool vector)
  {
    int comp[dim];
    for(int i=0; i<dim; i++) comp[i] = i;
    DATAINFO dinf = { name.c_str() , name.c_str() , 0 , (vector) ? dim : 1 , (int *) &comp };
    addData(func,&dinf,time);
  }

  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  addData(DiscFuncType &func , const DATAINFO * dinf, double time )
  {
    typedef typename DiscFuncType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscFuncType::LocalFunctionType LocalFuncType;

    assert(dinf);
    std::string name(dinf->name);
    assert( dinf->dimVal > 0);
    bool vector = (dinf->dimVal > 1) ? true : false;

    bool already=false;
    int size = vecFdata_.size();

    // add function wether is exists or not
    if(!already)
    {
      int num = (int) FunctionSpaceType::DimRange;
      if(vector) num = 1;

      vecFdata_.resize(size+num);
      for(int n=size; n < size+num; n++)
      {
        vecFdata_[n] = createDuneFunc();
        // set data components
        {
          DUNE_FDATA * data = vecFdata_[n];
          assert( data );

          // set the rigth evaluation functions
          data->evalDof =
            EvalDiscreteFunctions<GridType,DiscFuncType>::evalDof;

          data->evalCoord =
            EvalDiscreteFunctions<GridType,DiscFuncType>::evalCoord;

          data->mynum = n;

          data->allLevels = 0;

          data->discFunc = (void *) &func;
          data->indexSet = 0;
          data->polyOrd = func.getFunctionSpace().polynomOrder();
          data->continuous = (func.getFunctionSpace().continuous() == true ) ? 1 : 0;
          if(data->polyOrd == 0) data->continuous = 0;

          int dimVal = dinf->dimVal;
          int * comp = new int [dimVal];
          data->comp = comp;
          if(vector)
          {
            for(int j=0; j<dimVal; j++) comp[j] = dinf->comp[j];
            data->compName = -1;
          }
          else
          {
            comp[0] = n-size;
            data->compName = n-size;
          }

          if(data->compName >= 0)
          {
            std::stringstream str;
            str << name << "[" << data->compName << "]";
            data->name = str.str();
          }
          else
            data->name = name;

          data->dimVal   = dimVal;
          data->dimRange = FunctionSpaceType::DimRange;

          // set grid part selection methods
          typedef typename FunctionSpaceType :: GridPartType GridPartType;
          data->gridPart = ((void *) &func.getFunctionSpace().gridPart());
          data->setGridPartIterators = &SetIter<GridPartType>::setGPIterator;
        }

        GrapeInterface<dim,dimworld>::addDataToHmesh(this->hmesh_,vecFdata_[n]);
      }
    }
  }

  template<class GridType>
  template<class VectorType,class IndexSetType>
  inline void GrapeDataDisplay<GridType>::
  displayVector(const std::string name,
                const VectorType & data,
                const IndexSetType & indexSet,
                const int polOrd,
                const int dimRange,
                bool continuous )
  {
    // polOrd < 0 makes no sense
    assert( polOrd >= 0 );

    // only polord 0 or 1 supported
    assert( polOrd < 2 );

    // add to display
    this->addVector(name,data,indexSet,0.0,polOrd,dimRange,continuous);

    double min = data[0];
    double max = data[0];

    int cd = (polOrd == 0) ? 0 : dim ;

    for(int i=0; i<indexSet.size(cd)*dimRange; ++i)
    {
      min = std::min(data[i],min);
      max = std::max(data[i],max);
    }
    setMinMaxValue(min,max);

    // display
    GrapeInterface<dim,dimworld>::handleMesh ( this->hmesh_ );
    return;
  }

  template<class GridType>
  template<class VectorType,class IndexSetType>
  inline void GrapeDataDisplay<GridType>::
  addVector(const std::string name,
            const VectorType & data,
            const IndexSetType & indexSet,
            const double time,
            const int polOrd,
            const int dimRange,
            bool continuous )
  {
    int * comp = new int [dimRange];
    for(int i=0; i<dimRange; ++i) comp[i] = i;

    DATAINFO dinf = { name.c_str() , name.c_str() , 0 , 1 , comp };
    /* add function data */
    this->addVector(data,indexSet,&dinf,time,polOrd,dimRange,continuous);

    delete [] comp;
    return;
  }

  template<class GridType>
  template<class VectorType, class IndexSetType >
  inline void GrapeDataDisplay<GridType>::
  addVector(const VectorType & func , const IndexSetType & indexSet,
            const DATAINFO * dinf, const double time ,
            const int polOrd , const int dimRange, bool continuous )
  {
    assert(dinf);
    const char * name = dinf->name;
    assert( dinf->dimVal > 0);

    bool vector = (dinf->dimVal > 1) ? true : false;

    bool already=false;
    int size = vecFdata_.size();

    // add function wether is exists or not
    if(!already)
    {
      int num = dimRange;
      if(vector) num = 1;

      vecFdata_.resize(size+num);
      for(int n=size; n < size+num; n++)
      {
        vecFdata_[n] = createDuneFunc();
        {
          DUNE_FDATA * data = vecFdata_[n]->all;
          assert(data);

          // set the rigth evaluation functions
          data->evalDof =
            EvalVectorData<GridType,VectorType,IndexSetType>::evalDof;

          data->evalCoord =
            EvalVectorData<GridType,VectorType,IndexSetType>::evalCoord;

          data->mynum = n;
          data->nameValue = name;
          data->name = data->nameValue.c_str();
          data->allLevels = 0;

          data->discFunc = (void *) &func;
          data->indexSet = (void *) &indexSet;
          data->polyOrd  = polOrd;
          data->continuous = (continuous == true ) ? 1 : 0;
          if(data->polyOrd == 0) data->continuous = 0;

          int dimVal = dinf->dimVal;
          int * comp = new int [dimVal];
          data->comp = comp;
          if(vector)
          {
            for(int j=0; j<dimVal; j++) comp[j] = dinf->comp[j];
            data->compName = -1;
          }
          else
          {
            comp[0] = n-size;
            data->compName = n-size;
          }
          data->dimVal   = dimVal;
          data->dimRange = dimRange;

          data->gridPart = 0;
          data->setGridPartIterators = 0;
        }

        GrapeInterface<dim,dimworld>::addDataToHmesh(this->hmesh_,vecFdata_[n]);
      }
    }

    getFdataVec();
  }

  template<class GridType>
  inline void GrapeDataDisplay<GridType>::
  setMinMaxValue(const double minValue, const double maxValue) const
  {
    GrapeInterface<dim,dimworld>::colorBarMinMax(minValue,maxValue);
  }

} // end namespace Dune
