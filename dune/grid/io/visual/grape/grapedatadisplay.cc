// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if HAVE_GRAPE
#include "geldesc.hh"
#endif
#include <dune/geometry/referenceelements.hh>

namespace Dune
{

#if HAVE_GRAPE
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

  template<class EvalImp>
  inline void EvalFunctionData<EvalImp>::
  getMinMaxValues(DUNE_FDATA *df, double* min, double *max)
  {
    if(!df->valuesSet)
    {
      EvalImp::calcMinMax(df);
    }
    *min = df->minValue;
    *max = df->maxValue;
    return ;
  }

  //*******************************************************************
  //  --EvalDiscreteFunctions
  //*******************************************************************
  // evaluate scalar functions, means val has length 1
  template <class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalScalar (const EntityType &en, int geomType,
              DiscreteFunctionType & func, LocalFunctionType &lf,
              const int * comp, int localNum, double * val)
  {
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };
    static const GrapeLagrangePoints<ctype,dim,dimworld,polynomialOrder> lagrangePoints;
    const FieldVector<ctype,dim> & localPoint =
      lagrangePoints.getPoint(geomType,polynomialOrder,localNum);

    RangeType tmp_;
    // evaluate local function on local lagrange point
    lf.evaluate(localPoint,tmp_);

    // dimval == 1 here
    // 0 because we only have one value (dimVal == 1)
    val[0] = tmp_[comp[0]];
    return;
  }

  template <class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalVector (const EntityType &en, int geomType,
              DiscreteFunctionType & func, LocalFunctionType &lf,
              const int * comp, int vlength , int localNum, double * val)
  {
    enum { dim = EntityType::dimension };
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType SpaceType;
    enum { polOrd = SpaceType :: polynomialOrder };
    assert( comp );
    // dimval == dimension here
    // in that case we have only one dof that has to be returned for all
    // corners , kind of hack, but works for the moment
    if(polOrd == 0)
    {
      enum { polynomialOrder = SpaceType :: polynomialOrder };
      static const GrapeLagrangePoints<ctype,dim,dimworld,polynomialOrder> lagrangePoints;
      const FieldVector<ctype,dim> & localPoint =
        lagrangePoints.getPoint(geomType,polynomialOrder,localNum);

      RangeType tmp_;
      // evaluate local function on local lagrange point
      lf.evaluate(localPoint,tmp_);

      for(int i=0; i<vlength; ++i)
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
  evalDofNow (const EntityType &en, int geomType, DUNE_FDATA *df , int localNum, double * val)
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
      if (dim==1 || dimVal==1)
      {
        evalScalar(en,geomType, func,lf,comp,localNum,val);
      }
      else if (dim!=1 && dimVal==dim)
      {
        evalVector(en,geomType,func,lf,df->comp,dimVal,localNum,val);
      }
      else
      {
        assert(false);
        evalVector(en,geomType,func,lf,df->comp,dimVal,localNum,val);
      }
      return;
    }
  }

  template<class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  evalCoordNow ( const EntityType &entity, DUNE_FDATA *df, const double *coord, double * val )
  {
    const int dim = GridType::dimension;

    assert( df );
    assert( coord );

    const DiscreteFunctionType *function = static_cast< const DiscreteFunctionType * >( df->discFunc );
    assert( function );

    const int *comp = df->comp;
    assert( comp );

    // get local function
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunction;
    const LocalFunction localFunction = function->localFunction( entity );

    // convert double to FieldVector
    typename EntityType::Geometry::LocalCoordinate x;
    for( int i = 0; i < dim; ++i )
      x[ i ] = coord[ i ];

    // evaluate local function in local (on reference element) point x
    RangeType tmp;
    localFunction.evaluate( x, tmp );

    const int dimVal = df->dimVal;
    for( int i = 0; i < dimVal; ++i )
      val[ i ] = tmp[ comp[ i ] ];
  }


  template<class GridType, class DiscreteFunctionType>
  inline void EvalDiscreteFunctions<GridType,DiscreteFunctionType>::
  calcMinMax(DUNE_FDATA * df)
  {
    double minValue = 0.0;
    double maxValue = 1.0;

    assert( df->discFunc );
    const DiscreteFunctionType & func = *((const DiscreteFunctionType *) (df->discFunc));

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType ;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dimension = GridType :: dimension };

    if(df->dimVal == 1)
    {
      const DiscreteFunctionSpaceType & space = func.space();
      bool initialized = false;
      const IteratorType end = space.end();
      for(IteratorType it = space.begin(); it != end; ++it)
      {
        const EntityType & en = *it;
        int geomType = convertToGrapeType ( en.type() , dimension );
        double val = 0.0;
        for(int i=0; i<en.template count<dimension>(); ++i)
        {
          evalDofNow ( en , geomType, df , i , &val );
          if(!initialized)
          {
            minValue = maxValue = val;
            initialized = true;
          }
          minValue = std::min(minValue,val);
          maxValue = std::max(maxValue,val);
        }
      }
    }
    else
    {
      derr << "EvalDiscreteFunctions::calcMinMax: method not implemented for vectorial data! \n";
    }

    if((maxValue-minValue) < 1e-10)
    {
      std::cout << "WARNING: min="<<minValue << " and max="<<maxValue<<" values of function almost identical! \n";
      maxValue += 0.01*maxValue;
      minValue -= 0.01*minValue;
    }

    df->minValue = minValue;
    df->maxValue = maxValue;
    df->valuesSet= true;
  }



  // EvalGrapeFunction
  // -----------------

  template< class GV, int dimR, int polOrd >
  inline void EvalGrapeFunction< GV, dimR, polOrd >
  ::evalCoordNow ( const Entity &entity, DUNE_FDATA *fdata, const double *coord, double *val )
  {
    assert( (fdata != 0) && (coord != 0) && (val != 0) );

    const GrapeFunction *function = (const GrapeFunction *)(fdata->discFunc);
    assert( function != 0 );

    const DomainVector &x = reinterpret_cast< const DomainVector & >( *coord );
    RangeVector &y = reinterpret_cast< RangeVector & >( *val );

    function->evaluate( entity, x, y );
  }


  template< class GV, int dimR, int polOrd >
  inline void EvalGrapeFunction< GV, dimR, polOrd >
  ::evalDofNow ( const Entity &entity, int geomType, DUNE_FDATA *fdata, int localNum, double *val )
  {
    assert( (fdata != 0) && (val != 0) );

    const GrapeFunction *function = (const GrapeFunction *)(fdata->discFunc);
    assert( function != 0 );

    static const GrapeLagrangePoints< typename GridView::Grid::ctype, dimDomain, dimWorld, polOrd > lagrangePoints;
    const DomainVector &x = lagrangePoints.getPoint( geomType, polOrd, localNum );

    RangeVector &y = reinterpret_cast< RangeVector & >( *val );
    function->evaluate( entity, x, y );
  }


  template< class GV, int dimR, int polOrd >
  inline void EvalGrapeFunction< GV, dimR, polOrd >
  ::calcMinMax ( DUNE_FDATA *fdata )
  {
    assert( (fdata != NULL) && (fdata->discFunc != NULL) );

    double minValue = std::numeric_limits< double >::infinity();
    double maxValue = -std::numeric_limits< double >::infinity();

    const GrapeFunction *function = (const GrapeFunction *)(fdata->discFunc);
    const GridView &gridView = function->gridView();

    if( dimR == 1 )
    {
      typedef typename GridView::template Codim< 0 >::Iterator Iterator;
      typedef ReferenceElement< typename GridView::Grid::ctype, dimDomain > ReferenceElement;
      typedef ReferenceElements< typename GridView::Grid::ctype, dimDomain > ReferenceElements;

      const Iterator end = gridView.template end< 0 >();
      for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
      {
        const Entity &entity = *it;
        const ReferenceElement &refElement = ReferenceElements::general( entity.type() );

        const int numCorners = refElement.size( dimDomain );
        for( int i = 0; i < numCorners; ++i )
        {
          RangeVector y;
          function->evaluate( entity, refElement.position( i, dimDomain ), y );

          minValue = std::min( minValue, y[ 0 ] );
          maxValue = std::max( maxValue, y[ 0 ] );
        }
      }
    }
    else
    {
      std::cerr << "EvalGrapeFunction::calcMinMax not implemented for dimR > 1." << std::endl;
      minValue = 0.0;
      maxValue = 1.0;
    }

    if( (maxValue - minValue) < 1e-10 )
    {
      std::cout << "Warning: min ("<< minValue << ") and max (" << maxValue << ") values are almost identical." << std::endl;
      maxValue += 0.01 * maxValue;
      minValue -= 0.01 * minValue;
    }

    fdata->minValue = minValue;
    fdata->maxValue = maxValue;
    fdata->valuesSet= true;
  }



  //*******************************************************************
  //  --EvalVectorData
  //*******************************************************************
  template <class GridType, class VectorType, class IndexSetImp >
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>
  ::evalVectorLinear ( const EntityType &entity, int geomType,
                       VectorType &func, const IndexSetImp &indexSet,
                       const int *comp, int vlength, int localNum, double *val )
  {
    if( indexSet.contains( entity ) )
    {
      //int idx = vlength * indexSet.template subIndex<dim> (entity,localNum);
      int idx = vlength * indexSet.subIndex( entity, localNum, dim );
      val[ 0 ] = func[ idx + comp[ 0 ] ];
    }
  }

  template <class GridType, class VectorType, class IndexSetImp >
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>
  ::evalVectorConst ( const EntityType &entity, int geomType,
                      VectorType &func, const IndexSetImp &indexSet,
                      const int *comp, int vlength, int localNum, double *val )
  {
    if( indexSet.contains( entity ) )
    {
      int idx = vlength * indexSet.index( entity );
      val[ 0 ] = func[ idx + comp[ 0 ] ];
    }
  }

  template <class GridType, class VectorType, class IndexSetImp >
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>::
  evalDofNow (const EntityType &en, int geomType, DUNE_FDATA *df, int localNum, double * val)
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
  evalCoordNow(const EntityType &en, DUNE_FDATA *df , const double *coord, double * val)
  {
    assert( false );
    abort();
  }

  template <class GridType, class VectorType,class IndexSetImp>
  inline void EvalVectorData<GridType,VectorType,IndexSetImp>::
  calcMinMax(DUNE_FDATA * df)
  {
    double minValue,maxValue;
    assert( df->discFunc );
    VectorType & vector = *((VectorType *) (df->discFunc));

    IndexSetImp * set = ((IndexSetImp *) df->indexSet);
    assert( set );

    int size = (df->polyOrd > 0) ? set->size(dim) : set->size(0);

    int comp = df->comp[0];
    int dimVal = df->dimVal;

    // get first value of vector to set min and max
    if(size < comp) return;
    minValue = vector[comp];
    maxValue = vector[comp];

    for(int i=0; i<size; ++i)
    {
      if((i%dimVal) != comp) continue;
      if( vector[i] < minValue) minValue = vector[i];
      if( vector[i] > maxValue) maxValue = vector[i];
    }

    if((maxValue-minValue) < 1e-10)
    {
      std::cout << "WARNING: min="<<minValue << " and max="<<maxValue<<" values of function almost identical! \n";
      maxValue += 0.01*maxValue;
      minValue -= 0.01*minValue;
    }

    df->minValue = minValue;
    df->maxValue = maxValue;
    df->valuesSet= true;
  }
#endif

  //****************************************************************
  //
  // --GrapeDataDisplay, GrapeDataDisplay for given grid
  //
  //****************************************************************
  template <class GridType>
  inline GrapeDataDisplay<GridType>::
  GrapeDataDisplay (const GridType &grid, const int myrank ) :
    GrapeGridDisplay < GridType > (grid,myrank)
#if HAVE_GRAPE
    , vecFdata_ (0)
#endif
  {}

  template <class GridType>
  template <class GridPartType>
  inline GrapeDataDisplay<GridType>::
  GrapeDataDisplay (const GridPartType &gridPart, const int myrank ) :
    GrapeGridDisplay < GridType > (gridPart,myrank)
#if HAVE_GRAPE
    , vecFdata_ (0)
#endif
  {}

  template <class GridType>
  inline GrapeDataDisplay<GridType>::~GrapeDataDisplay()
  {
#if HAVE_GRAPE
    GrapeInterface<dim,dimworld>::deleteFunctions(this->hmesh_);

    for(size_t i=0 ; i<vecFdata_.size(); i++)
    {
      if( vecFdata_[i] ) deleteDuneFunc(vecFdata_[i]);
      vecFdata_[i] = 0;
    }
#endif
  }

#if HAVE_GRAPE
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
#endif


  //****************************************************************
  //
  // --GrapeDataDisplay, Some Subroutines needed in display
  //
  //****************************************************************
  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  dataDisplay(const DiscFuncType &func, bool vector)
  {
#if HAVE_GRAPE
    /* add function data */
    this->addData(func,func.name(),0.0,vector);

    /* display mesh */
    GrapeInterface<dim,dimworld>::handleMesh ( this->hmesh_ );
#endif
    return ;
  }

  template<class GridType>
  inline void GrapeDataDisplay<GridType>::
  display()
  {
#if HAVE_GRAPE
    /* display mesh without grid mode */
    GrapeInterface<dim,dimworld>::handleMesh ( this->hmesh_ );
#endif
    return ;
  }

  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  addData(const DiscFuncType &func, double time, bool vector )
  {
    this->addData(func,func.name(),time,vector);
  }

  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  addData(const DiscFuncType &func , std::string name , double time , bool vector)
  {
#if HAVE_GRAPE
    int comp[dim];
    for(int i=0; i<dim; i++) comp[i] = i;
    DATAINFO dinf = { name.c_str() , name.c_str() , 0 , (vector) ? dim : 1 , (int *) &comp };
    addData(func,&dinf,time);
#endif
  }

#if HAVE_GRAPE
  template<class GridType>
  template<class DiscFuncType>
  inline void GrapeDataDisplay<GridType>::
  addData(const DiscFuncType &func , const DATAINFO * dinf, double time )
  {
    typedef typename DiscFuncType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };
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
      int num = (int) DiscreteFunctionSpaceType::dimRange;
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

          data->getMinMaxValues =
            EvalDiscreteFunctions<GridType,DiscFuncType>::getMinMaxValues;

          data->mynum = n;

          data->allLevels = 0;

          data->discFunc = (void *) &func;
          data->indexSet = 0;
          data->polyOrd = (int) polynomialOrder;
          data->continuous = (func.space().continuous() == true ) ? 1 : 0;
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
          data->dimRange = DiscreteFunctionSpaceType::dimRange;

          // set grid part selection methods
          typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
          data->gridPart = ((void *) &func.space().gridPart());
          data->setGridPartIterators =
            &BaseType::template SetIter<GridPartType>::setGPIterator;
        }

        GrapeInterface<dim,dimworld>::addDataToHmesh(this->hmesh_,vecFdata_[n]);
      }

      // make grid part iterator default
      GrapeInterface<dim,dimworld>::setDefaultIterator(g_GridPart);
    }
  }
#endif

  template< class GridType >
  template< class GV, int dimR, int polOrd >
  inline void GrapeDataDisplay< GridType >
  ::addData ( const GrapeFunction< GV, dimR, polOrd > &function )
  {
#if HAVE_GRAPE
    DUNE_FDATA *fdata = createDuneFunc();
    assert( fdata != 0 );

    const unsigned int n = vecFdata_.size();
    vecFdata_.push_back( fdata );

    fdata->evalDof = EvalGrapeFunction< GV, dimR, polOrd >::evalDof;
    fdata->evalCoord = EvalGrapeFunction< GV, dimR, polOrd >::evalCoord;
    fdata->getMinMaxValues = EvalGrapeFunction< GV, dimR, polOrd >::getMinMaxValues;

    fdata->mynum = n;
    fdata->allLevels = 0;

    fdata->discFunc = (void *)&function;
    fdata->indexSet = 0;
    fdata->polyOrd = polOrd;
    fdata->continuous = false;

    const int dimRange = GrapeFunction< GV, dimR, polOrd >::dimRange;
    fdata->dimVal = dimRange;
    fdata->dimRange = dimRange;
    fdata->comp = new int[ dimRange ];
    for( int j = 0; j < dimRange; ++j )
      fdata->comp[ j ] = j;
    fdata->compName = -1;
    fdata->name = function.name();

    fdata->gridPart = (void*)&(function.gridView());
    fdata->setGridPartIterators
      = &BaseType::template GridViewIterators< typename GV::Traits >::set;

    GrapeInterface< dim, dimworld >::addDataToHmesh( this->hmesh_, vecFdata_[ n ] );

    // make grid view iterator default
    GrapeInterface< dim, dimworld >::setDefaultIterator( g_GridPart );
#endif
  }



  template<class GridType>
  template<class VectorType,class IndexSetType>
  inline void GrapeDataDisplay<GridType>::
  displayVector(const std::string name,
                const VectorType & data,
                const IndexSetType & indexSet,
                const int polOrd,
                const unsigned int dimRange,
                bool continuous )
  {
#if HAVE_GRAPE
    typedef typename IndexSetType::IndexType IndexType;

    // polOrd < 0 makes no sense
    assert( polOrd >= 0 );

    // only polord 0 or 1 supported
    assert( polOrd < 2 );

    // add to display
    this->addVector(name,data,indexSet,0.0,polOrd,dimRange,continuous);

    double min = data[0];
    double max = data[0];

    const int codim = (polOrd == 0) ? 0 : dim;
    const IndexType dataSize = indexSet.size( codim ) * dimRange;
    for( IndexType i = 0; i < dataSize; ++i )
    {
      min = std::min(data[i],min);
      max = std::max(data[i],max);
    }
    if(std::abs(max-min) < 1e-10)
    {
      max += 0.01*max;
      min -= 0.01*min;
    }
    // GrapeInterface<dim,dim>::colorBarMinMax(min,max);

    // display
    GrapeInterface<dim,dimworld>::handleMesh ( this->hmesh_ );
#endif
    return;
  }

#if HAVE_GRAPE
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
    DATAINFO dinf = { name.c_str() , name.c_str() , 0 , 1 , 0 };
    /* add function data */
    this->addVector(data,indexSet,&dinf,time,polOrd,dimRange,continuous);
    return;
  }
#endif

#if HAVE_GRAPE
  template<class GridType>
  template<class VectorType, class IndexSetType >
  inline void GrapeDataDisplay<GridType>::
  addVector(const VectorType & func , const IndexSetType & indexSet,
            const DATAINFO * dinf, const double time ,
            const int polOrd , const int dimRange, bool continuous )
  {
    assert(dinf);
    const char * name = dinf->name;

    // only add data as scalar data
    assert( dinf->dimVal == 1);

    bool already=false;
    int size = vecFdata_.size();

    // add function wether is exists or not
    if(!already)
    {
      int num = dimRange;

      vecFdata_.resize(size+num);
      for(int n=size; n < size+num; n++)
      {
        vecFdata_[n] = createDuneFunc();
        {
          DUNE_FDATA * data = vecFdata_[n];
          assert(data);

          // set the rigth evaluation functions
          data->evalDof =
            EvalVectorData<GridType,VectorType,IndexSetType>::evalDof;

          data->evalCoord =
            EvalVectorData<GridType,VectorType,IndexSetType>::evalCoord;

          data->getMinMaxValues =
            EvalVectorData<GridType,VectorType,IndexSetType>::getMinMaxValues;

          data->mynum = n;

          data->allLevels = 0;

          data->discFunc = (void *) &func;
          data->indexSet = (void *) &indexSet;
          data->polyOrd  = polOrd;
          data->continuous = (continuous == true ) ? 1 : 0;
          if(data->polyOrd == 0) data->continuous = 0;

          int dimVal = dinf->dimVal;
          assert( dimVal == 1 );
          int * comp = new int [dimVal];
          data->comp = comp;
          comp[0] = n-size;
          data->compName = n-size;

          if(data->compName >= 0)
          {
            std::stringstream str;
            str << name << "[" << data->compName << "]";
            data->name = str.str();
          }
          else
            data->name = name;

          data->dimVal   = dimVal;
          data->dimRange = dimRange;

          data->gridPart = 0;
          data->setGridPartIterators = 0;
        }

        GrapeInterface<dim,dimworld>::addDataToHmesh(this->hmesh_,vecFdata_[n]);
      }
    }
  }
#endif

} // end namespace Dune
