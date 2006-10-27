// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRAPE_DATA_DISPLAY_HH
#define DUNE_GRAPE_DATA_DISPLAY_HH

//- system includes
#include <vector>

//- local includes
#include "grapegriddisplay.hh"

/** @file
   @author Robert Kloefkorn
   @brief Provides a DataDisplay class using the GridDisplay and
   dune-fem module for class DiscreteFubctions support.
 */

namespace Dune
{

  template <typename ctype, int dim, int dimworld, int polOrd>
  class GrapeLagrangePoints;

  template <class EvalImpTraits>
  struct EvalFunctionData
  {
    typedef typename EvalImpTraits :: GridType GridType;
    typedef typename EvalImpTraits :: EvalImp EvalImp;

    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;

    // for the data visualization, call implementations evalCoordNow
    inline static void evalCoordNow (EntityType &en, DUNE_FDATA *fdata, const double *coord, double * val)
    {
      EvalImp::evalCoordNow(en,fdata,coord,val);
    }

    // for the data visualization, call implementations evalDofNow
    inline static void evalDofNow (EntityType &en, int geomType, DUNE_FDATA *fdata , int localNum, double * val)
    {
      EvalImp::evalDofNow(en,geomType,fdata,localNum,val);
    }

    // evaluate at given local coord
    inline static void evalCoord (DUNE_ELEM *he, DUNE_FDATA *df,
                                  const double *coord, double * val);

    // evaluate at dof
    inline static void evalDof (DUNE_ELEM *he, DUNE_FDATA *df, int localNum, double * val);

    // get min and max value for colorbar
    inline static void getMinMaxValues(DUNE_FDATA *df, double* min, double* max );
  };

  template <class GridImp, class DiscreteFunctionType>
  struct EvalDiscreteFunctions;

  template <class GridImp, class DiscreteFunctionType>
  struct EvalDiscreteFunctionsTraits
  {
    typedef GridImp GridType;
    typedef EvalDiscreteFunctions <GridImp, DiscreteFunctionType > EvalImp;
  };

  template <class GridImp, class DiscreteFunctionType>
  struct EvalDiscreteFunctions
    : public EvalFunctionData< EvalDiscreteFunctionsTraits <GridImp, DiscreteFunctionType > >
  {
    typedef GridImp GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType :: FunctionSpaceType FunctionSpaceType;

    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: DomainType DomainType;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;

    // for the data visualization
    inline static void evalCoordNow (EntityType &en, DUNE_FDATA *fdata, const double *coord, double * val);

    // for the data visualization
    inline static void evalDofNow (EntityType &en, int geomType, DUNE_FDATA *fdata , int localNum, double * val);

    // for the data visualization
    inline static void evalScalar (EntityType &en, int geomType,
                                   DiscreteFunctionType & func, LocalFunctionType &lf,
                                   const int * comp , int localNum, double * val);

    // for the data visualization
    inline static void evalVector (EntityType &en, int geomType,
                                   DiscreteFunctionType & func, LocalFunctionType &lf,
                                   const int * comp, int vend, int localNum, double * val);

    // calculate min and max value of function
    inline static void calcMinMax(DUNE_FDATA * df);
  };

  template <class GridImp, class VectorType, class IndexSetImp >
  struct EvalVectorData;

  template <class GridImp, class VectorType , class IndexSetImp >
  struct EvalVectorDataTraits
  {
    typedef GridImp GridType;
    typedef EvalVectorData <GridImp, VectorType, IndexSetImp > EvalImp;
  };

  template <class GridImp, class VectorType, class IndexSetImp >
  struct EvalVectorData
    : public EvalFunctionData< EvalVectorDataTraits <GridImp, VectorType, IndexSetImp > >
  {
    typedef GridImp GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;

    // for the data visualization
    inline static void evalCoordNow (EntityType &en, DUNE_FDATA *fdata, const double *coord, double * val);

    // for the data visualization
    inline static void evalDofNow (EntityType &en, int geomType, DUNE_FDATA *fdata , int localNum, double * val);

    // for the data visualization, evaluate linear funcs
    inline static void evalVectorLinear (EntityType &en, int geomType,
                                         VectorType & func, const IndexSetImp & set,
                                         const int * comp, int vend, int localNum, double * val);

    // for the data visualization, evaluate const funcs
    inline static void evalVectorConst (EntityType &en, int geomType,
                                        VectorType & func, const IndexSetImp & set,
                                        const int * comp, int vend, int localNum, double * val);

    // calculate min and max value of function
    inline static void calcMinMax(DUNE_FDATA * df);
  };


  template<class GridType>
  class GrapeDataDisplay : public GrapeGridDisplay < GridType >
  {
    typedef GrapeDataDisplay < GridType > MyDisplayType;

    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_DAT DUNE_DAT;
    typedef typename GrapeInterface<dim,dimworld>::F_DATA F_DATA;

  public:
    typedef GridType MyGridType;

    //! Constructor, make a GrapeDataDisplay for given grid
    inline GrapeDataDisplay(const GridType &grid, const int myrank = -1);

    //! Constructor, make a GrapeDataDisplay for given grid
    template <class GridPartType>
    inline GrapeDataDisplay(const GridPartType & gridPart, const int myrank = -1);

    //! Desctructor
    inline ~GrapeDataDisplay();

    /*! display data stored in vector
       @param name Name of data (i.e. solution)
       @param data Data vector storing data to display
       @param indexSet The corresponding index set related to the data
       @param dinf GrapeDataDisplay internal data
       @param polOrd polynominal order of Lagrangespace, only 0 and 1 allowed
       at the momnent
       @ param continuous continuous or not (i.e polOrd = 0 ==> discontinuous) default is discontinuous
     */
    template<class VectorType, class IndexSetType >
    inline void displayVector(const std::string name,
                              const VectorType & data, const IndexSetType & indexSet,
                              const int polOrd , int dimRange, bool continuous = false);

    //! Calls the display of the grid and draws the discrete function
    //! if discretefunction is NULL, then only the grid is displayed
    template <class DiscFuncType>
    inline void dataDisplay(const DiscFuncType &func, bool vector = false);

    //! display grid and data without grid mode
    inline void display();

    //! add discrete function to display
    template <class DiscFuncType>
    inline void addData(const DiscFuncType &func, double time = 0.0);

    //! add discrete function to display
    template <class DiscFuncType>
    inline void addData(const DiscFuncType &func, const DATAINFO * , double time );

    //! add discrete function to display
    template <class DiscFuncType>
    inline void addData(const DiscFuncType &func, std::string name , double time , bool vector = false );

    // retrun whether we have data or not
    bool hasData () { return (vecFdata_.size() > 0); }

    // return vector for copying in combined display
    std::vector < DUNE_FDATA * > & getFdataVec () { return vecFdata_; }

  private:
    /*! add vector to display
       @param data Data vector storing data to display
       @param indexSet The corresponding index set related to the data
       @param dinf GrapeDataDisplay internal data
       @param time simulation time of data
       @param polOrd polynominal order of Lagrangespace, only 0 and 1 allowed
       at the momnent
       @ param continuous continuous or not (i.e polOrd = 0 ==> discontinuous)
     */
    template<class VectorType, class IndexSetType >
    inline void addVector(const std::string name,
                          const VectorType & data, const IndexSetType & indexSet,
                          const double time , const int polOrd ,
                          const int dimRange, bool continuous );

    /*! add vector to display
       @param data Data vector storing data to display
       @param indexSet The corresponding index set related to the data
       @param dinf GrapeDataDisplay internal data
       @param time simulation time of data
       @param polOrd polynominal order of Lagrangespace, only 0 and 1 allowed
       at the momnent
       @ param continuous continuous or not (i.e polOrd = 0 ==> discontinuous)
     */
    template<class VectorType, class IndexSetType >
    inline void addVector(const VectorType & data, const IndexSetType & indexSet,
                          const DATAINFO * dinf, double time ,
                          const int polOrd , const int dimRange, bool continuous );

    //! hold the diffrent datas on this mesh
    std::vector < DUNE_FDATA * > vecFdata_;

    enum { polynomialOrder = 1 };
    // store lagrange points for evaluation
    GrapeLagrangePoints<ctype,dim,dimworld,polynomialOrder> lagrangePoints_;

    typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type;
    typedef void evalCoord_t (EntityCodim0Type &, DUNE_FDATA *, const double *, double * );
    typedef void evalDof_t (EntityCodim0Type &,int , DUNE_FDATA * , int , double * );

  protected:
    template <class GridPartType>
    struct IterationMethodsGP
    {
      // wrapper methods for first_item and next_item
      inline static int fst_item (DUNE_ELEM * he)
      {
        assert( he->display );
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template first_item<GridPartType>(he);
      }
      inline static int nxt_item (DUNE_ELEM * he)
      {
        assert( he->display );
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        return disp.template next_item<GridPartType>(he);
      }

      // delete iterators
      inline static void del_iter (DUNE_ELEM * he)
      {
        assert( he->display );
        MyDisplayType & disp = *((MyDisplayType *) he->display);
        typedef typename GridPartType :: template Codim<0> :: IteratorType IteratorType;
        disp.template delete_iterators<IteratorType> (he);
      }
    };

    template <class GridPartImp>
    struct SetIter
    {
      static void setGPIterator (DUNE_DAT * dune ,void * gridPart)
      {
        assert( gridPart );
        dune->gridPart = gridPart;
        dune->first_macro = &IterationMethodsGP<GridPartImp>::fst_item;
        dune->next_macro  = &IterationMethodsGP<GridPartImp>::nxt_item;
        dune->delete_iter = &IterationMethodsGP<GridPartImp>::del_iter;

        dune->first_child = 0;
        dune->next_child = 0;
      }
    };

  public:
    // create object DUNE_FDATA
    static DUNE_FDATA * createDuneFunc ();
    // delete object DUNE_FDATA
    static void deleteDuneFunc (DUNE_FDATA *);
  };

  template <typename ctype, int dim, int dimworld, int polOrd>
  class GrapeLagrangePoints
  {
    enum { maxPoints = 20 };
    enum { numberOfTypes = (dim == 2) ? 2 : 6 };

    std::vector < FieldMatrix<ctype,maxPoints,dim> > points_;
  public:
    //! create lagrange points for given polyOrder and dim,dimworld
    GrapeLagrangePoints ()
    {
      for(int type=0; type<numberOfTypes; type++)
      {
        FieldMatrix<ctype,maxPoints,dim> coords(0.0);
        int nvx = numberOfVertices(type);

        for(int i=0; i<nvx; i++)
        {
          const double * p = getCoordinate(type,i);
          for(int j=0; j<dimworld; j++)
          {
            assert( p );
            coords[i][j] = p[j];
          }
        }
        points_.push_back( coords );
      }
    }

    //! return lagrange point with localNum
    //! for given element type and polyOrder
    const FieldVector<ctype,dim> &
    getPoint (int geomType, int polyOrder , int localNum ) const
    {
      assert( polOrd == polyOrder );
      assert( geomType >= 0 );
      assert( geomType < numberOfTypes );
      return points_[geomType][localNum];
    }

  private:
    int numberOfVertices( int type )
    {
      if(type < 2)
        return GrapeInterface_two_two::getElementDescription(type)->number_of_vertices;
      else
        return GrapeInterface_three_three::getElementDescription(type)->number_of_vertices;
    }

    const double * getCoordinate( int type, int i )
    {
      if(type < 2)
        return GrapeInterface_two_two::getElementDescription(type)->coord[i];
      else
        return GrapeInterface_three_three::getElementDescription(type)->coord[i];
    }
  };

} // end namespace Dune

#include "grape/grapedatadisplay.cc"
#endif
