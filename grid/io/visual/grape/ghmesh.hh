// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __GRAPE_HMESH_H__
#define __GRAPE_HMESH_H__

enum { MAX_NAME_LENGTH = 32 };

typedef struct dune_elem DUNE_ELEM;
typedef struct dune_fdata DUNE_FDATA;
typedef struct dune_dat DUNE_DAT;

typedef void evalDof_t (DUNE_ELEM *, DUNE_FDATA *, int , double *);
typedef void evalCoord_t (DUNE_ELEM *, DUNE_FDATA *, const double *, double * );

/* interface element */
struct dune_elem
{

  // default constructor
  dune_elem()
    : type(127)
      , eindex(-1)
      , level(-1)
      , level_of_interest(-1)
      , has_children(0)
      , liter(0)
      , enditer(0)
      , hiter(0)
      , actElement(0)
      , gridPart(0)
      , display(0)
      , mesh(0)
  {
    // default set all coordinates to zero
    for(int i=0; i<MAX_EL_DOF; ++i)
    {
      vindex [i] = -1;
      vpointer[i] = (double *) coordinates[i];
      for(int j=0; j<3; ++j)
      {
        vpointer[i][j] = 0.0;
      }
    }
    for(int i=0; i<MAX_EL_FACE; ++i)
    {
      bnd [i] = -1;
    }
  }

  /*
   *  see g_eldesc.h for ElementType
   */
  int type;

  double *        vpointer [MAX_EL_DOF];
  double coordinates [MAX_EL_DOF][3];
  int vindex [MAX_EL_DOF] ;
  int bnd [MAX_EL_FACE] ;
  int eindex;
  int level;
  int level_of_interest;
  int has_children;

  /* is the pointer to LevelIterator or to LeafIterator */
  void          * liter;
  void          * enditer;

  // pointer fo hierarchic iterator */
  void          * hiter;

  /* points to actual iterator to compare an get type */
  /* down cast to EntityPointer */
  void          * actElement;

  /* actual choosen gridPart */
  void          * gridPart;

  // pointer to my display class
  void          * display;

  // pointer to mesh
  void          * mesh;
};

struct dune_fdata
{
  static std::set<DUNE_FDATA*>& dataList ()
  {
    static std::set<DUNE_FDATA*> dList;
    return dList;
  }

  // default constructor
  dune_fdata()
    : mynum (-1)
      , name()
      , evalCoord(0)
      , evalDof(0)
      , discFunc(0)
      , indexSet(0)
      , allLevels(0)
      , dimVal(0)
      , dimRange(0)
      , comp(0)
      , polyOrd(0)
      , continuous(0)
      , compName(0)
      , gridPart(0)
      , setGridPartIterators(0)
      , f_data (0)
      , minValue(0.0)
      , maxValue(1.0)
      , valuesSet(false)
      , valCache(0.0)
      , getMinMaxValues(0)
  {
    // add this data to list of dune data funcs
    dataList().insert(this);
  }

  // default destructor
  ~dune_fdata()
  {
    dataList().erase(this);
  }

  /* my number in the data vector */
  int mynum;

  /* name of data */
  std::string name;

  // functions to evaluate
  evalCoord_t * evalCoord;
  evalDof_t   * evalDof;

  /* pointer to object of discrete function or vector */
  const void *discFunc;

  /* pointer to index set of underlying datas */
  const void *indexSet;

  /* are all Levels occupied? */
  int allLevels;

  /* dimension of value, i.e. the length of the vector  */
  int dimVal;

  /* dimension of data, when vectorial data is interpreted as scalar data */
  int dimRange;

  /* index of current component */
  /* for scalar this vec has length 1 and contains the component number */
  /* for vector this contains the number of each component */
  int * comp;

  /* polynonial order of basis functions */
  int polyOrd;

  /* continuous or not */
  int continuous;

  /* max number of components */
  int compName;

  /* the corresponding gridPart */
  void * gridPart;

  /* function pointer to choose grid part iterators */
  void (*setGridPartIterators)(DUNE_DAT * , void * gridPart);

  /* pointer to f_data */
  void * f_data;

  /* minValue of function, for colorbar */
  double minValue;
  /* maxValue of function, for colorbar */
  double maxValue;

  /* true if min and max values have been calculated */
  bool valuesSet;

  /* cache for polOrd zero functions */
  double valCache;

  /* returns min and max values of function */
  void (*getMinMaxValues)(DUNE_FDATA *, double * min, double * max );
};

/* dune_dat */
struct dune_dat
{
  // default constructor
  dune_dat()
    : first_macro(0)
      , next_macro(0)
      , delete_iter(0)
      , first_child(0)
      , next_child(0)
      , copy(0)
      , check_inside(0)
      , wtoc(0)
      , ctow(0)
      , setIterationModus(0)
      , partition(-1)
      , iteratorType(-1) // g_LeafIterator
      , partitionIteratorType(-1)
      , gridPart(0)
      , all (0)
      , get_stackentry(0)
      , free_stackentry(0) {}

  /* the actual first and next macro for Iteration  */
  int (* first_macro)(DUNE_ELEM *) ;
  int (* next_macro)(DUNE_ELEM *) ;

  /* method to delete iterators */
  void (* delete_iter)(DUNE_ELEM *) ;

  /* first and next child , if 0, then no child iteration */
  int (* first_child)(DUNE_ELEM *) ;
  int (* next_child)(DUNE_ELEM *) ;

  void * (* copy)(const void *) ;

  int (* check_inside)(DUNE_ELEM *, const double * ) ;
  int (* wtoc)(DUNE_ELEM *, const double *, double * ) ;
  void (* ctow)(DUNE_ELEM *, const double *, double * ) ;


  /* selects the iterators, like leaf iterator .. */
  void (* setIterationModus)(DUNE_DAT *, DUNE_FDATA *);

  /* to which processor partition the element belongs */
  int partition;

  /* type of choosen iterator */
  int iteratorType;

  /* type of partition to iterate */
  int partitionIteratorType;

  /* actual gridPart */
  void * gridPart;

  DUNE_ELEM * all;

  /* get HELEMENT */
  void * (*get_stackentry)(DUNE_DAT * );
  /* free HELEMENT */
  void (*free_stackentry)(DUNE_DAT * , void *);
};

/* setup hmesh with given data */
extern void *setupHmesh(const int noe, const int nov,
                        const int maxlev, DUNE_DAT * dune);

/* delete given hmesh pointer */
extern void deleteHmesh( void * hmesh );
extern void deleteFunctions( void * hmesh );

extern void displayTimeScene(INFO * info);
extern void handleMesh (void *hmesh, bool gridMode );

extern DUNE_FDATA * extractData (void *hmesh , int num );

/* setup TimeScene Tree  */
extern void timeSceneInit(INFO *info, const int n_info, const int procs);
extern void addDataToHmesh(void  *hmesh, DUNE_FDATA * data);

extern void addHmeshToTimeScene(void * timescene, double time, void  *hmesh , int proc);

extern void addHmeshToGlobalTimeScene(double time, void  *hmesh , int proc);
extern void tsc_timebar(void *timescene, double t_start, double t_end);
extern void colorBarMinMax(const double min, const double max);

#endif
