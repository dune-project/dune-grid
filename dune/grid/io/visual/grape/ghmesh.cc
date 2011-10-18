// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __GRAPE_HMESH_C__
#define __GRAPE_HMESH_C__

#include "geldesc.hh"

/*****************************************************************************
* Globale defines                  *                     **
*****************************************************************************/

static const double EPSILON  = 1.0e-2;
#define STA(el) ((STACKENTRY *) (el))
static const double INFTY = 999999.;
#define FVAL_ON_NODE(b) (jump_data[b])

#define MINIMUM(a,b) (((a) > (b)) ? (b) : (a))
#define MAXIMUM(a,b) (((a) < (b)) ? (b) : (a))

GENMESHnD *genmesh3d_switch_part_light_model_on_off();
GENMESHnD *genmesh3d_switch_part_displaybar_on_off();
GENMESHnD *genmesh3d_partition_disp();

/*****************************************************************************
* Verwendete Strukturen                *                     **
*****************************************************************************/
typedef struct stackentry
{
  HELEMENT hel;

  int ref_flag ;

  double hmax;

  // default constructor
  stackentry () : hel()
                  , ref_flag(0), hmax(-1.0) {}
} STACKENTRY;

/* definition of dune_dat in g_eldesc.h */
/* stored as user_data in the mesh pointer */

namespace FctSelector {

  /* function list management */
#define FCT_NONE "none"
#define FCT_ZERO "zero"
#define SLT_DEFAULT "default"

#define SYMBOL_LENGTH 8

  typedef enum {
    fmtNone,
    fmtLog,
    fmtT_Fun
  } function_modifier_index;

  typedef struct function_modifier {
    function_modifier_index type;
    char name[MAX_STRLEN];
    GENMESH_FDATA *orig_f_data, *new_f_data;
    int new_f_size;
    void (*origfun)(void *, int, double *, double *, void *);
    void (*origbounds)(void *, double *, double *, void *);
    int number_of_functions;
    struct {
      GENMESH_FDATA *fdata;
      char name[MAX_STRLEN], symbol[SYMBOL_LENGTH];
      void (*fun)(void *, int, double *, double *, void *);
      int active;
    } *functions;
    ITEM *item;
    void *data;
    void (*modfun)(void *, int, double *, double *, void *);
    void (*modbounds)(void *, double *, double *, void *);
    GENMESH_FDATA *(*init)(struct function_modifier*);
    char *(*info)(struct function_modifier*);
    bool_t (*xdr)(XDR *, struct function_modifier*);
    void (*deleteFct)(struct function_modifier*);
  } function_modifier;

  /* function list management */
  typedef struct {
    char *slot, *name;
    function_modifier modifier;
  } chosen_function;

  typedef ::std::pair< const char * , const char * > DataFunctionName_t;

  // returns slot and name of selected function
  inline DataFunctionName_t getCurrentFunctionName(GRAPEMESH * mesh)
  {
    DataFunctionName_t fctName( (const char *)0, (const char *)0 );

    chosen_function * current = (chosen_function *)
                                g_list_current (mesh->current_function);
    if (current)
    {
      GENMESH_FDATA* fun = (GENMESH_FDATA *)
                           GRAPE_CALL(mesh, "get-function") (NULL, current->slot,NULL);
      fctName.first  = current->slot;
      fctName.second = fun->name;
    }
    return fctName;
  }

#undef FCT_NONE
#undef FCT_ZERO
#undef SLT_DEFAULT
#undef SYMBOL_LENGTH

} // end namespace fct selector

/**************************************************************************
*
* froward declarations
*
* ************************************************************************/

/* switch from leaf to level iteration for mesh */
int switchMethods( GENMESHnD *actHmesh);

/*****************************************************************************
******************************************************************************
**                      **
**  Memory Management Routinen              **
**                      **
******************************************************************************
*****************************************************************************/
inline static void gFreeElement(ELEMENT *el)
{
  assert( el );
  assert( el->mesh );
  DUNE_DAT * dat = (DUNE_DAT *) el->mesh->user_data;
  assert(dat);
  dat->free_stackentry(dat,(void *)el);
  return ;
}

/* complete element does nothing at the moment */
inline static ELEMENT * complete_element(ELEMENT *el, MESH_ELEMENT_FLAGS flags)
{
  return el;
}

/*****************************************************************************
* little help routines
*****************************************************************************/
inline static double dist(const double *x,  const double *y)
{
  double dist=0.0;
  int i;

  for(i=0; i<GRAPE_DIMWORLD; i++)
  {
    dist += (x[i]-y[i])*(x[i]-y[i]);
  }
  return (sqrt(dist));
}

inline static double calc_hmax(HELEMENT *el)
{
  return ( dist(el->vertex[0], el->vertex[1]) );
}

/*****************************************************************************
******************************************************************************
**
**  The functions to which each HMesh has a pointer, i.e. first_macro
**  here first_macro, next_macro, first_child, next_child
**
******************************************************************************
*****************************************************************************/
// update helement pointers
inline static void helementUpdate( DUNE_ELEM *elem, HELEMENT *grapeEl )
{
  // set pointers
  grapeEl->vertex       = elem->vpointer;
  grapeEl->vindex       = elem->vindex ;
  grapeEl->eindex       = elem->eindex ;
  grapeEl->level        = elem->level;
  grapeEl->has_children = elem->has_children;
  grapeEl->user_data    = (void *)elem ;

  // select appropriate element description
  grapeEl->descr = (ELEMENT_DESCRIPTION *) getElementDescription(elem->type);
  return ;
}

inline static HELEMENT * first_macro (GENMESHnD *mesh, MESH_ELEMENT_FLAGS flag)
{
  assert(mesh);
  DUNE_DAT * dat = (DUNE_DAT *) mesh->user_data;
  assert( dat );

  HELEMENT * el = (HELEMENT *) dat->get_stackentry(dat);
  DUNE_ELEM * elem = (DUNE_ELEM *) el->user_data;

  assert(el);
  assert(elem);

  // set iterator type depending on value of button
  GRAPE_CALL(iteratorButton,      "get-value") (&(dat->iteratorType));
  GRAPE_CALL(partitionTypeButton, "get-value") (&(dat->partitionIteratorType));

  F_DATA * f_data = (F_DATA*)GRAPE_CALL(mesh,"get-function")
                      ("scalar","scalar","vector","default", NULL);

  assert( dat->setIterationModus );
  DUNE_FDATA * func = 0;
  if(f_data)
  {
    DUNE_FDATA * data = (DUNE_FDATA *) f_data->function_data;
    if( DUNE_FDATA::dataList().find(data) !=
        DUNE_FDATA::dataList().end() )
    {
      func = data;
    }
  }
  dat->setIterationModus(dat,func);

  // note this pointer can be NULL
  elem->gridPart = dat->gridPart;

  /* store level of interest for LeafIterator */
  if(maxlevelButton->on_off == OFF) /* dont know why wrong, but it works */
    elem->level_of_interest = -1;
  else
    elem->level_of_interest = mesh->level_of_interest;

  el->present = hefAll;
  el->parent = NULL;

  elem->display = dat->all->display;
  elem->hiter = NULL;

  {
    /* call first macro and check for first element */
    int test = dat->first_macro(elem) ;
    // means no element exits at all
    if(!test) return NULL;
  }

  el->level = 0;
  el->mesh  = (GENMESHnD *)mesh ;

  helementUpdate(elem,el);
  ((STACKENTRY *)el)->hmax = calc_hmax(el);

  el->vinh    = NULL ;
  ((STACKENTRY *)el)->ref_flag   = -1;

  /***************************************************************/
  // is this assertion is thrown then something with the geometry types is
  // wrong
  assert( el->descr != 0 );

  return(el);
}

/* go next macro element */
inline static HELEMENT * next_macro(HELEMENT * el, MESH_ELEMENT_FLAGS flag)
{
  int mflag=0;
  assert(el) ;

  el->present = (MESH_ELEMENT_FLAGS) (hefAll & ! hefVinh);
  mflag = (*(((struct dune_dat *)(el->mesh->user_data))->next_macro))((DUNE_ELEM *)el->user_data);
  if(mflag)
  {
    helementUpdate( ((DUNE_ELEM *)el->user_data) , el);
    ((STACKENTRY *)el)->hmax = calc_hmax(el);
    return(el) ;
  }
  else
  {
    /*printf("next macro: bin draussen flag = %i \n",mflag);*/
    gFreeElement((ELEMENT *)el) ;
    return NULL ;
  }
}

/***********************************************************/
/* first_child, go to first child of current element */
/************************************************************/
inline static HELEMENT * first_child (HELEMENT * ael, MESH_ELEMENT_FLAGS flag)
{
  int actlevel = ael->level;

  DUNE_DAT * dat = (DUNE_DAT *) ael->mesh->user_data;
  assert( dat );

  // if pointer is zero then no child iteration
  if(dat->first_child)
  {
    if ( actlevel < ((HMESH *)ael->mesh)->level_of_interest )
    {
      HELEMENT * el = (HELEMENT *) dat->get_stackentry(dat);
      assert(el);
      // set mesh, needed for removal
      el->mesh = ael->mesh;

      DUNE_ELEM * elem = (DUNE_ELEM *)el->user_data;
      assert(elem);

      el->present = (MESH_ELEMENT_FLAGS) (hefAll & !hefVinh);

      DUNE_ELEM * aelem = (DUNE_ELEM *)ael->user_data;

      elem->display = aelem->display;
      elem->liter   = aelem->liter;
      elem->hiter   = aelem->hiter;

      /* call the dune method */
      if(dat->first_child(elem))
      {
        el->level = actlevel+1;
        el->mesh  = ael->mesh ;

        helementUpdate(elem,el);

        el->parent    = ael;
        ((STACKENTRY *)el)->hmax = ((STACKENTRY *)ael)->hmax * 0.5;

        el->vinh  = NULL ;
        ((STACKENTRY *)el)->ref_flag   = -1;
        /****************************************************/
        // is this assertion is thrown then something with the geometry types is
        // wrong
        assert( el->descr != 0 );
        return(el);
      }
      else
      {
        gFreeElement((ELEMENT *)el) ;
        return NULL ;
      }
    }
  }
  return NULL ;
}

/* go to next child of the current element */
inline static HELEMENT * next_child(HELEMENT * el, MESH_ELEMENT_FLAGS flag)
{
  assert(el) ;
  el->present = (MESH_ELEMENT_FLAGS) (hefAll & !hefVinh);
  DUNE_ELEM * elem = ((DUNE_ELEM *)el->user_data);
  assert( elem );

  DUNE_DAT * dat = (DUNE_DAT *) el->mesh->user_data;
  assert( dat );

  if(dat->next_child)
  {
    if(dat->next_child(elem))
    {
      DUNE_ELEM * elem = ((DUNE_ELEM *)el->user_data);
      ((STACKENTRY *)el)->ref_flag++;
      helementUpdate(elem,el);

      // only implemented on triangles
      if(elem->type != gr_triangle)
      {
        el->vinh = NULL;
      }

      return(el) ;
    }
    else
    {
      gFreeElement((ELEMENT *)el) ;
      return NULL ;
    }
  }

  return NULL;
}

/* first_child, go to first child of current element */
inline static HELEMENT * select_child (HELEMENT * ael, double *parent_coord,
                                       double *child_coord, MESH_ELEMENT_FLAGS flag)
{
  HELEMENT *child = NULL;
  HMESH *mesh = (HMESH *) ael->mesh;
  const ELEMENT_DESCRIPTION *descr = ael->descr;
  double coord [GRAPE_DIMWORLD];
  int inside = 0;

  /* in child_ccord the local coordinates of the point in child can be
   * stored */
  descr->coord_to_world(ael,parent_coord,(double *) &coord);

  if( mesh->first_child )
    child = mesh->first_child(ael, flag);
  if( !child ) return NULL;

  /* tranform to local coords of child */
  descr->world_to_coord(child, (double *) &coord, child_coord);

  while ( inside != -1 )
  {
    inside = descr->check_inside( child , child_coord );
    if( inside == -1)
    {
      return child;
    }

    child = mesh->next_child( child , flag);
    if ( ! child ) return NULL;
    descr->world_to_coord(child, (double *) &coord, child_coord);
  }
  return NULL;
}

inline static ELEMENT * first_element (GRAPEMESH *mesh, MESH_ELEMENT_FLAGS flag)
{
  return first_macro((GENMESHnD * )mesh,flag);
}

/* go next macro element */
inline static ELEMENT * next_element(ELEMENT * el, MESH_ELEMENT_FLAGS flag)
{
  return next_macro(el,flag);
}
/***************************************************************************
*
*  f_data function
*
***************************************************************************/
inline void f_bounds(HELEMENT *el, double* min, double* max,
                     void *function_data)
{
  (*min) =  1.0E+308;
  (*max) = -1.0E+308;
  return;
}
/****************************************************************************/
inline void grape_get_vertex_estimate(HELEMENT *el, double *value,
                                      void *function_data)
{
  *value = 1.0E+308;
  return;
}

/****************************************************************************/

inline double grape_get_element_estimate(HELEMENT *el, void *function_data)
{
  return 1.0E+308;
}

/***************************************************************************/
inline static void dune_function_info(HELEMENT *el, F_EL_INFO *f_el_info,
                                      void *function_data)
{
  assert( function_data );

  DUNE_FDATA * df = (DUNE_FDATA *) function_data;
  assert( df );

  // at the moment f_el_info only contains polynomial_degree
  f_el_info->polynomial_degree = df->polyOrd;
  return;
}

inline static void level_function_info(HELEMENT *el, F_EL_INFO *f_el_info,
                                       void *function_data)
{
  // the level function has polynomial order 0
  f_el_info->polynomial_degree = 0;
  return;
}

/***************************************************************************/

/* print DUNE_FDATA STRUCT */
inline void printfFdata(DUNE_FDATA *fem)
{
  printf("Dune Fdata %p \n",fem);
  printf("comp %d      | DiscFunc   %p \n",fem->comp[0],fem->discFunc);
  printf("-------------------------------------------\n");
}
inline void printDuneFunc(DUNE_FDATA *df)
{
  printf("DUNE_FDATA %p \n",df);
  printf("discFunc %p \n",df->discFunc);
  //printf("lf       %p \n",df->lf);
  printf("comp     %d \n",df->comp[0]);
}

/* call the function to evaluate, depending on index of local vertex or
 * local coordinate */
static inline void evaluateFunction(DUNE_ELEM * elem, DUNE_FDATA * fem, int ind,
                                    double G_CONST *coord, double * val)
{
  assert(elem);
  assert(fem);

  if(coord)
  {
    fem->evalCoord(elem,fem,coord,val);
  }
  else
  {
    fem->evalDof(elem,fem,ind,val);
  }
}

/* default function calling function to evaluate */
static inline void f_real(HELEMENT *el, int ind, double G_CONST *coord,
                          double *val, void *function_data)
{
  assert(el);
  DUNE_ELEM * elem = (DUNE_ELEM *)el->user_data;
  assert(elem != NULL);
  DUNE_FDATA *fem = (DUNE_FDATA *) function_data;

  assert(fem != NULL);
  assert(fem->discFunc != NULL);

  evaluateFunction(elem,fem,ind,coord,val);
  return;
}

/* function calling data to evaluate for picewise constant data */
static inline void f_real_polOrd_zero(HELEMENT *el, int ind, double G_CONST *coord,
                                      double *val, void *function_data)
{
  assert(el);
  DUNE_ELEM * elem = (DUNE_ELEM *)el->user_data;
  assert(elem != NULL);
  DUNE_FDATA *fem = (DUNE_FDATA *) function_data;

  assert(fem != NULL);
  assert(fem->discFunc != NULL);

  // if polynomial order zero
  // then only evaluate for ind = 0
  if(ind == 0)
  {
    evaluateFunction(elem,fem,ind,coord,val);
    // cache function value
    fem->valCache = val[0];
  }
  else
  {
    // use cached value
    val[0] = fem->valCache;
  }
  return;
}

/* function displaying the level of the element */
inline void f_level(HELEMENT *el, int ind, double G_CONST *coord,
                    double *val, void *function_data)
{
  assert(el);
  val[0] = (double) el->level;
  return;
}

/***************************************************************************/
inline void grapeInitScalarData(GRAPEMESH *grape_mesh, DUNE_FDATA * dfunc)
{
  assert( grape_mesh );

  if(dfunc)
  {
    F_DATA * f_data = (F_DATA *) dfunc->f_data;
    if (f_data)
    {
      f_data->next = 0;
      f_data->last = 0;

      int length = 2*strlen(dfunc->name.c_str());
      f_data->name = (char *) malloc(length*sizeof(char));
      sprintf(f_data->name,"%s",dfunc->name.c_str());
      printf("generate data for discrete function '%s'!\n",f_data->name);

      f_data->dimension_of_value = dfunc->dimVal;

      f_data->continuous_data = dfunc->continuous;

      f_data->f                   = (dfunc->polyOrd == 0) ? f_real_polOrd_zero : f_real;
      f_data->f_el_info           = dune_function_info;

      f_data->function_data = (void *) dfunc;

      f_data->get_bounds      = f_bounds;
      f_data->get_vertex_estimate   = grape_get_vertex_estimate;
      f_data->get_element_estimate  = grape_get_element_estimate;
      f_data->threshold     = 0.0;
#if GRAPE_DIM == 3
      f_data->geometry_threshold     = 0.0;
#else
      // if pointer 0, nothing done with this functions
      f_data->get_element_p_estimates = 0;
      f_data->get_edge_p_estimates    = 0;
#endif
      f_data->hp_threshold    = 0.0;
      f_data->hp_maxlevel     = grape_mesh->max_level;

      grape_mesh = (GRAPEMESH *) GRAPE_CALL(grape_mesh,"add-function") (f_data);
    }
    else if (grape_mesh->f_data != (GENMESH_FDATA *)f_data)
    {
      printf("select f_data for \n");
      grape_mesh->f_data = (GENMESH_FDATA *)f_data;
    }
  }
  else
  {
    printf("no dfunc, or no vec\n");
  }

  return;
}

/* generates the function to display the level of an element */
inline void grapeAddLevelFunction(GRAPEMESH *grape_mesh)
{
  assert( grape_mesh );

  /* function info for level display */
  /* the variables are only needed once, therefore static */

  static std::string level_name("level");
  static F_DATA levelData;

  static F_DATA * f_data = 0;

  if (!f_data)
  {
    f_data = &levelData;
    assert( f_data );

    f_data->name = (char *) level_name.c_str();
    f_data->dimension_of_value = 1;
    f_data->continuous_data    = 0;

    f_data->f                   = f_level;
    f_data->f_el_info           = &level_function_info;

    // no function data here
    f_data->function_data = 0;

    f_data->get_bounds      = f_bounds;
    f_data->get_vertex_estimate   = grape_get_vertex_estimate;
    f_data->get_element_estimate  = grape_get_element_estimate;
    f_data->threshold     = 0.0;
#if GRAPE_DIM == 3
    f_data->geometry_threshold     = 0.0;
#else
    f_data->get_element_p_estimates = 0;
    f_data->get_edge_p_estimates    = 0;
#endif
    f_data->hp_threshold    = 0.0;
    f_data->hp_maxlevel     = grape_mesh->max_level;
  }

  assert( f_data );
  grape_mesh = (GRAPEMESH *) GRAPE_CALL(grape_mesh,"add-function") (f_data);
  return;
}

/* add data to Hmesh  */
inline void addDataToHmesh(void  *hmesh, DUNE_FDATA * data)
{
  GRAPEMESH *mesh = (GRAPEMESH *) hmesh;
  assert(mesh);

  if(data)
  {
    assert( data->f_data );
    assert( data->discFunc );

    /* setup dune data */
    grapeInitScalarData (mesh, data );
  }
  else
  {
    fprintf(stderr,"ERROR: no function data for setup in addDataToHmesh! \n");
  }
}

/*****************************************************************************
******************************************************************************
**                      **
**  Die MESH Routinen "copy_element"                              **
**                      **
******************************************************************************
*****************************************************************************/

inline static ELEMENT * copy_element(ELEMENT *el, MESH_ELEMENT_FLAGS flag)
{
  DUNE_DAT * dat = (DUNE_DAT *) el->mesh->user_data;
  HELEMENT * cel = (HELEMENT *) dat->get_stackentry(dat);

  assert(el) ;
  assert(cel) ;

  DUNE_ELEM * hexa_elem = 0, * chexa_elem = 0;

  hexa_elem  = (DUNE_ELEM *) el->user_data;
  chexa_elem = (DUNE_ELEM *) dat->copy(hexa_elem) ;

  assert(chexa_elem) ;
  cel->mesh              = el->mesh ;
  cel->vertex            = chexa_elem->vpointer;
  cel->vindex            = el->vindex ;
  cel->eindex            = el->eindex ;
  cel->descr             = el->descr ;
  cel->parent             = ((HELEMENT *)el)->parent ;
  cel->user_data         = (void *)chexa_elem ;
  cel->level = el->level;

  return ( (ELEMENT *)cel );
}

inline static void get_geometry_vertex_estimate(HELEMENT* helement, double* results)
{
  /* planar mesh -> all geometry-estimates 0*/
  int i;
  for(i=0; i<3; i++)
    results[i] = 1e5;
  return;
}


inline static double get_geometry_element_estimate(HELEMENT* helement)
{
  /*planar mesh -> estimators 0*/
  return(1e5);
}

/* method to get partition number from mesh */
inline HMESH * get_partition_number (int * partition)
{
  HMESH * hmesh = (HMESH *) START_METHOD (G_INSTANCE);
  assert(hmesh != 0);
  DUNE_DAT * dunedata = (DUNE_DAT *) hmesh->user_data;
  assert(dunedata != 0);
  *partition = dunedata->partition;

  END_METHOD(hmesh);
}

inline GRAPEMESH *getMesh ( const char *meshName = "Dune Mesh" )
{
  return (GRAPEMESH *)GRAPE_CALL( GrapeMesh, "new-instance" ) ( meshName );
}

// if parameter is not 0 ,then mesh is freed (i.e. pushed to stack)
inline GRAPEMESH * getAndFreeMesh( GRAPEMESH * mesh = 0 )
{
  static std::stack< GRAPEMESH * > meshStack;

  if(mesh)
  {
    // free mesh
    meshStack.push(mesh);
    return 0;
  }
  else
  {
    // get mesh
    //if(meshStack.empty())
    {
      return getMesh();
    }
    /*
       // not working yet
       else
       {
       GRAPEMESH * m = meshStack.top();
       meshStack.pop();
       return m;
       }
     */
  }
}

/*****************************************************************************
******************************************************************************
**                      **
**  Die Hautroutine zum Initialisieren und Aufrufen eines HMESH"      **
**                      **
**  --setupHmesh
**  --hmesh
**
******************************************************************************
*****************************************************************************/
inline void * setupHmesh(const int noe, const int nov,
                         const int maxlev, DUNE_DAT * dune,
                         const char *meshName )
{
  GRAPEMESH *mesh = getMesh( meshName );

  assert(mesh != NULL);

  mesh->first_macro = first_macro ;
  mesh->next_macro  = next_macro ;

  mesh->first_child = first_child ;
  mesh->next_child  = next_child ;

  mesh->select_child  = select_child;

  mesh->copy_element  = copy_element ;
  mesh->free_element  = gFreeElement ;
  mesh->complete_element = complete_element;

  mesh->max_number_of_vertices = MAX_EL_DOF ;
  mesh->max_eindex = noe ;
  mesh->max_vindex = nov ;

#if GRAPE_DIM==2
  mesh->dimension_of_world = GRAPE_DIMWORLD;
#endif
  mesh->max_dimension_of_coord = GRAPE_DIMWORLD;
  mesh->max_dindex = 20;

  mesh->max_level = maxlev;
  mesh->level_of_interest = maxlev;

  mesh->get_geometry_vertex_estimate  = get_geometry_vertex_estimate;
  mesh->get_geometry_element_estimate = get_geometry_element_estimate;
  mesh->get_lens_element_estimate = 0;

  mesh->threshold              = 1.0;

  mesh->user_data = (void *) dune;

  mesh->set_time = 0;
  mesh->get_time = 0;
  mesh->f_data   = 0;

  grapeAddLevelFunction(mesh);

  return ((void *) mesh);
}

// delete Hmesh , not really working yet
inline void deleteHmesh( void * hmesh )
{
  getAndFreeMesh( (GRAPEMESH *) hmesh);
}

// delete Hmesh , not really working yet
inline void deleteFunctions( void * hmesh )
{
  //assert( hmesh );
  //GRAPEMESH * mesh = (GRAPEMESH *) hmesh;

  /*
     GENMESH_FDATA * f_data = mesh->f_data;
     while( f_data )
     {
     char * name = f_data->name;
     if(name)
     {
      std::cout << "Got " << name << " func \n";
      mesh = (GRAPEMESH *) GRAPE_CALL(mesh,"remove-function")(name);
     }
     f_data = mesh->f_data;
     }
   */

  //mesh->f_data = 0;
  //mesh->user_data = 0;
}

static inline void addProjectUIF()
{
  static int firstCall = 1;

  // only call this once otherwise
  // grape cannot be runed twice with the same program
  if(firstCall)
  {
    char p_name[32];
    sprintf(p_name,"uif-m%d",GRAPE_DIM);
    g_project_add(p_name);
    firstCall = 0;
  }
}

extern "C" {
  extern MESH2D   * mesh2d_isoline_disp();
  extern MESH2D   * mesh2d_isoline_select_disp();
  extern HPMESH2D * hpmesh2d_isoline_disp();
  extern GENMESH2D* genmesh2d_isoline_disp();
  extern GENMESH2D* genmesh2d_isoline_select_disp();
  extern GENMESH2D* genmesh2d_geometry_graph_disp();
  extern MESH2D *mesh2d_isoline_disp();
  extern MESH2D *mesh2d_isoline_select_disp();
  extern MESH2D *mesh2d_vect_disp();
  extern HPMESH2D* hpmesh2d_pdegfine_disp();
  extern GENMESH2D* genmesh2d_chess_disp();
  extern MESH2D *mesh2d_flic_disp();

  extern GENMESH3D* genmesh3d_bnd_isoline_disp();
  extern GENMESH3D* genmesh3d_bnd_isoline_select_disp();
  extern GENMESH3D* genmesh3d_clip_isoline_multi_disp ();
  extern GENMESH3D* genmesh3d_volume_disp ();
  extern MESH3D* mesh3d_bnd_isoline_select_disp();
  extern MESH3D* mesh3d_bnd_isoline_disp();
}

// methodName is set by scene set_min_max_value
static std::string grapeMethodName;

static inline void setMinMaxValuesToColorbars(const char* meshName,
                                              const double min,
                                              const double max)
{
  typedef std::list< COLORBAR * > ColorBarListType;
  static bool firstCall = true;
  static ColorBarListType colorBarList;

  if( firstCall )
  {
    ///////////////////////////////
    // all functions from mesh2d
    ///////////////////////////////
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh2d_isoline_disp,"isoline-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh2d_isoline_select_disp,"isoline-select-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (hpmesh2d_isoline_disp,"isoline-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh2d_geometry_graph_disp,"geometry-graph-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (mesh2d_isoline_select_disp,"mesh2d-isoline-select");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (mesh2d_isoline_disp,"mesh2d-isoline");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (mesh2d_vect_disp,"mesh2d-vect");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (hpmesh2d_pdegfine_disp,"pdegfine-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh2d_chess_disp, "GenMesh2d::chess-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    //{
    //  COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar")(mesh2d_flic_disp,"flic");
    //  if ( colorBar ) colorBarList.push_back( colorBar );
    //}
    ///////////////////////////////
    // all functions from mesh3d
    ///////////////////////////////
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh3d_bnd_isoline_disp,"bnd-isoline-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh3d_bnd_isoline_select_disp,"bnd-isoline-select-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh3d_clip_isoline_multi_disp,"genmesh3d-clip-isoline-multi-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh3d_volume_disp, "GenMesh3d-volume");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (mesh3d_bnd_isoline_select_disp,"bnd-isoline-select-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (mesh3d_bnd_isoline_disp,"bnd-isoline-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    {
      COLORBAR * colorBar = (COLORBAR *) GRAPE_CALL(Colorbar,"get-stdcolorbar") (genmesh3d_bnd_isoline_disp,"bnd-isoline-disp");
      if ( colorBar ) colorBarList.push_back( colorBar );
    }
    //////////////////////////////
    //////////////////////////////
    firstCall = false;
  }

#ifdef GRAPE_GENMESH_WITH_DIFFERENT_COLORBARS
  if ( meshName )
  {
    COLORBAR * cb = get_colorbar_from_handle( meshName );
    if( cb )
    {
      cb->min = min;
      cb->max = max;
    }
  }
#endif

  typedef ColorBarListType :: iterator iterator;

  const iterator end = colorBarList.end();
  for(iterator it = colorBarList.begin(); it != end; ++it)
  {
    COLORBAR * cb = (*it);
    cb->min = min;
    cb->max = max;
  }
}

inline GRAPEMESH * setMinMaxValue()
{
  GRAPEMESH * mesh = 0;
  mesh = (GRAPEMESH *) START_METHOD(G_INSTANCE);
  ASSURE (mesh, "No HMESH in setMinMaxValue! \n", END_METHOD (0));

  F_DATA * f_data = (F_DATA*)GRAPE_CALL(mesh,"get-function")
                      ("scalar","scalar","vector","default", 0);
  DUNE_FDATA * func = 0;
  if(f_data)
  {
    DUNE_FDATA * data = (DUNE_FDATA *) f_data->function_data;
    if( DUNE_FDATA::dataList().find(data) !=
        DUNE_FDATA::dataList().end() )
    {
      func = data;
    }
  }

  double min=0.0, max=1.0;
  if(func)
  {
    assert(func->getMinMaxValues);
    func->getMinMaxValues(func,&min,&max);
  }

  setMinMaxValuesToColorbars(mesh->name,min,max);

  END_METHOD(mesh);
}

inline void colorBarMinMax(const double min, const double max)
{
  GRAPE_CALL(Colorbar,"set-default-min-max") (min,max);
}

/* forward declaration */
static void grape_add_remove_methods(void);

inline void handleMesh(void *hmesh, bool gridMode )
{
  GRAPEMESH *mesh = (GRAPEMESH *) hmesh;
  assert(mesh != NULL);

  // remember last selected function
  static std::string lastFunctionName("zero");
  static std::string lastSlotName("default");

  // static stack to keep scenes
  static std::stack< SCENE *> sceneStack;

  MANAGER * mgr = (MANAGER *)GRAPE_CALL(Manager,"get-stdmgr") ();

  SCENE * sc = 0 ;
  if(sceneStack.empty())
  {
    sc = (SCENE *)GRAPE_CALL(Scene,"new-instance") ("dune hmesh");
  }
  else
  {
    sc = sceneStack.top();
    sceneStack.pop();
  }

  addProjectUIF();

#ifdef GRID_MODE
  if(gridMode)
  {
    /* if no data then switch to grid mode */
    GRAPHICDEVICE *grdev;

    grdev = (GRAPHICDEVICE *)GRAPE_CALL(GraphicDevice,"get-stddev") ();
    grdev->clear();
    if (grdev && (grdev->grid_patch != G_GRID)) {
      GRAPE_CALL(grdev,"grid-patch") (G_GRID);
    }
  }
#endif

  sc->object = (TREEOBJECT *)mesh;

  if((!maxlevelButton)) setupLeafButton(mgr,sc,0);

  grape_add_remove_methods();

  // if function name exists, this function is selected again
  if( (lastFunctionName != "") && (lastSlotName != "") )
  {
    GRAPE_CALL(mesh,"select-function") (lastSlotName.c_str(),lastFunctionName.c_str());
  }

  GRAPE_CALL(mgr,"handle") (sc);  // grape display call

  removeLeafButton(mgr,sc);

  FctSelector :: DataFunctionName_t fctName =
    FctSelector :: getCurrentFunctionName(mesh);

  // remember last selected function
  lastSlotName     = fctName.first;
  lastFunctionName = fctName.second;

  // remove obj
  sc->object = 0;
  // preserve sc
  sceneStack.push(sc);

  return ;
}

/*
 * setup TimeScene Tree  */
inline void addHmeshToTimeScene(void * timescene, double time, void  *hmesh, int proc)
{
  TIMESCENE *tsc = (TIMESCENE *) timescene;
  GRAPEMESH *mesh = (GRAPEMESH *) hmesh;
  int i=0;
  assert(tsc != NULL); assert ( mesh != NULL);

  for(i=0; i<proc; i++)
  {
    tsc = (TIMESCENE *) tsc->next_scene;
  }
  assert(tsc);

  if(tsc->dynamic)
  {
    tsc->dynamic = (G_SCENE_OBJECT *)GRAPE_CALL(tsc->dynamic,"put") (mesh, mesh, time);
  }
  else
    tsc->dynamic = (G_SCENE_OBJECT *)GRAPE_CALL(TimeStep,"put") (mesh, mesh, time);

  return;
}

/*
 * setup TimeScene Tree  */
inline void addHmeshToGlobalTimeScene(double time, void  *hmesh, int proc)
{
  TIMESCENE *tsc = globalTsc;
  GRAPEMESH *mesh = (GRAPEMESH *) hmesh;
  assert(tsc  != NULL);
  assert(mesh != NULL);

  if(tsc->dynamic)
  {
    tsc->dynamic = (G_SCENE_OBJECT *)GRAPE_CALL(tsc->dynamic,"put") (mesh, mesh, time);
  }
  else
    tsc->dynamic = (G_SCENE_OBJECT *)GRAPE_CALL(TimeStep,"put") (mesh, mesh, time);

  return;
}

inline DUNE_FDATA * extractData ( void * hmesh , int num )
{
  HMESH *mesh = (HMESH *) hmesh;
  int count = 0;
  assert ( mesh != NULL );

  printf("actual mesh %p \n",mesh);

  GENMESH_FDATA *next_data = mesh->f_data;
  while ( (count < num) && (next_data != NULL) )
  {
    next_data = next_data->next;
    count++;
  }

  if( next_data )
  {
    DUNE_FDATA * df = (DUNE_FDATA *) next_data->function_data;
    std::cout << "df->name = " << df->name << std::endl;
    return df;
  }

  return NULL;
}


/* copy function data */
inline void copyFdata(F_DATA *copy, F_DATA *org)
{
  copy->name = org->name;
  copy->last = org->last;
  copy->next = org->next;

  copy->dimension_of_value  = org->dimension_of_value;
  copy->continuous_data = org->continuous_data;
  copy->function_data  = org->function_data;

  copy->f = org->f;
  copy->f_el_info = org->f_el_info;

  copy->get_bounds  = org->get_bounds;
  copy->get_vertex_estimate = org->get_vertex_estimate;
  copy->get_element_estimate  = org->get_element_estimate;

  copy->threshold = org->threshold;
#if GRAPE_DIM == 3
  copy->geometry_threshold = org->geometry_threshold;
#else
  copy->get_element_p_estimates = org->get_element_p_estimates;
  copy->get_edge_p_estimates    = org->get_edge_p_estimates;
#endif
  copy->hp_threshold = org->hp_threshold;
  copy->hp_maxlevel = org->hp_maxlevel;
}

inline static void copyHmeshes(GRAPEMESH *orgMesh, GRAPEMESH * self)
{
  assert( self );
  assert( orgMesh );

  if( (!self->f_data) && (orgMesh->f_data) )
  {
    self->level_of_interest = orgMesh->level_of_interest;
    F_DATA *next_data = (F_DATA *) orgMesh->f_data;
    /* to keep the same order we have to go backward */
    while (next_data)
    {
      if(next_data->next) next_data = (F_DATA *) next_data->next;
      else break;
    }
    while (next_data)
    {
      F_DATA * f_data = (F_DATA *) malloc(sizeof(F_DATA));
      assert(f_data != NULL);
      copyFdata(f_data,next_data);

      self = (GRAPEMESH *) GRAPE_CALL(self,"add-function") (f_data);
      next_data = (F_DATA *) next_data->last;
    }
  }

  self->max_dimension_of_coord = orgMesh->max_dimension_of_coord;
  self->max_eindex = orgMesh->max_eindex;
  self->max_vindex = orgMesh->max_vindex;
  self->max_dindex = orgMesh->max_dindex;
  self->max_number_of_vertices = orgMesh->max_number_of_vertices;

  self->access_mode = orgMesh->access_mode;
  self->access_capability = orgMesh->access_capability;

  /* we have to do that, GRAPE sucks  */
  /* set other function_data pointers */
  if(orgMesh->f_data)
  {
    GENMESH_FDATA * sf = self->f_data;
    while(sf != NULL)
    {
      const char * sfname = sf->name;
      GENMESH_FDATA * nf = orgMesh->f_data;
      int length = strlen(sfname);
      while( (nf != NULL) )
      {
        /* compare the real function name, ha, not with me */
        const char * nfname = nf->name;
        if( strncmp(sfname,nfname,length) == 0 )
        {
          sf->function_data = nf->function_data;
          break;
        }

        /* GRAPE sucks, sucks, sucks
         * Robert after debugin' for this shit more than one day */
        if( nf != nf->last )
          nf = nf->next;
        else
          break;
      }

      /* go next f_data */
      if( sf != sf->last )
        sf = sf->next;
    }
  }

  /* copy current function selections to orgMesh */
  self = (GRAPEMESH *) GRAPE_CALL(self, "copy-function-selector") (orgMesh);

  self->user_data = orgMesh->user_data;

  self->copy_element = orgMesh->copy_element;
  self->free_element = orgMesh->free_element;

  self->complete_element = orgMesh->complete_element;
  self->set_time = orgMesh->set_time;
  self->get_time = orgMesh->get_time;

  self->first_macro = orgMesh->first_macro;
  self->next_macro = orgMesh->next_macro;
  self->first_child = orgMesh->first_child;
  self->next_child = orgMesh->next_child;
  self->select_child = orgMesh->select_child;

  self->max_level = orgMesh->max_level;
  self->level_of_interest = orgMesh->level_of_interest;
  /* do not set level_of_interest, because is set by user during run time */

  self->get_geometry_vertex_estimate  = orgMesh->get_geometry_vertex_estimate;
  self->get_geometry_element_estimate = orgMesh->get_geometry_element_estimate;
  self->get_lens_element_estimate     = orgMesh->get_lens_element_estimate;
  self->threshold                     = orgMesh->threshold;

#if GRAPE_DIM==2
  self->dimension_of_world = orgMesh->dimension_of_world;
#endif
}

/* interpol method for timescence, just constant interpolation */
inline static GRAPEMESH *grape_mesh_interpol(GRAPEMESH *mesh1, GRAPEMESH *mesh2,
                                             double factor)
{
  GRAPEMESH *self=NULL;
  GRAPEMESH *org =NULL;

  self = (GRAPEMESH *)START_METHOD(G_INSTANCE);
  ASSURE (self, "No HMESH in method interpol! \n", END_METHOD (NULL));

  if (factor < 0.5)
    org = mesh1;
  else
    org = mesh2;

  // copy meshes
  copyHmeshes(org,self);

  END_METHOD(self);
}

/****************************************************************************/
/* handling of multiple functions (selection by next/last)                  */
/****************************************************************************/

inline static HMESH *next_f_data_send(void)
{
  HMESH *self;
  printf("next_f_data_send called! \n");

  self = (HMESH *)START_METHOD(G_INSTANCE);
  ASSURE(self, "", END_METHOD(NULL));

  if (self->f_data && self->f_data->next)
  {
    self->f_data->next->last = self->f_data;  /*only to be sure...*/
    self->f_data = self->f_data->next;
  }
  if (self->f_data)
    printf("new f_data is: %s\n", self->f_data->name);

  END_METHOD(self);
}

inline static HMESH *prev_f_data_send(void)
{
  HMESH *self;
  printf("prev_f_data_send called! \n");

  self = (HMESH *)START_METHOD(G_INSTANCE);
  ASSURE(self, "", END_METHOD(NULL));

  if (self->f_data && self->f_data->last)
  {
    self->f_data->last->next = self->f_data;  /*only to be sure...*/
    self->f_data = self->f_data->last;
  }
  if (self->f_data)
    printf("new f_data is: %s\n", self->f_data->name);

  END_METHOD(self);
}

inline SCENE* scene_maxlevel_on_off ()
{
  SCENE* sc = (SCENE*) START_METHOD (G_INSTANCE);
  GRAPE_ALERT(sc,"maxlevel-on-off: No hmesh!",END_METHOD(NULL));

  if( maxlevelButton->on_off == ON )
  {
    GRAPE_CALL(maxlevelButton,"set-state") (UNPRESSED);
    maxlevelButton->on_off = OFF;
  }
  else
  {
    GRAPE_CALL(maxlevelButton,"set-state") (PRESSED);
    maxlevelButton->on_off = ON;
  }
  END_METHOD (sc);
}

inline SCENE* scene_set_min_max_values ()
{
  SCENE* sc = (SCENE*) START_METHOD (G_INSTANCE);
  GRAPE_ALERT( sc, "set-min-max-values: No hmesh!", END_METHOD(NULL));

  // only if method different from display
  if(sc->method_name)
  {
    // remember method name
    grapeMethodName = sc->method_name;
    GRAPE_CALL(sc, "universal") ("value-min-max");
  }
  END_METHOD (sc);
}

GENMESH3D * genmesh3d_switch_iterateLeafs_on_off();

static int ruler_bnd_id = 0;
inline static HELEMENT * bnd_next_macro (HELEMENT * prevEl, MESH_ELEMENT_FLAGS flag)
{
  HELEMENT * el = next_macro(prevEl,flag);
  if(el)
  {
    DUNE_ELEM* elem = (DUNE_ELEM*) el->user_data;
    assert( elem );
    for(int i=0; i<MAX_EL_FACE; ++i)
      if(ruler_bnd_id == elem->bnd[i]) return el;
  }
  else
  {
    return 0;
  }
  return bnd_next_macro(el,flag);
}

inline static HELEMENT * bnd_first_macro (GENMESHnD *mesh, MESH_ELEMENT_FLAGS flag)
{
  HELEMENT * el = first_macro(mesh,flag);

  if(el)
  {
    DUNE_ELEM* elem = (DUNE_ELEM*) el->user_data;
    assert( elem );
    for(int i=0; i<MAX_EL_FACE; ++i)
      if(ruler_bnd_id == elem->bnd[i]) return el;
  }
  else
  {
    return 0;
  }
  return bnd_next_macro(el,flag);
}

inline HMESH* genmesh_boundary_disp ()
{
  HMESH * hmesh = (HMESH*) START_METHOD (G_INSTANCE);
  GRAPE_ALERT(hmesh, "genmesh-boundary-id: No hmesh!", END_METHOD(NULL));

  MANAGER* mgr = (MANAGER *) GRAPE_CALL(Manager,"get-stdmgr") ();
  assert( mgr );

  hmesh->first_macro = bnd_first_macro ;
  hmesh->next_macro  = bnd_next_macro ;

  static RULER * bndIdRuler = 0;

  if(!bndIdRuler)
  {
    bndIdRuler = (RULER *)
                 new_item(Ruler, I_Name, "boundary id",
                          I_Instance, hmesh,
                          I_Var, &ruler_bnd_id, dfINT,
                          I_ColorRGB,0.42,0.42,0.0,
                          I_End);
  }

  // if called first time, add ruler
  if(GRAPE_CALL(mgr,"new-handle") (genmesh_boundary_disp,1))
  {
    GRAPE_CALL(mgr,"add-inter") (bndIdRuler);
  }

  GRAPE_CALL(hmesh,"display") ();

  hmesh->first_macro = first_macro ;
  hmesh->next_macro  = next_macro ;

  END_METHOD (hmesh);
}

// make hard copy of hmesh, i.e. create new object and copy data
inline HMESH* newHmeshHardCopy()
{
  HMESH * hmesh = (HMESH*) START_METHOD (G_INSTANCE);
  GRAPE_ALERT(hmesh, "hmesh-hardcopy: No hmesh!", END_METHOD(NULL));

  std::string newName("H: ");
  newName += hmesh->name;

  // get new hmesh
  HMESH* copy = (HMESH *) GRAPE_CALL(HMesh,"new-instance") (newName.c_str());
  GRAPE_ALERT(copy, "hmesh_hardcopy: No new instance!", END_METHOD(NULL));

  // copy mesh
  copyHmeshes(hmesh,copy);

  END_METHOD (copy);
}

/* add some usefull methods */
inline static void grape_add_remove_methods(void)
{
  // if true returned, then grape was already initialized
  if( GRAPE_CALL(HMesh,"find-method") ("next-f-data-send") ) return ;

  if( (GRAPE_CALL(Scene,"find-method") ("set-min-max-values")) )
  {
    GRAPE_CALL(Scene,"delete-method") ("set-min-max-values");
  }
  GRAPE_CALL(Scene,"add-method") ("set-min-max-values",scene_set_min_max_values);

  printf("add-method 'next-f-data-send' on HMesh%dd!\n",GRAPE_DIM);
  GRAPE_CALL(HMesh,"add-method") ("next-f-data-send",&next_f_data_send);
  printf("add-method 'prev-f-data-send' on HMesh%dd!\n",GRAPE_DIM);
  GRAPE_CALL(HMesh,"add-method") ("prev-f-data-send",&prev_f_data_send);
  GRAPE_CALL(GrapeMesh,"add-method") ("interpol", &grape_mesh_interpol);

  printf("add-method 'value-min-max' on HMesh%dd!\n",GRAPE_DIM);

  GRAPE_CALL(HMesh,"add-method") ("value-min-max", &setMinMaxValue);

  // overload hardcopy
  printf("add-method 'hardcopy' on HMesh%dd!\n",GRAPE_DIM);
  GRAPE_CALL(HMesh,"add-method") ("hardcopy",&newHmeshHardCopy);

#if GRAPE_DIM == 3
  GRAPE_CALL(GenMesh3d,"add-method") ("get-partition-number",get_partition_number);
  GRAPE_CALL(HMesh,"add-method") ("boundary-id-disp",genmesh_boundary_disp);
#endif

  if( ! (GRAPE_CALL(Scene,"find-method") ("maxlevel-on-off")) )
    GRAPE_CALL(Scene,"add-method") ("maxlevel-on-off",scene_maxlevel_on_off);

  {
    char p_name[32];
    sprintf(p_name,"uif-m%d",GRAPE_DIM);
    g_project_add(p_name);
  }
}
#endif
