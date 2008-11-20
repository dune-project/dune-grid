// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*  Header--File for extra Albert Functions                                 */
/****************************************************************************/
#ifndef DUNE_ALBERTAEXTRA_HH
#define DUNE_ALBERTAEXTRA_HH

#include <algorithm>
#include <cstring>

#include <dune/grid/albertagrid/albertaheader.hh>

#if HAVE_ALBERTA

#ifdef __ALBERTApp__
namespace Albert {
#endif

#define ALBERTA_ERROR          ALBERTA print_error_funcname(funcName, __FILE__, __LINE__),\
  ALBERTA print_error_msg
#define ALBERTA_ERROR_EXIT     ALBERTA print_error_funcname(funcName, __FILE__, __LINE__),\
  ALBERTA print_error_msg_exit

#define ALBERTA_TEST_EXIT(test) if ((test)) ;else ALBERTA_ERROR_EXIT

#define getDofVec( vec, drv ) \
  (assert(drv != 0); (vec = (drv)->vec); assert(vec != 0));

#if DUNE_ALBERTA_VERSION < 0x200
//! recompute setting of neighbours, because macro_el_info of ALBERTA does
//! that wrong for the Dune context.
inline void computeNeigh(const MACRO_EL *mel, EL_INFO *elinfo, int neigh)
{
  // set right neighbour element
  elinfo->neigh[neigh]      = mel->neigh[neigh]->el;
  // get vertex of opposite coord
  int oppvx = mel->opp_vertex[neigh];
  elinfo->opp_vertex[neigh] = oppvx;

  // copy to opp_coord
  REAL_D *coord  = elinfo->opp_coord;
  const REAL * const * neighcoord  = mel->neigh[neigh]->coord;
  std::memcpy(coord[neigh],neighcoord[oppvx],sizeof(REAL_D));
}
#endif

//! if level iterator is used macro_el_info does not the right thing
inline void fillMacroInfo(TRAVERSE_STACK *stack,
                          const MACRO_EL *mel, EL_INFO *elinfo, int level)
{
  /* Alberta version */
  fill_macro_info(stack->traverse_mesh,mel,elinfo);

#if (DUNE_ALBERTA_VERSION < 0x200) && (DIM == 2)
  // only works for dim 2 at the moment
  // because there we have a different fill_elinfo method
  // quick solution, the method fill_macro_info has to be rewritten
  // not now, dont have the time
  if(level == elinfo->level)
  {
    for(int i=0 ; i<N_VERTEX(stack->traverse_mesh); ++i)
    {
      if(mel->neigh[i])
      {
        computeNeigh(mel,elinfo,i);
      }
      else
      {
        elinfo->neigh[i] = 0;
        elinfo->opp_vertex[i] = 0;
      }
    }
  }
#endif
}


// provides the element number generation and management
#include "agelementindex.cc"

// This three function are used by albertgrid.hh ~.cc
// but not defined in the regular albert.h
//extern void free_leaf_data(void *leaf_data, MESH *mesh);
//extern void free_dof(DOF *dof, MESH *mesh, int position);

inline void enlargeTraverseStack(TRAVERSE_STACK *stack);
inline static TRAVERSE_STACK *getTraverseStack(void);
inline static TRAVERSE_STACK *freeTraverseStack(TRAVERSE_STACK *stack);
inline void printTraverseStack(const TRAVERSE_STACK *stack);

//! organize the TRAVERSE_STACK Management, so we can use the nice Albert
//! functions get_traverse_stack and free_traverse_stack
//! this count the copy made of this class and call free_traverse_stack
//! only if no more copies left
class ManageTravStack
{
  //! traverse stack for mesh traverse, see Albert Docu
  TRAVERSE_STACK * stack_;

  //! number of copies that exist from this stack_
  mutable int *refCount_;

  mutable bool owner_;

public:
  //! initialize the member variables
  ManageTravStack() : stack_ (0) , refCount_ (0) , owner_(false) {}

  //! if a copy is made, the refcout is increased
  ManageTravStack(const ManageTravStack & copy)
  {
    stack_ = 0;
    refCount_ = 0;
    if(copy.stackExists())
    {
      stack_ = copy.stack_;
      refCount_ = copy.refCount_;
      ++(*refCount_);
      copy.owner_ = false;
      owner_ = true;
    }
  }

  //! get new TRAVERSE_STACK using the original Albert Routine
  //! get_traverse_stack, which get an new or free stack
  void create ()
  {
    // remove existing stack, does nothing if no stack exists
    remove();

    assert( stack_ == 0 );
    assert( refCount_ ==  0 );
    stack_ = getTraverseStack();
    refCount_ = new int (1);
    owner_ = true;
  }

  //! set Stack free, if no more refences exist
  ~ManageTravStack()
  {
    remove();
  }

  bool stackExists() const
  {
    return stack_ != 0;
  }

  //! return the TRAVERSE_STACK pointer for use
  TRAVERSE_STACK * getStack() const
  {
    // if this assertion is thrown then either the stack = 0
    // or we want to uese the pointer but are not the owner
    assert( stack_ );
    assert( (!owner_) ? (std::cerr << "\nERROR:The feature of copying iterators is not supported by AlbertaGrid at the moment! \n\n", 0) : 1);
    return stack_;
  }

private:
  //! if copy is made than one more Reference exists
  ManageTravStack & operator = (const ManageTravStack & copy)
  {
    remove();
    // do not use this method
    if(copy.stack_ != 0)
    {
      stack_ = copy.stack_;
      refCount_ = copy.refCount_;
      ++(*refCount_);
      copy.owner_ = false;
      owner_ = true;
    }
    assert(false);
    return (*this);
  }

  void remove()
  {
    if(refCount_ && stack_)
    {
      (*refCount_)--;
      if((*refCount_) <= 0)
      {
        // in free_traverse_stack stack != 0 is checked
        if(stack_)
        {
          stack_ = freeTraverseStack(stack_);
          owner_ = false;
        }
        if(refCount_)
        {
          delete refCount_;
          refCount_ = 0;
        }
      }
    }
    stack_ = 0;
    refCount_ = 0;
  }
};


//***********************************************************
// Traverse Stacks
//***********************************************************
static inline void initTraverseStack(TRAVERSE_STACK *stack);
static inline void resetTraverseStack(TRAVERSE_STACK *stack);

inline static TRAVERSE_STACK *getTraverseStack(void)
{
  TRAVERSE_STACK * stack = get_traverse_stack();

#if DUNE_ALBERTA_VERSION >= 0x200
  initTraverseStack(stack);
#endif

  assert( stack );
  // if we use copyTraverseStack we should only create stacks with
  // stack_size > 0 otherwise we get errors in TreeIterator
  if(stack->stack_size <= 0) enlargeTraverseStack( stack );
  return stack;
}

inline static TRAVERSE_STACK *freeTraverseStack(TRAVERSE_STACK *stack)
{
  // reset stack, i.e set pointer to mesh to 0 ...
  resetTraverseStack(stack);
  free_traverse_stack(stack);
  return 0;
}

inline void copyTraverseStack( TRAVERSE_STACK* stack, const TRAVERSE_STACK* org )
{
  const int & used = stack->stack_size;
  // we have problems to copy stack of length 0
  assert( used > 0 );

  if(stack->elinfo_stack) MEM_FREE(stack->elinfo_stack,used, EL_INFO);
  if(stack->info_stack) MEM_FREE(stack->info_stack,used, U_CHAR );
  if(stack->save_elinfo_stack) MEM_FREE(stack->save_elinfo_stack,used,EL_INFO );
  if(stack->save_info_stack) MEM_FREE(stack->save_info_stack,used,U_CHAR);

  // NOTE: at this point also the used value changes
  // because stack->stack_size is changed
  memcpy( stack, org, sizeof(TRAVERSE_STACK));

  stack->elinfo_stack = 0;
  stack->elinfo_stack = MEM_ALLOC(used, EL_INFO);

  // here we have to copy all EL_INFOs seperately, the normal way does not
  // work, unfortunately
  if (used > 0)
  {
    for (int i=0; i<used; i++)
    {
      memcpy(&(stack->elinfo_stack[i]),&(org->elinfo_stack[i]),sizeof(EL_INFO));
    }
  }

  assert( used == org->stack_size );

  // the pointer have to be created new
  stack->info_stack        = 0;
  stack->info_stack        = MEM_ALLOC(used, U_CHAR);
  stack->save_elinfo_stack = 0;
  stack->save_elinfo_stack = MEM_ALLOC(used, EL_INFO);
  stack->save_info_stack   = 0;
  stack->save_info_stack   = MEM_ALLOC(used, U_CHAR);

  memcpy(stack->elinfo_stack     ,org->elinfo_stack,     used * sizeof(EL_INFO));
  memcpy(stack->info_stack       ,org->info_stack,       used * sizeof(U_CHAR));
  memcpy(stack->save_elinfo_stack,org->save_elinfo_stack,used * sizeof(EL_INFO));
  memcpy(stack->save_info_stack  ,org->save_info_stack,  used * sizeof(U_CHAR));

  /*
     std::cout << "Print original\n";
     printTraverseStack(org);
     std::cout << "Print copy \n";
     printTraverseStack(stack);
   */
  return;
}

static inline void resetTraverseStack(TRAVERSE_STACK *stack)
{
  stack->traverse_mesh = 0;
  stack->traverse_level = 0;
  stack->traverse_mel = 0;
  stack->stack_used = 0;
  stack->save_stack_used = 0;
  stack->el_count = 0;
  return;
}

static inline void initTraverseStack(TRAVERSE_STACK *stack)
{
  // integers and pointers
  stack->traverse_mesh = 0;
  stack->traverse_level = 0;
  stack->traverse_mel = 0;
  stack->el_count = 0;
  stack->stack_used = 0;
  stack->save_stack_used = 0;

  // pointers to stacks and sizes
  stack->elinfo_stack = 0;
  stack->stack_size = 0;
  stack->info_stack = 0;
  stack->save_elinfo_stack = 0;
  stack->save_info_stack = 0;

  // pointer to next stack
  stack->next = 0;
  return;
}

inline void enlargeTraverseStack(TRAVERSE_STACK *stack)
{
  int i;
  int new_stack_size = stack->stack_size + 10;

  stack->elinfo_stack = MEM_REALLOC(stack->elinfo_stack, stack->stack_size,
                                    new_stack_size, EL_INFO);

  if (stack->stack_size > 0)
    for (i=stack->stack_size; i<new_stack_size; i++)
      stack->elinfo_stack[i].fill_flag = stack->elinfo_stack[0].fill_flag;

  stack->info_stack = MEM_REALLOC(stack->info_stack, stack->stack_size,
                                  new_stack_size, U_CHAR);
  stack->save_elinfo_stack = MEM_REALLOC(stack->save_elinfo_stack,
                                         stack->stack_size,
                                         new_stack_size, EL_INFO);
  stack->save_info_stack   = MEM_REALLOC(stack->save_info_stack,
                                         stack->stack_size,
                                         new_stack_size, U_CHAR);

  stack->stack_size = new_stack_size;
}

inline void printTraverseStack(const TRAVERSE_STACK *stack)
{
  FUNCNAME("printTraverseStack");
  MSG("****************************************************\n");
  MSG("current stack %8X | size %d \n", stack,stack->stack_size);
  MSG("traverse_level %d \n",stack->traverse_level);
  MSG("traverse_mesh  %8X \n",stack->traverse_mesh);
  MSG("elinfo_stack      = %8X\n",stack->elinfo_stack);
  MSG("info_stack        = %8X\n",stack->info_stack);
  MSG("save_elinfo_stack = %8X\n",stack->save_elinfo_stack);
  MSG("save_info_stack   = %8X\n\n",stack->save_info_stack);

  MSG("stack_used        = %d\n",stack->stack_used);
  MSG("save_stack_used   = %d\n",stack->save_stack_used);

  MSG("Current elements :\n");
  for(int i=0; i<stack->stack_used+1; ++i)
  {
    //int no = stack->elinfo_stack[i].el->dof[6][0];
    MSG("have element %p \n",stack->elinfo_stack[i].el);
  }

  MSG("****************************************************\n");
}

inline void printElInfo(const EL_INFO *elf)
{
  FUNCNAME("printElInfo");

  MSG("Element %d | level %d  | ",INDEX(elf->el),elf->level);
  printf("Neighs: ");
  for(int i=0; i<N_VERTEX(elf->mesh); i++)
  {
    ALBERTA EL* el = elf->neigh[i];
    printf(" %p |",el);
  }
  printf("\n");


  for(int i=0; i<N_VERTEX(elf->mesh); i++)
    printf("%d %f %f \n",i,elf->coord[i][0],elf->coord[i][1]);


  printf("\n******************************************\n");

}


//****************************************************************
//
//  Wrapper for ALBERTA refine and coarsen routines.
//  Calling direct refine in the grid.refine() method leads to
//  infinite loop. Donno wy?
//  This wrappers solved the problem.
//
//****************************************************************

// wrapper for Albert refine routine
inline static U_CHAR AlbertRefine ( MESH * mesh )
{
  return refine ( mesh );
}

// wrapper for Albert coarsen routine
inline static U_CHAR AlbertCoarsen ( MESH * mesh )
{
  U_CHAR flag = coarsen ( mesh );
  // is mesh was really coarsend, then make dof_compress
  if(flag == MESH_COARSENED) dof_compress ( mesh );
  return flag;
}

//*********************************************************************
//
//  Help Routines for the ALBERTA Mesh
//
//*********************************************************************
namespace AlbertHelp
{

  template <int mydim, int cdim>
  inline void makeEmptyElInfo(EL_INFO * elInfo)
  {
    elInfo->mesh = 0;
    elInfo->el = 0;
    elInfo->parent = 0;
    elInfo->macro_el = 0;
    elInfo->level = 0;
#if DIM > 2
    elInfo->orientation = 0;
    elInfo->el_type = 0;
#endif

    for(int i =0; i<mydim+1; i++)
    {
      for(int j =0; j< cdim; j++)
      {
        elInfo->coord[i][j] = 0.0;
        elInfo->opp_coord[i][j] = 0.0;
      }
      VERTEX_BOUNDARY_ID(elInfo,i) = 0;
    }
  }

  static EL_INFO * getFatherInfo(TRAVERSE_STACK * stack, EL_INFO * elInfo, int level)
  {
    //assert( level == elInfo->level );
    EL_INFO * fatherInfo = 0;

    //std::cout << elInfo->el << " element \n";

    // if this level > 0 return father = elInfoStack -1,
    // else return father = this
    assert(stack != 0);

    //std::cout << "get father of level "<< level << "\n";

    if(level > 0)
    {
      fatherInfo = stack->elinfo_stack + level;
      //std::cout << fatherInfo->el << " father \n";
      //std::cout << fatherInfo->el->child[0] << " ch 0 | ch 1 " << fatherInfo->el->child[1] << "\n";
    }
    else
    {
      assert( (true) ? (printf("No Father for macro element, return macro element\n"),1) : 1);
      return elInfo;
    }
    return fatherInfo;
  }

  //**************************************************************************
  //  calc Maxlevel of AlbertGrid and remember on wich level an element lives
  //**************************************************************************

  static int Albert_MaxLevel_help=-1;

  // function for mesh_traverse, is called on every element
  inline static void calcmxl (const EL_INFO * elf)
  {
    int level = elf->level;
    if(Albert_MaxLevel_help < level) Albert_MaxLevel_help = level;
  }

  // remember on which level an element realy lives
  inline int calcMaxLevel ( MESH * mesh , DOF_INT_VEC * levelVec )
  {
    Albert_MaxLevel_help = -1;

    // see ALBERTA Doc page 72, traverse over all hierarchical elements
    meshTraverse(mesh,-1, CALL_LEAF_EL|FILL_NOTHING, calcmxl);

    // check if ok
    assert(Albert_MaxLevel_help != -1);
    return Albert_MaxLevel_help;
  }



  //**************************************************************************
  inline static void printNeighbour (const EL_INFO * elf)
  {
    int i;
    printf("%d EL \n",INDEX(elf->el));
    for(i=0; i<3; i++)
      if(elf->neigh[i])
        printf("%d Neigh \n",INDEX(elf->neigh[i]));
      else printf("%d Neigh \n",-1);
    printf("----------------------------------\n");
  }

  //*********************************************************************

  static int AlbertaLeafDataHelp_processor = -1;
  // Leaf Data for Albert, only the leaf elements have this data set
  template <int cdim, int vertices>
  struct AlbertLeafData
  {
#ifdef LEAFDATACOORDS
    typedef Dune::FieldMatrix<double,vertices,cdim> CoordinateMatrixType;
    typedef Dune::FieldVector<double,cdim> CoordinateVectorType;
#endif
    // type of stored data
    typedef struct {
#ifdef LEAFDATACOORDS
      CoordinateMatrixType coord;
#endif
      double determinant;
      int processor;
    } Data;

    // keep element numbers
    inline static void AlbertLeafRefine(EL *parent, EL *child[2])
    {
      Data * ldata;
      int i, processor=-1;

      ldata = (Data *) parent->child[1];
      assert(ldata != 0);

      //std::cout << "Leaf refine for el = " << parent << "\n";

      double childDet = 0.5 * ldata->determinant;
      processor = ldata->processor;

      /* bisection ==> 2 children */
      for(i=0; i<2; i++)
      {
        Data *ldataChi = (Data *) child[i]->child[1];
        assert(ldataChi != 0);
        ldataChi->determinant = childDet;
        ldataChi->processor = processor;

#ifdef LEAFDATACOORDS
        // calculate the coordinates
        {
          const CoordinateMatrixType &oldCoord = ldata->coord;
          CoordinateMatrixType &coord = ldataChi->coord;
          for (int j = 0; j < cdim; ++j)
          {
            coord[2][j] = 0.5 * (oldCoord[0][j] + oldCoord[1][j]);
            coord[i  ][j] = oldCoord[2][j];
            coord[1-i][j] = oldCoord[i][j];
          }
          //    std::cout << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
        }
#endif
      }
    }

    inline static void AlbertLeafCoarsen(EL *parent, EL *child[2])
    {
      Data *ldata;
      int i;

      ldata = (Data *) parent->child[1];
      assert(ldata != 0);
      ldata->processor = -1;
      double & det = ldata->determinant;
      det = 0.0;

      //std::cout << "Leaf coarsen for el = " << parent << "\n";

      /* bisection ==> 2 children */
      for(i=0; i<2; i++)
      {
        Data *ldataChi = (Data *) child[i]->child[1];
        assert(ldataChi != 0);
        det += ldataChi->determinant;
        if(ldataChi->processor >= 0)
          ldata->processor = ldataChi->processor;
      }
    }

#if DUNE_ALBERTA_VERSION < 0x200
    // we dont need Leaf Data
    inline static void initLeafData(LEAF_DATA_INFO * linfo)
    {
      linfo->leaf_data_size = sizeof(Data);
      linfo->refine_leaf_data  = &AlbertLeafRefine;
      linfo->coarsen_leaf_data = &AlbertLeafCoarsen;
      return;
    }
#endif

    // function for mesh_traverse, is called on every element
    inline static void setLeafData(const EL_INFO * elf)
    {
      assert( elf->el->child[0] == 0 );
      Data *ldata = (Data *) elf->el->child[1];
      assert(ldata != 0);

#ifdef LEAFDATACOORDS
      for(int i=0; i<vertices; ++i)
      {
        CoordinateVectorType & c = ldata->coord[i];
        const ALBERTA REAL_D & coord = elf->coord[i];
        for(int j=0; j<cdim; ++j) c[j] = coord[j];
        //      std::cout << c << " coord \n";
      }
#endif

      ldata->determinant = ALBERTA el_det(elf);
      ldata->processor   = AlbertaLeafDataHelp_processor;
    }

    // remember on which level an element realy lives
    inline static void initLeafDataValues( MESH * mesh, int proc )
    {
      AlbertaLeafDataHelp_processor = proc;

      // see ALBERTA Doc page 72, traverse over all hierarchical elements
      ALBERTA meshTraverse(mesh,-1, CALL_LEAF_EL|FILL_COORDS,setLeafData);

      AlbertaLeafDataHelp_processor = -1;
    }

  }; // end of AlbertLeafData

  // struct holding the needed DOF_INT_VEC for AlbertGrid
  typedef struct dofvec_stack DOFVEC_STACK;
  struct dofvec_stack
  {
    // storage of unique element numbers
    DOF_INT_VEC * elNumbers[numOfElNumVec];
    // contains information about refine status of element
    DOF_INT_VEC * elNewCheck;
#ifndef CALC_COORD
    // vector that stores the coordinates
    DOF_REAL_D_VEC * coords;
#endif
  };

  static DOF_INT_VEC * elNumbers[numOfElNumVec];
  static DOF_INT_VEC * elNewCheck = 0;
#ifndef CALC_COORD
  static DOF_REAL_D_VEC * coordVec = 0;
#endif

  // return pointer to created elNumbers Vector to mesh
  template <int codim>
  inline DOF_INT_VEC * getElNumbersCodim()
  {
    int * vec = 0;
    GET_DOF_VEC(vec,elNumbers[codim]);
    FOR_ALL_DOFS(elNumbers[codim]->fe_space->admin, vec[dof] = getElementIndexForCodim<codim>() );
    return elNumbers[codim];
  }
  // return pointer to created elNumbers Vector to mesh
  inline DOF_INT_VEC * getElNumbers(const int codim)
  {
    switch (codim)
    {
    case 0 : return getElNumbersCodim<0> ();
    case 1 : return getElNumbersCodim<1> ();
    case 2 : return getElNumbersCodim<2> ();
    case 3 : return getElNumbersCodim<3> ();
    }
    // should not get here
    assert( false );
    abort();
    return 0;
  }

  // return pointer to created elNumbers Vector to mesh
  inline static int calcMaxIndex(DOF_INT_VEC * drv)
  {
    int maxindex = 0;
    int * vec=0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, if(vec[dof] > maxindex) { maxindex = vec[dof] } );
    // we return +1 because this means a size
    return maxindex+1;
  }

  // return pointer to created elNewCheck Vector to mesh
  inline DOF_INT_VEC * getElNewCheck()
  {
    int * vec=0;
    GET_DOF_VEC(vec,elNewCheck);
    FOR_ALL_DOFS(elNewCheck->fe_space->admin, vec[dof] = 0 );
    return elNewCheck;
  }

#ifndef CALC_COORD
  template <int dimworld>
  inline static void setLocalCoords(const EL_INFO * elf)
  {
    const EL * element = elf->el;
    assert(element);
    const DOF_ADMIN  * admin = coordVec->fe_space->admin;
    REAL_D * vec = 0;
    const int nv = admin->n0_dof[VERTEX];

    GET_DOF_VEC(vec,coordVec);
    assert(vec);

    for(int i=0; i<N_VERTEX(admin->mesh); ++i)
    {
      int dof = element->dof[i][nv];
      REAL_D & vecCoord = vec[dof];
      const REAL_D & coord = elf->coord[i];
      for(int j=0; j<dimworld; ++j)
      {
        vecCoord[j] = coord[j];
      }
    }
  }

  // return pointer to created elNewCheck Vector to mesh
  template <int dimworld>
  inline DOF_REAL_D_VEC * getCoordVec()
  {
    MESH * mesh = coordVec->fe_space->admin->mesh;
    assert(mesh);
    meshTraverse(mesh,-1, CALL_EVERY_EL_PREORDER|FILL_COORDS, & setLocalCoords<dimworld> );

    /*
       REAL_D * vec = 0;
       GET_DOF_VEC(vec,coordVec);
       FOR_ALL_DOFS(coordVec->fe_space->admin,
        std::cout << "vec["<<dof<<"] = {";
        for(int i=0; i<dimworld; ++i)
          std::cout << vec[dof][i] << ",";
        std::cout << "}\n";
       );
     */

    return coordVec;
  }
#endif

  template <int dimworld>
  inline void getDofVecs(DOFVEC_STACK * dofvecs)
  {
    for(int i=0; i<numOfElNumVec; i++)
    {
      dofvecs->elNumbers[i]  = getElNumbers(i);
      elNumbers[i]  = 0;
    }

    dofvecs->elNewCheck   = getElNewCheck();  elNewCheck = 0;
#ifndef CALC_COORD
    dofvecs->coords       = getCoordVec<dimworld> ();    coordVec   = 0;
#endif
  }

  struct MeshCallBack
  {
    template <class HandlerImp>
    struct Refinement
    {
      static void apply(void * handler, EL * el)
      {
        assert( handler );
        ((HandlerImp *) handler)->postRefinement(el);
      }
    };

    template <class HandlerImp>
    struct Coarsening
    {
      static void apply(void * handler, EL * el)
      {
        assert( handler );
        ((HandlerImp *) handler)->preCoarsening(el);
      }
    };

    typedef void callBackPointer_t (void * , EL * );

    // pointer to actual mesh, for checking only
    MESH * mesh_;
    // pointer to data handler
    void * dataHandler_;
    // method to cast back and call methods of data handler
    callBackPointer_t * postRefinement_;
    callBackPointer_t * preCoarsening_;

    void reset ()
    {
      mesh_           = 0;
      dataHandler_    = 0;
      postRefinement_ = 0;
      preCoarsening_  = 0;
    }

    const MESH * lockMesh () const { return mesh_; }
    const void * dataHandler () const { return dataHandler_; }

    template <class HandlerImp>
    inline void setPointers(MESH * mesh, HandlerImp & handler)
    {
      mesh_ = mesh; // set mesh pointer for checking
      dataHandler_ = (void *) &handler; // set pointer of data handler

      postRefinement_ = & Refinement<HandlerImp>::apply;
      preCoarsening_  = & Coarsening<HandlerImp>::apply;
    }

    inline void postRefinement( EL * el )
    {
      assert( preCoarsening_ != 0 );
      postRefinement_(dataHandler_,el);
    }
    inline void preCoarsening( EL * el )
    {
      assert( preCoarsening_ != 0 );
      preCoarsening_(dataHandler_,el);
    }

  private:
    MeshCallBack ()
      : mesh_(0), dataHandler_(0), postRefinement_(0), preCoarsening_(0) {}
  public:
    static MeshCallBack & instance()
    {
      static MeshCallBack inst;
      return inst;
    }
  };

#ifndef CALC_COORD
  // set entry for new elements to 1
  template <int dim>
  inline static void refineCoordsAndRefineCallBack ( DOF_REAL_D_VEC * drv , RC_LIST_EL *list, int ref)
  {
    static MeshCallBack & callBack = MeshCallBack::instance();

    const int nv = drv->fe_space->admin->n0_dof[VERTEX];
    REAL_D* vec = 0;
    GET_DOF_VEC(vec,drv);
    assert(ref > 0);

    const EL * el = GET_EL_FROM_LIST(*list);

    // refinement edge is alwyas between vertex 0 and 1
    const int dof0 = el->dof[0][nv];
    const int dof1 = el->dof[1][nv];

    assert( el->child[0] );
    // new dof has local number dim
    const int dofnew = el->child[0]->dof[dim][nv];

    // get coordinates
    const REAL_D & oldCoordZero = vec[dof0];
    const REAL_D & oldCoordOne  = vec[dof1];
    REAL_D & newCoord = vec[dofnew];

    // new coordinate is average between old on same edge
    // see ALBERTA docu page 159, where this method is described
    // as real_refine_inter
    for(int j=0; j<dim; ++j)
      newCoord[j] = 0.5*(oldCoordZero[j] + oldCoordOne[j]);

    if(callBack.dataHandler())
    {
      // make sure that mesh is the same as in MeshCallBack
      assert( drv->fe_space->admin->mesh == callBack.lockMesh() );
      for(int i=0; i<ref; ++i)
      {
        EL * elem = GET_EL_FROM_LIST(list[i]);

        //std::cout << "call refine for element " << elem << "\n";
        callBack.postRefinement(elem);
      }
    }
  }

  inline static void
  coarseCallBack ( DOF_REAL_D_VEC * drv , RC_LIST_EL *list, int ref)
  {
    static MeshCallBack & callBack = MeshCallBack::instance();

    if(callBack.dataHandler())
    {
      assert( drv->fe_space->admin->mesh == callBack.lockMesh() );
      assert(ref > 0);
      for(int i=0; i<ref; ++i)
      {
        EL * el = GET_EL_FROM_LIST(list[i]);
        //std::cout << "call coarse for element " << el << "\n";
        callBack.preCoarsening(el);
      }
    }
  }
#endif

  // set entry for new elements to 1
  inline static void refineElNewCheck ( DOF_INT_VEC * drv , RC_LIST_EL *list, int ref)
  {
    const DOF_ADMIN * admin = drv->fe_space->admin;
    const int nv = admin->n0_dof[CENTER];
    const int k  = admin->mesh->node[CENTER];
    int *vec = 0;

    GET_DOF_VEC(vec,drv);
    assert(ref > 0);

    for(int i=0; i<ref; i++)
    {
      const EL * el = GET_EL_FROM_LIST(list[i]);
      //std::cout << "refine elNewCheck for el = " << el << "\n";

      // get level of father which is the absolute value
      // if value is < 0 then this just means that element
      // was refined
      int level = std::abs( vec[el->dof[k][nv]] ) + 1;
      for(int ch=0; ch<2; ch++)
      {
        // set new to negative level of the element
        // which then can be used to check whether an element is new on not
        // also this vector stores the level of the element which is needed
        // for the face iterator
        vec[el->child[ch]->dof[k][nv]] = -level;
      }
    }
  }

  // clear Dof Vec
  inline static void clearDofVec ( DOF_INT_VEC * drv )
  {
    int * vec=0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, vec[dof] = 0 );
  }

  // calculate max absolute value of given vector
  inline int calcMaxAbsoluteValueOfVector ( const DOF_INT_VEC * drv )
  {
    const int * vec = 0;
    int maxi = 0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, maxi = std::max( maxi , std::abs(vec[dof]) ) );
    return maxi;
  }

  // set all values of vector to its positive value
  inline static void set2positive ( DOF_INT_VEC * drv )
  {
    int * vec=0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, vec[dof] = std::abs( vec[dof] ) );
  }

  // clear Dof Vec
  inline static void setDofVec ( DOF_INT_VEC * drv , int val )
  {
    int * vec=0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, vec[dof] = val );
  }

  // clear Dof Vec
  inline static void copyOwner ( DOF_INT_VEC * drv , int * ownvec )
  {
    int * vec=0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, vec[dof] = ownvec[dof] );
  }


  inline DOF_INT_VEC * getDofNewCheck(const FE_SPACE * espace,
                                      const char * name)
  {
    DOF_INT_VEC * drv = get_dof_int_vec(name,espace);
    int * vec=0;
    drv->refine_interpol = &refineElNewCheck;
    drv->coarse_restrict = 0;
    GET_DOF_VEC(vec,drv);
    FOR_ALL_DOFS(drv->fe_space->admin, vec[dof] = 0 );
    return drv;
  }

  // setup of DOF_INT_VEC for reading
  template <int dimworld>
  inline void makeTheRest(DOFVEC_STACK * dofvecs)
  {
    dofvecs->elNewCheck = getDofNewCheck(dofvecs->elNumbers[0]->fe_space,"el_new_check");
#ifndef CALC_COORD
    {
      MESH * mesh = dofvecs->elNumbers[0]->fe_space->admin->mesh;
      assert( mesh );

      enum { dim = dimworld };
      int vdof[dim+1]; // add at each vertex one dof for vertex numbering

      for(int i=0; i<dim+1; i++)
      {
        vdof[i] = 0;
      }
      vdof[0] = 1;

      const FE_SPACE * vSpace = get_fe_space(mesh, "vertex_dofs", vdof, 0);

      // coords should not exist at this point
      assert( !dofvecs->coords );

      coordVec = get_dof_real_d_vec("coordinates",vSpace);
      coordVec->refine_interpol = &refineCoordsAndRefineCallBack<dimworld>;
      coordVec->coarse_restrict = &coarseCallBack; // coords don't need to be coarsened

      dofvecs->coords = getCoordVec<dimworld> ();
      coordVec = 0;
    }
#endif
  }

#if DUNE_ALBERTA_VERSION < 0x200
  template <int dim>
  static inline const FE_SPACE* getFeSpace(MESH *mesh, const char *name,
                                           const int(&ndof)[dim+1])
  {
    // dont delete dofs on higher levels
    // needed for element numbering
    mesh->preserve_coarse_dofs = 1;
    return get_fe_space( mesh, name, ndof , 0 );
  }
#else
  template <int dim>
  static inline const FE_SPACE* getFeSpace(MESH *mesh, const char *name,
                                           const int(&ndof)[4])
  {
    // the last 1 stands for preserve_coarse_dofs
    return get_fe_space( mesh, name, ndof , 0 , 1);
  }
#endif

  // initialize dofAdmin for vertex and element numbering
  // and a given dim
  // --initDofAdmin
  template <int dim>
  static inline void initDofAdmin(MESH *mesh)
  {
#if DUNE_ALBERTA_VERSION < 0x200
    enum { nNodes = dim+1 };
    enum { vertex = 0, center = dim, face = dim-1, edge = 1 };
#else
    enum { nNodes = N_NODE_TYPES };
    enum { vertex = VERTEX, center = CENTER,
           face = (dim == 3) ? FACE : EDGE, edge = EDGE };
#endif

    int edof[nNodes]; // add one dof at element for element numbering
    int vdof[nNodes]; // add at each vertex one dof for vertex numbering
    int fdof[nNodes]; // add at each face one dof for face numbering
    int edgedof[nNodes]; // add at each edge one dof for edge numbering

    for(int i=0; i<nNodes; ++i)
    {
      vdof[i] = 0; fdof[i]    = 0;
      edof[i] = 0; edgedof[i] = 0;
    }

    vdof[vertex] = 1;

    if(dim == 3)
    {
      edgedof[edge] = 1; // edge dof only in 3d
    }

    fdof[face] = 1; // means edges in 2d and faces in 3d
    edof[center] = 1;

    {
      // !! NOTE:
      // the order in which the refine routines are called is jsut turned
      // arround, which mean that actual the coordinates are refined at last
      // and therefore the element call back is done in this function
#ifndef CALC_COORD
      const FE_SPACE * vSpace =
#endif
      getFeSpace<dim> (mesh, "vertex_dofs", vdof);

#ifndef CALC_COORD
      coordVec = get_dof_real_d_vec("coordinates",vSpace);
      coordVec->refine_interpol = &refineCoordsAndRefineCallBack<dim>;
      coordVec->coarse_restrict = &coarseCallBack; // coords don't need to be coarsened
#endif
    }

    //**********************************************************************
    // all the element vectors
    //**********************************************************************
    // space for center dofs , i.e. element numbers
    const FE_SPACE * elemSpace = getFeSpace<dim> (mesh, "center_dofs", edof);

    // the first entry is called at last for refinement and coarsening
    // the methods for the adaptation call back are put to elNewCheck
    // refine and coarsening procedures
    elNewCheck = get_dof_int_vec("el_new_check",elemSpace);
    elNewCheck->refine_interpol = &refineElNewCheck;
    elNewCheck->coarse_restrict = 0;

    // the element numbers, ie. codim = 0
    elNumbers[0] = get_dof_int_vec("element_numbers",elemSpace);
    elNumbers[0]->refine_interpol = &RefineNumbering<dim,0>::refineNumbers;
    elNumbers[0]->coarse_restrict = &RefineNumbering<dim,0>::coarseNumbers;

    //**********************************************************************

    {
      // the face number space , i.e. codim == 1
      const FE_SPACE * fSpace = getFeSpace<dim> (mesh, "face_dofs", fdof);

      // the face numbers, i.e. codim = 1
      elNumbers[1] = get_dof_int_vec("face_numbers",fSpace);
      elNumbers[1]->refine_interpol = &RefineNumbering<dim,1>::refineNumbers;
      elNumbers[1]->coarse_restrict = &RefineNumbering<dim,1>::coarseNumbers;
    }

    if(dim == 3)
    {
      // the edge number space , i.e. codim == 2
      const FE_SPACE * eSpace = getFeSpace<dim> (mesh, "edge_dofs", edgedof);

      // the edge numbers, i.e. codim = 2
      elNumbers[2] = get_dof_int_vec("edge_numbers",eSpace);
      elNumbers[2]->refine_interpol = &RefineNumbering<dim,2>::refineNumbers;
      elNumbers[2]->coarse_restrict = &RefineNumbering<dim,2>::coarseNumbers;
    }

    return;
  }

  static std::stack < BOUNDARY * > * Alberta_tmpBndStack = 0;

  inline static void initBndStack( std::stack < BOUNDARY * > * bndStack )
  {
    Alberta_tmpBndStack = bndStack;
  }
  inline static void removeBndStack ()
  {
    Alberta_tmpBndStack = 0;
  }

#if DUNE_ALBERTA_VERSION < 0x200
  // initialize boundary for mesh
  inline const BOUNDARY *initBoundary(MESH * Spmesh, int bound)
  {
    BOUNDARY *b = (BOUNDARY *) new BOUNDARY ();
    assert(b != 0);

    assert(Alberta_tmpBndStack);
    Alberta_tmpBndStack->push( b );

    // bound is of type signed char which goes from -127 to 128
    if((bound < -127) && (bound > 128))
    {
      std::cerr << "Got boundary id = " << bound << "\n";
      std::cerr << "Wrong boundary id: range is only from -127 to 128 !\n";
      std::cerr << "Correct your macro grid file!\n";
      abort();
    }

    b->param_bound = 0;
    b->bound = bound;

    return b;
  }
  // create Mesh fof Version 1.2
  template <int dim>
  static MESH* createMesh(const char * name, const char * filename)
  {
    typedef AlbertLeafData<dim,dim+1> LeafDataType;
    MESH * mesh = get_mesh(name,
                           initDofAdmin<dim>,
                           LeafDataType::initLeafData);
    read_macro(mesh, filename, initBoundary);
    return mesh;
  }

#else
  typedef NODE_PROJECTION* initBoundary_t (MESH *,MACRO_EL *,int);
  static initBoundary_t* initBoundary = 0;

  // create Mesh fof Version 2.0
  template <int dim>
  static MESH* createMesh(const char * name, const char * filename)
  {
    // get macro data
    MACRO_DATA* mdata = read_macro(filename);

    // create mesh
#if DUNE_ALBERTA_VERSION >= 0x201
    MESH* mesh = GET_MESH(dim, name, mdata, NULL, NULL );
#else
    MESH* mesh = GET_MESH(dim, name, mdata, NULL );
#endif

    // free macro data
    free_macro_data(mdata);

    // init dof admins
    initDofAdmin<dim> ( mesh );

    //! type of leaf data
    typedef AlbertLeafData<dim,dim+1> LeafDataType;
    init_leaf_data(mesh, sizeof(typename LeafDataType :: Data),
                   LeafDataType :: AlbertLeafRefine,
                   LeafDataType :: AlbertLeafCoarsen);

    return mesh;
  }
#endif

  // mark elements that not belong to my processor
  inline void partitioning ( MACRO_EL * mel, int proc, int mynumber  )
  {
    if(proc == mynumber)
    {
      mel->el->mark = 0;
    }
    else
    {
      mel->el->mark = 1;
    }
  }

  inline void printMacroData(MACRO_DATA * mdata)
  {
    FUNCNAME("printMacroData");
    MSG("noe %d , nvx %d \n",mdata->n_macro_elements,mdata->n_total_vertices);
    for(int i=0; i<mdata->n_total_vertices; i++)
      MSG("coords [%f | %f ]\n",mdata->coords[i][0],mdata->coords[i][1]);

#if DUNE_ALBERTA_VERSION < 0x200
    for(int i=0; i<mdata->n_macro_elements; i++)
      MSG("bound [%d | %d | %d ]\n",mdata->boundary[i][0],mdata->boundary[i][1],mdata->boundary[i][2]);
#endif

  }

  // function for mesh_traverse, is called on every element
  inline static void storeLevelOfElement(const EL_INFO * elf)
  {
    const DOF_ADMIN * admin = elNewCheck->fe_space->admin;
    const int nv = admin->n0_dof[CENTER];
    const int k  = admin->mesh->node[CENTER];
    int *vec = 0;
    const EL * el   = elf->el;

    int level = elf->level;
    if( level <= 0 ) return;

    assert(el);
    GET_DOF_VEC(vec,elNewCheck);

    vec[el->dof[k][nv]] = level;
    return ;
  }

  // remember on which level an element realy lives
  inline void restoreElNewCheck( MESH * mesh, DOF_INT_VEC * elNChk )
  {
    elNewCheck = elNChk;
    assert(elNewCheck != 0);

    // see ALBERTA Doc page 72, traverse over all hierarchical elements
    meshTraverse(mesh,-1,CALL_EVERY_EL_PREORDER|FILL_NEIGH,storeLevelOfElement);
    elNewCheck = 0;
    return ;
  }

} // end namespace AlbertHelp

#ifdef __ALBERTApp__
} // end namespace Albert
#endif

#endif // HAVE_ALBERTA

#endif  /* !_ALBERTAEXTRA_H_ */
