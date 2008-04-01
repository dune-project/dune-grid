// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __ALBERTGRID_ELMEM_CC__
#define __ALBERTGRID_ELMEM_CC__

namespace AlbertHelp {

  // for vertices not needed now
  // dim = 2 means, we have a vector for element and edge numbering
#ifdef ALBERTA_VERSION_12
  enum { numOfElNumVec = DIM };
#else
  enum { numOfElNumVec = DIM_OF_WORLD };
#endif

  // IndexManagerType defined in albertgrid.hh
  // not thread save !!!
  static IndexManagerType * tmpIndexStack[numOfElNumVec];

  inline static void initIndexManager_elmem_cc(IndexManagerType (&newIm)[numOfElNumVec])
  {
    for(int i=0; i<numOfElNumVec; i++)
    {
      tmpIndexStack[i] = &newIm[i];
      assert(tmpIndexStack[i] != 0);
    }
  }

  inline static void removeIndexManager_elmem_cc(int numOfVec)
  {
    for(int i=0; i<numOfVec; i++) tmpIndexStack[i] = 0;
  }

  template <int codim>
  inline static int getElementIndex()
  {
    assert(tmpIndexStack[codim] != 0);
    return (*tmpIndexStack[codim]).getIndex();
  }

  template <int codim>
  inline static int getElementIndexForCodim()
  {
    assert((codim >= 0) && (codim < DIM+1));
    return getElementIndex<codim> ();
  }

  template <int codim>
  inline static void freeElementIndex(int idx)
  {
    assert(tmpIndexStack[codim] != 0);
    (*tmpIndexStack[codim]).freeIndex(idx);
  }

  // codim to ALBERTA Dof Type translator
  template <int dim, int codim> struct AlbertaDofType {
    enum { type  = VERTEX };
  }; // dof located on vertices

  // element dofs
  template <int dim> struct AlbertaDofType<dim,0>
  {
    enum { type  = CENTER };
  }; // dofs located inside an element

  // dim = 2 here, codim = 1
  template <> struct AlbertaDofType<2,1>
  {
    enum { type = EDGE };
  }; // dofs located on edges

  // faces in 3d  (dim = 3 , codim = 1 )
  template <> struct AlbertaDofType<3,1>
  {
    enum { type  = FACE };
  }; // dofs located on faces

  // edges in 3d ( dim = 3 , codim = 2 )
  template <> struct AlbertaDofType<3,2> {
    enum { type  = EDGE };
  }; // dofs located on edges

  //****************************************************************************
  //
  //  if the grid is refined then the two methods on this class
  //  organize new element numbers for all codima and set them free if an
  //  element is removed
  //
  //****************************************************************************

  // default is doing nothing
  template <int codim>
  inline void preserveDofs (int * vec, const int k, const int nv, const EL * father, const int split_face )
  {
    assert(false);
    abort();
  }

  // create new element numbers for children
  template <>
  inline void preserveDofs<0> (int * vec, const int k, const int nv, const EL * el, const int split_face )
  {
    enum { codim = 0 };
    // create two new element numbers
    for(int i=0; i<2; ++i)
    {
      assert( el->child[i] );
      const int dof = el->child[i]->dof[k][nv];
      vec[dof] = getElementIndex<codim>();
    }
  }

  // preserve dofs for faces
  template <>
  inline void preserveDofs<1> (int * vec, const int k, const int nv, const EL * el, const int split_face )
  {
    enum { codim = 1 };

    // face number stays the same
    for(int i=0; i<2; ++i)
    {
      // child 0 gets face number 1 as face 2
      // child 1 gets face number 0 as face 2, see Alberta doc page 4
      assert( el->child[i] );
      const int newdof  = el->child[i]->dof[k + split_face][nv];
      const int olddof  = el->dof[k + (1-i)][nv];
      vec[newdof] = vec[olddof];
    }

    // codim == 1 here
    // create two new face number for face 2 in every element
    for(int i=0; i<2; ++i)
    {
      const EL * ch = el->child[i];
      for(int m=0; m < split_face; ++m)
      {
        const int dof = ch->dof[k + m][nv];
        // the dofs could be already set by neighbours
        // therefore this check is necessary and cannot be optimised
        vec[dof] = getElementIndex<codim>();
      }
    }
  }

  // preserve dofs for edges
  template <>
  inline void preserveDofs<2> (int * vec, const int k, const int nv, const EL * el, const int split_face )
  {
    enum { codim = 2 };

    // to be revised, we need el type here other wise the edges of child 1
    // might be swaped

    // edge numbers that stays the same
    // we have to don nothing about them

    // only adjust the new edges
    // see Alberta doc page 147
    const int newmap[3] = {2 , 4 , 5};
    // codim == 1 here
    // create two new face number for face 2 in every element
    for(int i=0; i<2; ++i)
    {
      const EL * ch = el->child[i];
      for(int m=0; m <3; ++m)
      {
        const int dof = ch->dof[k + newmap[m] ][nv];
        // the dofs could be already set by neighbours
        // therefore this check is necessary and cannot be optimised
        vec[dof] = getElementIndex<codim>();
      }
    }
  }

  template <int dim , int codim>
  struct RefineNumbering
  {
    // get element index form stack or new number
    inline static void refineNumbers ( DOF_INT_VEC * drv , RC_LIST_EL *list, int ref)
    {
      const DOF_ADMIN * admin = drv->fe_space->admin;
      // nv is the number of dofs at one dof locate , i.e. the number of dofs
      // strored at one face , here nv= is always 0

      enum { dtype = AlbertaDofType<dim,codim>::type };
      assert( admin->n0_dof[dtype] == 0 );
      const int nv = 0;

      // k is the off set for the dofs, for example at which point we have
      // face dofs, the offset of vertex dofs is always 0
      const int k = admin->mesh->node[dtype];

      int *vec = 0;
      GET_DOF_VEC(vec,drv);

      assert( ref > 0 );
      assert( codim < 3 );

      for(int i=0; i<ref; ++i)
      {
        const EL * el = GET_EL_FROM_LIST(list[i]);

        // in 3d the old face is face 3 in the new element
        // in 2d its the face 2
        enum { split_face = (dim == 3) ? 3 : 2 };
        preserveDofs<codim> (vec,k,nv,el, split_face );
      }
    }

    // put element index to stack, if element is coarsend
    inline static void coarseNumbers ( DOF_INT_VEC * drv , RC_LIST_EL *list, int ref)
    {
      const DOF_ADMIN * admin = drv->fe_space->admin;
      const int nv = admin->n0_dof    [AlbertaDofType<dim,codim>::type];
      const int k  = admin->mesh->node[AlbertaDofType<dim,codim>::type];
      int *vec = 0;
      GET_DOF_VEC(vec,drv);

      assert(ref > 0);

      for(int i=0; i<ref; i++)
      {
        const EL * el = GET_EL_FROM_LIST(list[i]);
        for(int ch=0; ch<2; ++ch)
        {
          const int dof = el->child[ch]->dof[k][nv];

          // put element index to stack, see elmem.cc
          freeElementIndex<codim>( vec[dof] );
        }
      }
    }
  };

} // end namespace AlbertHelp

#endif
