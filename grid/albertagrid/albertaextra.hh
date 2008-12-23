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

// provides the element number generation and management
#include "agelementindex.cc"

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
    }
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

  // Leaf Data for Albert, only the leaf elements have this data set
  template <int cdim, int vertices>
  struct AlbertLeafData
  {
#ifdef LEAFDATACOORDS
    typedef Dune::FieldMatrix<double,vertices,cdim> CoordinateMatrixType;
    typedef Dune::FieldVector<double,cdim> CoordinateVectorType;
#endif
    // type of stored data
    typedef struct
    {
#ifdef LEAFDATACOORDS
      CoordinateMatrixType coord;
#endif
      double determinant;
    } Data;

    // keep element numbers
    inline static void AlbertLeafRefine( EL *parent, EL *child[2] )
    {
      Data * ldata;
      int i;

      ldata = (Data *) parent->child[1];
      assert(ldata != 0);

      //std::cout << "Leaf refine for el = " << parent << "\n";

      double childDet = 0.5 * ldata->determinant;

      /* bisection ==> 2 children */
      for(i=0; i<2; i++)
      {
        Data *ldataChi = (Data *) child[i]->child[1];
        assert(ldataChi != 0);
        ldataChi->determinant = childDet;

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
      double & det = ldata->determinant;
      det = 0.0;

      //std::cout << "Leaf coarsen for el = " << parent << "\n";

      /* bisection ==> 2 children */
      for(i=0; i<2; i++)
      {
        Data *ldataChi = (Data *) child[i]->child[1];
        assert(ldataChi != 0);
        det += ldataChi->determinant;
      }
    }

#if DUNE_ALBERTA_VERSION < 0x200
    // we dont need Leaf Data
    inline static void initLeafData(LEAF_DATA_INFO * linfo)
    {
      linfo->leaf_data_size = sizeof(Data);
      linfo->refine_leaf_data  = &AlbertLeafRefine;
      linfo->coarsen_leaf_data = &AlbertLeafCoarsen;
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
    }

    // remember on which level an element realy lives
    inline static void initLeafDataValues( MESH * mesh, int proc )
    {
      // see ALBERTA Doc page 72, traverse over all hierarchical elements
      ALBERTA meshTraverse(mesh,-1, CALL_LEAF_EL|FILL_COORDS,setLeafData);
    }

  }; // end of AlbertLeafData


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


#ifndef CALC_COORD
  template< int dimworld >
  inline static void setLocalCoords( const EL_INFO *elInfo )
  {
    const EL *element = elInfo->el;
    assert( element != NULL );

    const DOF_ADMIN *admin = coordVec->fe_space->admin;
    const int nv = admin->n0_dof[ VERTEX ];

    REAL_D *vec = NULL;
    GET_DOF_VEC( vec, coordVec );
    assert( vec != NULL );

    for( int i = 0; i < N_VERTEX(admin->mesh); ++i )
    {
      int dof = element->dof[ i ][ nv ];
      REAL_D &vecCoord = vec[ dof ];
      const REAL_D &coord = elInfo->coord[ i ];
      for( int j = 0; j < dimworld; ++j )
        vecCoord[ j ] = coord[ j ];
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


  // clear Dof Vec
  inline static void clearDofVec ( DOF_INT_VEC * drv )
  {
    int *vec = NULL;
    GET_DOF_VEC( vec, drv );
    FOR_ALL_DOFS( drv->fe_space->admin, vec[ dof ] = 0 );
  }


  inline DOF_INT_VEC * getDofNewCheck(const FE_SPACE * espace,
                                      const char * name)
  {
    DOF_INT_VEC * drv = get_dof_int_vec(name,espace);
    clearDofVec( drv );
    return drv;
  }


#if DUNE_ALBERTA_VERSION >= 0x201
  template< int dim >
  static inline const ALBERTA FE_SPACE *
  getFeSpace( MESH *mesh, const char *name, const int (&ndof)[ 4 ] )
  {
    return get_dof_space( mesh, name, ndof, ADM_PRESERVE_COARSE_DOFS );
  }
#elif DUNE_ALBERTA_VERSION == 0x200
  template< int dim >
  static inline const ALBERTA FE_SPACE *
  getFeSpace( MESH *mesh, const char *name, const int (&ndof)[ 4 ] )
  {
    // the last 1 stands for preserve_coarse_dofs
    return get_fe_space( mesh, name, ndof, 0, 1 );
  }
#else
  template <int dim>
  static inline const ALBERTA FE_SPACE *
  getFeSpace ( MESH *mesh, const char *name, const int (&ndof)[ dim+1 ] )
  {
    // don't delete dofs above the leaf level
    mesh->preserve_coarse_dofs = 1;
    return get_fe_space( mesh, name, ndof, 0 );
  }
#endif


#ifndef CALC_COORD
  // setup of DOF_INT_VEC for reading
  template< int dim, int dimworld >
  inline DOF_REAL_D_VEC *makeTheRest ( MESH *mesh )
  {
    assert( mesh != NULL );

#if DUNE_ALBERTA_VERSION < 0x200
    enum { nNodes = dim+1 };
    enum { vertex = 0, center = dim, face = dim-1, edge = 1 };
#else
    enum { nNodes = N_NODE_TYPES };
    enum { vertex = VERTEX, center = CENTER,
           face = (dim == 3) ? FACE : EDGE, edge = EDGE };
#endif

    int vdof[ nNodes ]; // add at each vertex one dof for vertex numbering
    for( int i = 0; i < nNodes; ++i )
      vdof[ i ] = 0;
    vdof[ vertex ] = 1;

    const FE_SPACE *vSpace = getFeSpace< dim >( mesh, "vertex_dofs", vdof );

    coordVec = get_dof_real_d_vec("coordinates",vSpace);
    coordVec->refine_interpol = &refineCoordsAndRefineCallBack<dimworld>;
    coordVec->coarse_restrict = &coarseCallBack; // coords don't need to be coarsened
    DOF_REAL_D_VEC *coords = getCoordVec< dimworld >();
    coordVec = 0;
    return coords;
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
