// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UG_WRAPPER_3D_HH
#define DUNE_UG_WRAPPER_3D_HH

/** \file
 * \brief Encapsulates some UG macros and functions
 */

/** \todo These macros have been copied from pargm.h to here, because currently there
    are problems with the double inclusion. */
#define PRIO2LISTPART(listtype,prio)                                         \
  ((listtype == ELEMENT_LIST) ? ((prio == PrioHGhost) ? 0 :                \
                                 (prio == PrioVGhost) ? 0 : (prio == PrioVHGhost) ? 0 :               \
                                 (prio == PrioMaster) ? 1 : -1) :                                 \
   ((prio == PrioHGhost) ? 0 : (prio ==PrioVGhost) ? 0 :            \
      (prio == PrioVHGhost) ? 0 :                                  \
      (prio == PrioBorder) ? 2 : (prio == PrioMaster) ? 2 : -1))

#define GRID_ATTR(g) g->level+32

namespace Dune {

  /** \brief Encapsulates a few UG methods and macros
   *
   * This class provides a wrapper to several methods and macros from
   * UG.  There are two reasons for doing this.  First, we don't want
   * to call UG macros directly from DUNE, because they pollute the
   * namespace and therefore we undefine them all.  Secondly,  UG methods
   * appear in the namespaces UG::D2 and UG::D3, but we need the dimension
   * as a template parameter.
   */
  template<>
  class UG_NS<3> {
  public:

    // //////////////////////////////////////////////
    //   Types exported by UG
    // //////////////////////////////////////////////

    typedef UG::D3::RefinementRule RefinementRule;

    typedef UG::D3::CoeffProcPtr CoeffProcPtr;

    typedef UG::D3::UserProcPtr UserProcPtr;

    typedef UG::D3::BndSegFuncPtr BndSegFuncPtr;

    /** \todo This type becomes obsolete as soon as UG implements faces and edges */
    typedef UG::D3::vector Vector;

    typedef UG::D3::multigrid MultiGrid;

    typedef UG::D3::grid Grid;

    typedef UG::D3::edge Edge;

    typedef UG::D3::node Node;

    typedef UG::D3::element Element;

    /** \brief Types of the subentities parametrized by the codimension.  Gets specialized below */
    template <int codim>
    class Entity;

    // //////////////////////////////////////////////
    //   Constants exported by UG
    // //////////////////////////////////////////////

    enum {GM_REFINE_NOT_CLOSED = UG::D3::GM_REFINE_NOT_CLOSED};

    enum {GM_COPY_ALL = UG::D3::GM_COPY_ALL};

    enum {GM_REFINE_TRULY_LOCAL = UG::D3::GM_REFINE_TRULY_LOCAL};

    enum {GM_REFINE_PARALLEL = UG::D3::GM_REFINE_PARALLEL};

    enum {GM_REFINE_NOHEAPTEST = UG::D3::GM_REFINE_NOHEAPTEST};

    enum {NEWEL_CE = UG::D3::NEWEL_CE};

    enum {COARSEN_CE = UG::D3::COARSEN_CE};

    enum {REFINECLASS_CE = UG::D3::REFINECLASS_CE};

    enum {ECLASS_CE = UG::D3::ECLASS_CE};

    enum {RED = UG::D3::RED};

    enum {YELLOW_CLASS = UG::D3::YELLOW_CLASS};

    enum {RED_CLASS = UG::D3::RED_CLASS};

    enum {COARSE = UG::D3::COARSE};

    enum {GM_OK = UG::D3::GM_OK};

    enum {MAX_SONS = UG::D2::MAX_SONS};

    /** \brief The PFIRSTNODE macro which returns the first node in a
     * grid even in a parallel setting.
     */
    static UG_NS<3>::Node* PFirstNode(const UG_NS<3>::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::NODE_LIST;
      using UG::ELEMENT_LIST;
      return PFIRSTNODE(grid);
    }

    /** \brief The FIRSTNODE macro which returns the first node in a
     * grid even in a parallel setting.
     */
    static UG_NS<3>::Node* FirstNode(const UG_NS<3>::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::NODE_LIST;
      using UG::ELEMENT_LIST;
      return FIRSTNODE(grid);
    }

    /** \brief The PFIRSTELEMENT macro which returns the first element in a
     * grid even in a parallel setting.
     */
    static UG_NS<3>::Element* PFirstElement(const UG_NS<3>::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::ELEMENT_LIST;
      return PFIRSTELEMENT(grid);
    }

    /** \brief The FIRSTELEMENT macro which returns the first element in a
     * grid even in a parallel setting.
     */
    static UG_NS<3>::Element* FirstElement(const UG_NS<3>::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::ELEMENT_LIST;
      return FIRSTELEMENT(grid);
    }

    /** \brief Returns pointers to the coordinate arrays of an UG element */
    static void Corner_Coordinates(UG_NS<3>::Element* theElement, double* x[]) {

      using UG::D3::TETRAHEDRON;
      using UG::D3::NODE;
      using UG::D3::PYRAMID;
      using UG::D3::PRISM;
      using UG::D3::HEXAHEDRON;
      using UG::D3::n_offset;
      using UG::UINT;
      int n;    // Dummy variable just to please the macro
      CORNER_COORDINATES(theElement, n, x);
    }

    static int GlobalToLocal(int n, const double** cornerCoords,
                             const double* EvalPoint, double* localCoord) {
      return UG::D3::UG_GlobalToLocal(n, cornerCoords, EvalPoint, localCoord);
    }

    //! return true if element has an exact copy on the next level
    static bool hasCopy (UG_NS<3>::Element* theElement) {
      using UG::D3::ELEMENT;
      using UG::D3::control_entries;
      using UG::UINT;
      return REFINECLASS(theElement) == YELLOW_CLASS;
    }

    //! return true if element has an exact copy on the next level
    static bool isRegular (UG_NS<3>::Element* theElement) {
      using UG::D3::ELEMENT;
      using UG::D3::control_entries;
      using UG::UINT;
      return ECLASS(theElement) == RED_CLASS;
    }

    //! \todo Please doc me!
    static int Sides_Of_Elem(UG_NS<3>::Element* theElement) {
      using UG::D3::element_descriptors;
      using UG::UINT;
      return SIDES_OF_ELEM(theElement);
    }

    //! Encapsulates the NBELEM macro
    static UG_NS<3>::Element* NbElem(UG_NS<3>::Element* theElement, int nb) {
      using UG::D3::ELEMENT;
      using UG::D3::nb_offset;
      using UG::UINT;
      return NBELEM(theElement, nb);
    }

    //! Returns true if the i-th side of the element is on the domain boundary
    static bool Side_On_Bnd(UG_NS<3>::Element* theElement, int i) {
      using UG::D3::BNDS;
      using UG::D3::BEOBJ;
      using UG::D3::side_offset;
      using UG::UINT;
      return OBJT(theElement)==BEOBJ && SIDE_ON_BND(theElement, i);
    }

    //! \todo Please doc me!
    static int Edges_Of_Elem(const UG_NS<3>::Element* theElement) {
      using UG::D3::element_descriptors;
      using UG::UINT;
      return EDGES_OF_ELEM(theElement);
    }

    //! \todo Please doc me!
    static int Corners_Of_Elem(const UG_NS<3>::Element* theElement) {
      using UG::D3::element_descriptors;
      using UG::UINT;
      return CORNERS_OF_ELEM(theElement);
    }

    //! \todo Please doc me!
    // Dummy implementation for vertices
    static int Corners_Of_Elem(const UG_NS<3>::Node* theElement) {
      return 1;
    }

    //! \todo Please doc me!
    static int Corners_Of_Side(const UG_NS<3>::Element* theElement, int side) {
      using UG::D3::element_descriptors;
      using UG::UINT;
      return CORNERS_OF_SIDE(theElement, side);
    }

    //! \todo Please doc me!
    static int Corner_Of_Side(const UG_NS<3>::Element* theElement, int side, int corner) {
      using UG::D3::element_descriptors;
      using UG::UINT;
      return CORNER_OF_SIDE(theElement, side, corner);
    }

    static int nSons(const UG::D3::element* element) {
      return UG::D3::ReadCW(element, UG::D3::NSONS_CE);
    }

    static int myLevel (const UG_NS<3>::Element* theElement) {
      using UG::D3::ELEMENT;
      using UG::UINT;
      return LEVEL(theElement);
    }

    static int myLevel (const UG_NS<3>::Node* theNode) {
      using UG::D3::NODE;
      using UG::UINT;
      return LEVEL(theNode);
    }

    //! get father element of vertex
    static UG_NS<3>::Element* NFather(UG_NS<3>::Node* theNode) {
      return theNode->myvertex->iv.father;
    }

    //! get father node of vertex
    static UG_NS<3>::Node* NodeNodeFather(UG_NS<3>::Node* theNode) {
      using UG::D3::NDOBJ;
      using UG::UINT;
      if (theNode->father==0)
        return 0;         // no father at all
      if (OBJT(theNode->father)==NDOBJ)
        return (UG::D3::node*) theNode->father;
      else
        return 0;         // may be edge or element
    }

    //! get father element of vertex
    static void PositionInFather(UG_NS<3>::Node* theNode, FieldVector<double, 3>& local) {
      local[0] = theNode->myvertex->iv.xi[0];
      local[1] = theNode->myvertex->iv.xi[1];
      local[2] = theNode->myvertex->iv.xi[2];
    }

    //! get father element of vertex
    static void NodePositionGlobal(UG_NS<3>::Node* theNode, FieldVector<double, 3>& global) {
      global[0] = theNode->myvertex->iv.x[0];
      global[1] = theNode->myvertex->iv.x[1];
      global[2] = theNode->myvertex->iv.x[2];
    }

    static int GetSons(const UG::D3::element* element, UG::D3::element* sonList[MAX_SONS]) {
      return UG::D3::GetSons(element, sonList);
    }

    static int GetNodeContext(const UG::D3::element* element, const UG::D3::node** context) {
      return UG::D3::GetNodeContext(element, const_cast<UG::D3::node**>(context));
    }

    //! Encapsulates the GRID_ATTR macro
    static unsigned char Grid_Attr(const UG_NS<3>::Grid* grid) {
      return GRID_ATTR(grid);
    }

    static int MarkForRefinement(UG::D3::element* element, int rule, int data) {
      return UG::D3::MarkForRefinement(element, (UG::D3::RefinementRule)rule, data);
    }

    //! Encapsulates the TAG macro
    static unsigned int Tag(const UG_NS<3>::Element* theElement) {
      using UG::UINT;
      return TAG(theElement);
    }

    //! Doesn't ever get called, but needs to be there to calm the compiler
    static unsigned int Tag(const UG_NS<3>::Node* theNode) {
      DUNE_THROW(GridError, "Called method Tag() for a vertex.  This should never happen!");
      return 0;
    }

    //! get corner in local coordinates, corner number in UG's numbering system
    template<class T>
    static void  getCornerLocal (const UG_NS<3>::Element* theElement, int corner, FieldVector<T, 3>& local)
    {
      using UG::D3::element_descriptors;
      using UG::UINT;
      local[0] = LOCAL_COORD_OF_TAG(TAG(theElement),corner)[0];
      local[1] = LOCAL_COORD_OF_TAG(TAG(theElement),corner)[1];
      local[2] = LOCAL_COORD_OF_TAG(TAG(theElement),corner)[2];
    }

    //! Next element in the UG element lists
    static UG_NS<3>::Element* succ(const UG_NS<3>::Element* theElement) {
      return theElement->ge.succ;
    }

    //! Next element in the UG nodes lists
    static UG_NS<3>::Node* succ(const UG_NS<3>::Node* theNode) {
      return theNode->succ;
    }

    //! Calm the compiler
    static void* succ(const void* theWhatever) {
      DUNE_THROW(NotImplemented, "No successor available for this kind of object");
      return 0;
    }

    //! Return true if the element is a leaf element
    static bool isLeaf(const UG_NS<3>::Element* theElement) {
      return UG::D3::EstimateHere(theElement);
    }

    //! Return true if the node is a leaf node
    static bool isLeaf(const UG_NS<3>::Node* theNode) {
#ifndef ModelP
      return !theNode->son;
#else
      DUNE_THROW(NotImplemented, "isLeaf for nodes in a parallel grid");
#endif
    }

    // /////////////////////////////////////////////
    //   Level indices
    // /////////////////////////////////////////////

    //! Gets the level index of a UG element
    static int& levelIndex(UG_NS<3>::Element* theElement) {
      return theElement->ge.levelIndex;
    }

    //! Gets the level index of a UG element
    static const int& levelIndex(const UG_NS<3>::Element* theElement) {
      return theElement->ge.levelIndex;
    }

    //! Gets the level index of a UG sidevector
    static int& levelIndex(Vector* theVector) {
      return reinterpret_cast<int&>(theVector->index);
    }

    //! Gets the level index of a UG sidevector
    static const int& levelIndex(const Vector* theVector) {
      return reinterpret_cast<const int&>(theVector->index);
    }

    //! Gets the level index of a UG edge
    static int& levelIndex(UG_NS<3>::Edge* theEdge) {
      return theEdge->levelIndex;
    }

    //! Gets the level index of a UG edge
    static const int& levelIndex(const UG_NS<3>::Edge* theEdge) {
      return theEdge->levelIndex;
    }

    //! Gets the level index of a UG node
    static int& levelIndex(UG_NS<3>::Node* theNode) {
      return theNode->levelIndex;
    }

    //! Gets the level index of a UG node
    static const int& levelIndex(const UG_NS<3>::Node* theNode) {
      return theNode->levelIndex;
    }

    // /////////////////////////////////////////////
    //   Leaf indices
    // /////////////////////////////////////////////

    //! Gets the leaf index of a UG element
    static int& leafIndex(UG_NS<3>::Element* theElement) {
      return theElement->ge.leafIndex;
    }

    //! Gets the leaf index of a UG element
    static const int& leafIndex(const UG_NS<3>::Element* theElement) {
      return theElement->ge.leafIndex;
    }

    //! Gets the level index of a UG sidevector
    static int& leafIndex(Vector* theVector) {
      return reinterpret_cast<int &>(theVector->skip);
    }

    //! Gets the level index of a UG sidevector
    static const int& leafIndex(const Vector* theVector) {
      return reinterpret_cast<const int &>(theVector->skip);
    }

    //! Gets the leaf index of a UG edge
    static int& leafIndex(UG_NS<3>::Edge* theEdge) {
      return theEdge->leafIndex;
    }

    //! Gets the leaf index of a UG edge
    static const int& leafIndex(const UG_NS<3>::Edge* theEdge) {
      return theEdge->leafIndex;
    }

    //! Gets the leaf index of a UG node
    static int& leafIndex(UG_NS<3>::Node* theNode) {
      return theNode->myvertex->iv.leafIndex;
    }

    //! Gets the leaf index of a UG node
    static const int& leafIndex(const UG_NS<3>::Node* theNode) {
      return theNode->myvertex->iv.leafIndex;
    }

    // /////////////////////////////////////////////
    //   IDs
    // /////////////////////////////////////////////

    //! Gets the index of a UG element
    static unsigned int id(const UG_NS<3>::Element* theElement) {
      return theElement->ge.id;
    }

    //! Gets the index of a UG node
    static unsigned int id(const UG_NS<3>::Node* theNode) {
      return theNode->myvertex->iv.id | 0xC0000000;
    }

    //! \todo Please doc me!
    static void Local_To_Global(int n, double** y,
                                const FieldVector<double, 3>& local,
                                FieldVector<double, 3>& global) {
      using UG::DOUBLE;
      LOCAL_TO_GLOBAL(n,y,local,global);
    }

    /**
     * \param n Number of corners of the element
     * \param x Coordinates of the corners of the element
     * \param local Local evaluation point
     *
     * \return The return type is int because the macro INVERSE_TRANSFORMATION
     *  returns 1 on failure.
     */
    static int Transformation(int n, double** x,
                              const FieldVector<double, 3>& local, FieldMatrix<double,3,3>& mat) {
      using UG::D3::DOUBLE_VECTOR;
      using UG::DOUBLE;
      double det;
#ifndef SMALL_D
      const double SMALL_D = DBL_EPSILON*10;
#endif
      INVERSE_TRANSFORMATION(n, x, local, mat, det);
      return 0;
    }

    /** \brief Compute the integration element on an element face
        \param nc Number of corners of the element face
        \param co_global Coordinates of the corners of the face
        \param ip_local Local coordinates where integration element is to be evaluated
     */
    static double SurfaceElement(int nc,
                                 const double co_global[MAX_CORNERS_OF_ELEM][3],
                                 const double ip_local[3]) {
      double result;
      if (UG::D3::SurfaceElement(3, nc, co_global, ip_local, &result))
        DUNE_THROW(GridError, "UG::D3::SurfaceElement returned error code!");
      return result;
    }

    //! Returns the i-th corner of a UG element
    static UG_NS<3>::Node* Corner(const UG_NS<3>::Element* theElement, int i) {
      using UG::D3::NODE;
      using UG::D3::n_offset;
      using UG::UINT;
      return CORNER(theElement, i);
    }

    //! get edge from node i to node j (in UG's numbering !
    static UG_NS<3>::Edge* GetEdge (UG_NS<3>::Node* nodei, UG_NS<3>::Node* nodej) {
      return UG::D3::GetEdge(nodei,nodej);
    }

    //! access side vector from element
    static UG_NS<3>::Vector* SideVector (UG_NS<3>::Element* theElement, int i)
    {
      using UG::D3::VECTOR;
      using UG::D3::svector_offset;
      using UG::UINT;
      return SVECTOR(theElement,i);
    }

    //! \todo Please doc me!
    static UG_NS<3>::Element* EFather(UG_NS<3>::Element* theElement) {
      using UG::D3::ELEMENT;
      using UG::D3::father_offset;
      using UG::UINT;
      return EFATHER(theElement);
    }

    static unsigned int ReadCW(void* obj, int ce) {
      return UG::D3::ReadCW(obj, ce);
    }

    static void WriteCW(void* obj, int ce, int n) {
      UG::D3::WriteCW(obj, ce, n);
    }

    //! \todo Please doc me!
    static int InitUg(int* argcp, char*** argvp) {
      return UG::D3::InitUg(argcp, argvp);
    }

    static void ExitUg() {
      UG::D3::ExitUg();
    }

    static void DisposeMultiGrid(UG::D3::multigrid* mg) {
      UG::D3::DisposeMultiGrid(mg);
    }

    //! \todo Please doc me!
    static void* CreateBoundaryValueProblem(const char* BVPname,
                                            int numOfCoeffFunc,
                                            UG::D3::CoeffProcPtr coeffs[],
                                            int numOfUserFct,
                                            UG::D3::UserProcPtr userfct[]) {
      return UG::D3::CreateBoundaryValueProblem(BVPname, 0, numOfCoeffFunc, coeffs,
                                                numOfUserFct, userfct);
    }

    static void* BVP_GetByName(const char* bvpName) {
      return UG::D3::BVP_GetByName(bvpName);
    }

    static void Set_Current_BVP(void** thisBVP) {
      UG::D3::Set_Current_BVP(thisBVP);
    }

    //! \todo Please doc me!
    static UG_NS<3>::MultiGrid* GetMultigrid(const char* name) {
      return UG::D3::GetMultigrid(name);
    }

    //! \todo Please doc me!
    static void SetSubdomain(UG_NS<3>::Element* theElement, int id) {
      using UG::D3::control_entries;
      using UG::D3::SUBDOMAIN_CE;
      using UG::UINT;
      SETSUBDOMAIN(theElement, id);
    }

    static int LBCommand(int argc, const char** argv) {
      /** \todo Can we remove the cast? */
      return UG::D3::LBCommand(argc, (char**)argv);
    }

    static int ConfigureCommand(int argc, const char** argv) {
      /** \todo Kann man ConfigureCommand so ‰ndern daﬂ man auch ohne den const_cast auskommt? */
      return UG::D3::ConfigureCommand(argc, (char**)argv);
    }

    static int NewCommand(int argc, char** argv) {
      return UG::D3::NewCommand(argc, argv);
    }

    static int CreateFormatCmd(int argc, char** argv) {
      return UG::D3::CreateFormatCmd(argc, argv);
    }

    static void* CreateDomain(const char* name, const double* midPoint, double radius,
                              int segments, int corners, int convex) {
      return UG::D3::CreateDomain(name, midPoint, radius, segments, corners, convex);
    }

    static void* InsertInnerNode(UG::D3::grid* grid, const double* pos) {
      return UG::D3::InsertInnerNode(grid, pos);
    }

    static void* CreateBoundarySegment(const char *name, int left, int right,
                                       int index, int res,
                                       int *point,
                                       const double *alpha, const double *beta,
                                       UG::D3::BndSegFuncPtr boundarySegmentFunction,
                                       void *userData) {
      return UG::D3::CreateBoundarySegment(name,            // internal name of the boundary segment
                                           left,                    //  id of left subdomain
                                           right,                    //  id of right subdomain
                                           index,                // Index of the segment
                                           UG::D3::NON_PERIODIC,   // I don't know what this means
                                           res,                    // Resolution, only for the UG graphics
                                           point,
                                           alpha,
                                           beta,
                                           boundarySegmentFunction,
                                           userData);
    }
  };

  template <>
  class UG_NS<3>::Entity<0> {
  public:
    typedef UG::D3::element T;
  };

  template <>
  class UG_NS<3>::Entity<2> {
  public:
    typedef UG::D3::edge T;
  };

  template <>
  class UG_NS<3>::Entity<3> {
  public:
    typedef UG::D3::node T;
  };


} // namespace Dune

#undef PRIO2LISTPART
#undef GRID_ATTR

#endif
