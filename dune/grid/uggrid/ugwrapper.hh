// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
 * \brief Encapsulates some UG macros and functions
 */

/** \todo Here only to provide the constant DBL_EPSILON.  There's maybe a better way? */
#include "float.h"

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
#if UG_DIM == 2
  template<int dim>
  class UG_NS {};
#endif

  template<>
  class UG_NS< UG_DIM > {
  public:

    // //////////////////////////////////////////////
    //   Types exported by UG
    // //////////////////////////////////////////////

#if UG_DIM == 2
#define UG_NAMESPACE UG::D2
#else
#define UG_NAMESPACE UG::D3
#endif

    enum Priorities {
      PrioNone = UG_NAMESPACE::PrioNone,
      PrioMaster = UG_NAMESPACE::PrioMaster,
      PrioBorder = UG_NAMESPACE::PrioBorder,

      PrioHGhost = UG_NAMESPACE::PrioHGhost,
      PrioVGhost = UG_NAMESPACE::PrioVGhost,
      PrioVHGhost = UG_NAMESPACE::PrioVHGhost
    };


    typedef UG_NAMESPACE ::RefinementRule RefinementRule;

    typedef UG_NAMESPACE ::CoeffProcPtr CoeffProcPtr;

    typedef UG_NAMESPACE ::UserProcPtr UserProcPtr;

    typedef UG_NAMESPACE ::BndSegFuncPtr BndSegFuncPtr;

    /** \brief This is actually a type of the UG algebra, not the grid.
     * We need it to implement face indices and ids in 3d, since UG
     * doesn't actually have objects for faces in 3d grids. */
    typedef UG_NAMESPACE ::vector Vector;

    /** \brief UG type for a hierarchical grid */
    typedef UG_NAMESPACE ::multigrid MultiGrid;

    /** \brief UG type for a level grid */
    typedef UG_NAMESPACE ::grid Grid;

    typedef UG_NAMESPACE ::edge Edge;

    typedef UG_NAMESPACE ::node Node;

    typedef UG_NAMESPACE ::element Element;

    typedef UG_NAMESPACE ::vertex Vertex;

    typedef UG_NAMESPACE ::BVP BVP;

    typedef UG_NAMESPACE ::BVP_DESC BVP_DESC;

    /** \brief Point on a UG boundary patch */
    typedef UG_NAMESPACE ::BNDP BNDP;

    /** \brief Types of the subentities parametrized by the codimension.  Gets specialized below */
    template <int codim>
    class Entity;

    // Type used for local and global ids
#ifdef ModelP
    typedef UG_NAMESPACE::DDD_GID UG_ID_TYPE;
#else
    typedef UG::INT UG_ID_TYPE;
#endif

#ifdef ModelP
    /* DDD Interfaces */
    typedef UG_NAMESPACE::DDD_IF_DIR DDD_IF_DIR;
    typedef UG_NAMESPACE::DDD_IF DDD_IF;
    typedef UG_NAMESPACE::DDD_OBJ DDD_OBJ;
    typedef UG_NAMESPACE::DDD_HEADER DDD_HEADER;

    static void DDD_IFOneway(
#if DUNE_UGGRID_HAVE_DDDCONTEXT
                             DDD::DDDContext& context,
#endif
                             DDD_IF dddIf,
                             DDD_IF_DIR dddIfDir,
                             size_t s,
#if DUNE_UGGRID_HAVE_DDDCONTEXT
                             UG_NAMESPACE::ComProcPtr2 gather,
                             UG_NAMESPACE::ComProcPtr2 scatter
#else
                             UG_NAMESPACE::ComProcPtr gather,
                             UG_NAMESPACE::ComProcPtr scatter
#endif
      )
    {
      UG_NAMESPACE::DDD_IFOneway(
#if DUNE_UGGRID_HAVE_DDDCONTEXT
        context,
#endif
        dddIf, dddIfDir, s, gather, scatter);
    }

#if DUNE_UGGRID_DDD_InfoProcListRange
    static auto DDD_InfoProcListRange(DDD::DDDContext& context, DDD_HEADER* hdr) noexcept
    {
      return UG_NAMESPACE::DDD_InfoProcListRange(context, hdr);
    }
#else
    static int *DDD_InfoProcList(
#if DUNE_UGGRID_HAVE_DDDCONTEXT
                             DDD::DDDContext& context,
#endif
      DDD_HEADER *hdr)
    {
#if DUNE_UGGRID_HAVE_DDDCONTEXT
      //return DDD::DDD_InfoProcList(context, hdr);
      return UG_NAMESPACE::DDD_InfoProcList(context, hdr);
#else
      return UG_NAMESPACE::DDD_InfoProcList(hdr);
#endif
    }
#endif

    static DDD_IF_DIR IF_FORWARD()
    {
      return UG_NAMESPACE::IF_FORWARD;
    }

    static DDD_IF_DIR IF_BACKWARD()
    {
      return UG_NAMESPACE::IF_BACKWARD;
    }

#if DUNE_UGGRID_HAVE_DDDCONTEXT
#  define DDD_CONTEXT_PARAM const DDD::DDDContext& context
#  define DDD_GET_IF(id) (UG_NAMESPACE::ddd_ctrl(context).id)
#else
#  define DDD_CONTEXT_PARAM
#  define DDD_GET_IF(id) (UG_NAMESPACE::id)
#endif

    /*! Master->HGhost/VHGhost */
    static DDD_IF ElementIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(ElementIF);
    }

    /*! ElementSymmIF: Master/HGhost/VHGhost */
    static DDD_IF ElementSymmIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(ElementSymmIF);
    }

    /*! ElementVIF: Master->VGhost/VHGhost */
    static DDD_IF ElementVIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(ElementVIF);
    }

    /*! ElementSymmVIF: Master/VGhost/VHGhost" */
    static DDD_IF ElementSymmVIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(ElementSymmVIF);
    }

    /*! Master->VGhost/HGhost/VHGhost */
    static DDD_IF ElementVHIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(ElementVHIF);
    }

    /*! ElementSymmVHIF: Master/VGhost/HGhost/VHGhost */
    static DDD_IF ElementSymmVHIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(ElementSymmVHIF);
    }

    /*! BorderNodeIF: Border->Master */
    static DDD_IF BorderNodeIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(BorderNodeIF);
    }

    /*! BorderNodeSymmIF: Border/Master */
    static DDD_IF BorderNodeSymmIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(BorderNodeSymmIF);
    }

    /*! OuterNodeIF: Master->HGhost/VGhost */
    static DDD_IF OuterNodeIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(OuterNodeIF);
    }


    /*! NodeVIF: Master->VGhost/VHGhost */
    static DDD_IF NodeVIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(NodeVIF);
    }

    /*! NodeIF: Master->VGhost/HGhost/VHGhost */
    static DDD_IF NodeIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(NodeIF);
    }

    /*! NodeAllIF: All/All */
    static DDD_IF NodeAllIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(NodeAllIF);
    }

    /*! Node_InteriorBorder_All_IF: Master/Border->All */
    static DDD_IF NodeInteriorBorderAllIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(Node_InteriorBorder_All_IF);
    }

    /*! BorderVectorIF: Border->Master */
    static DDD_IF BorderVectorIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(BorderVectorIF);
    }

    /*! BorderVectorSymmIF: Master/Border */
    static DDD_IF BorderVectorSymmIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(BorderVectorSymmIF);
    }

    /*! OuterVectorIF: Master->HGhost/VHGhost */
    static DDD_IF OuterVectorIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(OuterVectorIF);
    }

    /*! OuterVectorSymmIF: Master/Border/HGhost/VHGhost */
    static DDD_IF OuterVectorSymmIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(OuterVectorSymmIF);
    }

    /*! VectorVIF: Master->VGhost/VHGhost */
    static DDD_IF VectorVIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(VectorVIF);
    }

    /*! VectorVAllIF: Master/Border/VGhost/VHGhost->Master/Border */
    static DDD_IF VectorVAllIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(VectorVAllIF);
    }

    /*! VectorIF: Master->VGhost/VHGhost/HGhost */
    static DDD_IF VectorIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(VectorIF);
    }

    static DDD_IF FacetInteriorBorderAllIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(Facet_InteriorBorder_All_IF);
    }

    static DDD_IF FacetAllAllIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(Facet_All_All_IF);
    }

    /*! Master->HGhost/VHGhost */
    static DDD_IF EdgeIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(EdgeIF);
    }

    /*! EdgeSymmIF: Master/HGhost/VHGhost */
    static DDD_IF BorderEdgeSymmIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(BorderEdgeSymmIF);
    }

    /*! EdgeHIF: Master/HGhost/VHGhost */
    static DDD_IF EdgeHIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(EdgeHIF);
    }

    /*! EdgeVHIF: Master->VGhost/HGhost/VHGhost */
    static DDD_IF EdgeVHIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(EdgeVHIF);
    }

    /*! EdgeSymmVHIF: Master/VGhost/HGhost/VHGhost */
    static DDD_IF EdgeSymmVHIF(DDD_CONTEXT_PARAM)
    {
      return DDD_GET_IF(EdgeSymmVHIF);
    }

#undef DDD_CONTEXT_PARAM
#undef DDD_GET_IF

    /** \brief Encapsulates the UG EPRIO macro */
    static int EPriority(const UG_NS< UG_DIM >::Element* element)
    {
      return EPRIO(element);
    }

    /** \brief Returns the priority of the side vector */
    static int Priority(const UG_NS< UG_DIM >::Vector* side)
    {
      return PARHDR(side)->prio;
    }

    /** \brief Returns the priority of the edge (the UG EPRIO macro) */
    static int Priority(const UG_NS< UG_DIM >::Edge* edge)
    {
      return PARHDR(edge)->prio;
    }

    /** \brief Returns the priority of the node (the UG EPRIO macro) */
    static int Priority(const UG_NS< UG_DIM >::Node* node)
    {
      return PARHDR(node)->prio;
    }

    static DDD_HEADER* ParHdr(UG_NS< UG_DIM >::Vector *side)
    {
      return PARHDR(side);
    }

    static DDD_HEADER* ParHdr(UG_NS< UG_DIM >::Edge *edge)
    {
      return PARHDR(edge);
    }

    static DDD_HEADER* ParHdr(UG_NS< UG_DIM >::Node *node)
    {
      return PARHDR(node);
    }

    /** \brief This entry tells the UG load balancer what rank this particular element
     * is supposed to be sent to.
     */
    static UG::INT& Partition(UG_NS< UG_DIM >::Element* element)
    {
      return PARTITION(element);
    }
#endif

    // //////////////////////////////////////////////
    //   Constants exported by UG
    // //////////////////////////////////////////////

    enum {GM_REFINE_NOT_CLOSED = UG_NAMESPACE ::GM_REFINE_NOT_CLOSED};

    enum {GM_COPY_ALL = UG_NAMESPACE ::GM_COPY_ALL};

    enum {GM_REFINE_TRULY_LOCAL = UG_NAMESPACE ::GM_REFINE_TRULY_LOCAL};

    enum {GM_REFINE_PARALLEL = UG_NAMESPACE ::GM_REFINE_PARALLEL};

    enum {GM_REFINE_NOHEAPTEST = UG_NAMESPACE ::GM_REFINE_NOHEAPTEST};

    /** \brief Control word entries */
    enum {NEWEL_CE      = UG_NAMESPACE ::NEWEL_CE,
          COARSEN_CE    = UG_NAMESPACE ::COARSEN_CE,
          ECLASS_CE     = UG_NAMESPACE ::ECLASS_CE,
          MARK_CE       = UG_NAMESPACE ::MARK_CE,
          REFINE_CE     = UG_NAMESPACE ::REFINE_CE};

    /** \brief Refinement rules */
    enum {NO_REFINEMENT = UG_NAMESPACE ::NO_REFINEMENT,
          RED           = UG_NAMESPACE ::RED,
          COARSE        = UG_NAMESPACE ::COARSE};

    enum {RED_CLASS = UG_NAMESPACE ::RED_CLASS};

    enum {GM_OK = UG_NAMESPACE ::GM_OK};

    enum {MAX_SONS = UG_NAMESPACE ::MAX_SONS};

    /** \brief The PFIRSTNODE macro which returns the first node in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Node* PFirstNode(const UG_NS< UG_DIM >::Grid* grid) {
      using UG_NAMESPACE ::PrioHGhost;
      using UG_NAMESPACE ::PrioVGhost;
      using UG_NAMESPACE ::PrioVHGhost;
      using UG_NAMESPACE ::PrioMaster;
      using UG_NAMESPACE ::PrioBorder;
      using UG_NAMESPACE ::ELEMENT_LIST;
      using UG_NAMESPACE ::NODE_LIST;
      return PFIRSTNODE(grid);
    }

    /** \brief The FIRSTNODE macro which returns the first node in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Node* FirstNode(UG_NS< UG_DIM >::Grid* grid) {
      using UG_NAMESPACE ::PrioHGhost;
      using UG_NAMESPACE ::PrioVGhost;
      using UG_NAMESPACE ::PrioVHGhost;
      using UG_NAMESPACE ::PrioMaster;
      using UG_NAMESPACE ::PrioBorder;
      using UG_NAMESPACE ::ELEMENT_LIST;
      using UG_NAMESPACE ::NODE_LIST;
      return FIRSTNODE(grid);
    }

    /** \brief The PFIRSTELEMENT macro which returns the first element in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Element* PFirstElement(const UG_NS< UG_DIM >::Grid* grid) {
      using UG_NAMESPACE ::PrioHGhost;
      using UG_NAMESPACE ::PrioVGhost;
      using UG_NAMESPACE ::PrioVHGhost;
      using UG_NAMESPACE ::PrioMaster;
      using UG_NAMESPACE ::PrioBorder;
      using UG_NAMESPACE ::ELEMENT_LIST;
      using UG_NAMESPACE ::NODE_LIST;
      return PFIRSTELEMENT(grid);
    }

    /** \brief The FIRSTELEMENT macro which returns the first element in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Element* FirstElement(UG_NS< UG_DIM >::Grid* grid) {
      using UG_NAMESPACE ::PrioHGhost;
      using UG_NAMESPACE ::PrioVGhost;
      using UG_NAMESPACE ::PrioVHGhost;
      using UG_NAMESPACE ::PrioMaster;
      using UG_NAMESPACE ::PrioBorder;
      using UG_NAMESPACE ::ELEMENT_LIST;
      return FIRSTELEMENT(grid);
    }

    /** \brief Returns pointers to the coordinate arrays of a UG element */
    static int Corner_Coordinates(const UG_NS< UG_DIM >::Element* theElement, double* x[]) {
      using UG_NAMESPACE ::NODE;
      using UG_NAMESPACE ::TRIANGLE;
      using UG_NAMESPACE ::QUADRILATERAL;
      using UG_NAMESPACE ::TETRAHEDRON;
      using UG_NAMESPACE ::PYRAMID;
      using UG_NAMESPACE ::PRISM;
      using UG_NAMESPACE ::n_offset;
      using UG::UINT;
      int n;
      CORNER_COORDINATES(theElement, n, x);
      return n;
    }

    /** \brief Returns pointers to the coordinate arrays of a UG node */
    static int Corner_Coordinates(const UG_NS< UG_DIM >::Node* theNode, double* x[]) {
      x[0] = theNode->myvertex->iv.x;
      return 1;
    }

    /** \brief Returns pointers to the coordinate arrays of a UG edge */
    static int Corner_Coordinates(const UG_NS< UG_DIM >::Edge* theEdge, double* x[]) {
      x[0] = theEdge->links[0].nbnode->myvertex->iv.x;
      x[1] = theEdge->links[1].nbnode->myvertex->iv.x;
      return 2;
    }

    /** \brief Returns pointers to the coordinate arrays of a UG vector */
    static int Corner_Coordinates(const UG_NS< UG_DIM >::Vector* theVector, double* x[]) {
      UG_NS< UG_DIM >::Element* center;
      unsigned int side;
      UG_NS< UG_DIM >::GetElementAndSideFromSideVector(theVector, center, side);
      int n = Corners_Of_Side(center, side);
      for (int i = 0; i < n; i++)
      {
        unsigned idxInElem = Corner_Of_Side(center, side, i);
        x[i] = Corner(center, idxInElem)->myvertex->iv.x;
      }
      return n;
    }

    static int GlobalToLocal(int n, const double** cornerCoords,
                             const double* EvalPoint, double* localCoord) {
      if (UG_DIM==2)
        // in 2d we can call this only for triangles and quadrilaterals
        assert(n==3 or n==4);
      else
        // in 3d: tetrahedra, pyramids, prisms, hexahedra
        assert(n==4 or n==5 or n==6 or n==8);
      return UG_NAMESPACE ::UG_GlobalToLocal(n, cornerCoords, EvalPoint, localCoord);
    }

    /** \brief Computes the element volume */
    static double Area_Of_Element(int n, const double** cornerCoords) {
      double area = 0.0;
      using UG::DOUBLE;
      using UG_NAMESPACE ::DOUBLE_VECTOR;
#if UG_DIM == 2
      AREA_OF_ELEMENT_2D(n,cornerCoords,area);
#else
      AREA_OF_ELEMENT_3D(n,cornerCoords,area);
#endif
      return area;
    }

    static int myLevel (const UG_NS< UG_DIM >::Element* theElement) {
      using UG::UINT;
      return LEVEL(theElement);
    }

    static int myLevel (const UG_NS< UG_DIM >::Node* theNode) {
      using UG::UINT;
      return LEVEL(theNode);
    }

    static int myLevel (const UG_NS< UG_DIM >::Edge* theEdge) {
      using UG::UINT;
      return LEVEL(theEdge);
    }

    static int myLevel (const UG_NS< UG_DIM >::Vector* theVector) {
      return myLevel((UG_NS< UG_DIM >::Element*)VOBJECT(theVector));
    }

    //! return true if element has an exact copy on the next level
    static bool hasCopy (const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::control_entries;
      using UG::UINT;
      using UG_NAMESPACE ::REFINECLASS_CE;
      using UG_NAMESPACE ::YELLOW_CLASS;
      return REFINECLASS(theElement) == YELLOW_CLASS;
    }

    //! Returns true if element is on level 0 or has been created by red refinement
    static bool isRegular (const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::control_entries;
      using UG::UINT;
      return ECLASS(theElement) == RED_CLASS;
    }

    //! \todo Please doc me!
    static int Sides_Of_Elem(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return SIDES_OF_ELEM(theElement);
    }

    //! Encapsulates the NBELEM macro
    static UG_NS<UG_DIM>::Element* NbElem(const UG_NS< UG_DIM >::Element* theElement, int nb) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::nb_offset;
      using UG::UINT;
      return NBELEM(theElement, nb);
    }

    static size_t boundarySegmentIndex(const UG_NS< UG_DIM >::Element* theElement, int nb) {
      using UG_NAMESPACE ::BNDS;
      using UG::UINT;
      using UG_NAMESPACE ::side_offset;

      BNDS* bnds = ELEM_BNDS(theElement,nb);
      size_t id = UG_NAMESPACE ::GetBoundarySegmentId(bnds);
      return id;
    }

    //! Returns true if the i-th side of the element is on the domain boundary
    static bool Side_On_Bnd(const UG_NS< UG_DIM >::Element* theElement, int i) {
      using UG_NAMESPACE ::BNDS;
      using UG_NAMESPACE ::BEOBJ;
      using UG_NAMESPACE ::side_offset;
      using UG::UINT;
      using UG_NAMESPACE ::GM_OBJECTS;
      return OBJT(theElement)==BEOBJ && SIDE_ON_BND(theElement, i);
    }

    //! Returns true if at least one face of the element is a boundary face
    static bool isBoundaryElement(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::BEOBJ;
      using UG_NAMESPACE ::GM_OBJECTS;
      using UG::UINT;
      return OBJT(theElement)==BEOBJ;
    }

    //! \todo Please doc me!
    static int Edges_Of_Elem(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return EDGES_OF_ELEM(theElement);
    }

    //! \todo Please doc me!
    static int Corners_Of_Elem(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNERS_OF_ELEM(theElement);
    }

    /** \Brief the 'number of corners' of a vertex, i.e., 1.  Here for consistency
        \return 1
     */
    static int Corners_Of_Elem(const UG_NS< UG_DIM >::Node* theNode) {
      return 1;
    }

    //! Return number of corners of a given side
    static int Corners_Of_Side(const UG_NS< UG_DIM >::Element* theElement, int side) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNERS_OF_SIDE(theElement, side);
    }

    //! Return local number of a given corner of a given element side
    static int Corner_Of_Side(const UG_NS< UG_DIM >::Element* theElement, int side, int corner) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNER_OF_SIDE(theElement, side, corner);
    }

    //! Return local number of a given corner of a given element edge
    static int Corner_Of_Edge(const UG_NS< UG_DIM >::Element* theElement, int edge, int corner) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNER_OF_EDGE(theElement, edge, corner);
    }

    //! Return number of sons of an element
    static int nSons(const UG_NAMESPACE ::element* element) {
      return UG_NAMESPACE ::ReadCW(element, UG_NAMESPACE ::NSONS_CE);
    }

    static int GetSons(const UG_NAMESPACE ::element* element, UG_NAMESPACE ::element* sonList[MAX_SONS]) {
      return UG_NAMESPACE ::GetSons(element, sonList);
    }

    /** \todo Remove the const casts */
    static int GetNodeContext(const UG_NAMESPACE ::element* element, const UG_NAMESPACE ::node** context) {
      return UG_NAMESPACE ::GetNodeContext(element, const_cast<UG_NAMESPACE ::node**>(context));
    }

    //! Encapsulates the GRID_ATTR macro
    static int Grid_Attr(const UG_NS< UG_DIM >::Grid* grid) {
      return GRID_ATTR(grid);
    }

    static int MarkForRefinement(UG_NAMESPACE ::element* element, int rule, int data) {
      return UG_NAMESPACE ::MarkForRefinement(element, (UG_NAMESPACE ::RefinementRule)rule, data);
    }

    //! Encapsulates the TAG macro
    static unsigned int Tag(const UG_NS< UG_DIM >::Element* theElement) {
      using UG::UINT;
      return TAG(theElement);
    }

    //! Doesn't ever get called, but needs to be there to calm the compiler
    static unsigned int Tag(const UG_NS< UG_DIM >::Node* theNode) {
      DUNE_THROW(GridError, "Called method Tag() for a vertex.  This should never happen!");
      return 0;
    }

    //! get corner in local coordinates, corner number in UG's numbering system
    template<class T>
    static void  getCornerLocal (const UG_NS< UG_DIM >::Element* theElement, int corner, FieldVector<T, UG_DIM>& local)
    {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      for (int i=0; i<UG_DIM; i++)
        local[i] = LOCAL_COORD_OF_TAG(TAG(theElement),corner)[i];
    }

    //! Next element in the UG element lists
    static UG_NS< UG_DIM >::Element* succ(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.succ;
    }

    //! Next element in the UG nodes lists
    static UG_NS< UG_DIM >::Node* succ(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->succ;
    }

    //! Calm the compiler
    static void* succ(const void* theWhatever) {
      DUNE_THROW(NotImplemented, "No successor available for this kind of object");
      return 0;
    }

    //! Return true if the element is a ghost element
#ifdef ModelP
    static bool isGhost(const UG_NS< UG_DIM >::Element* theElement) {
      if (EPRIO(theElement) == PrioHGhost
          || EPRIO(theElement) == PrioVGhost
          || EPRIO(theElement) == PrioVHGhost)
        return true;
      else
        return false;
    }
#endif

    //! Return true if the element is a leaf element
    static bool isLeaf(const UG_NS< UG_DIM >::Element* theElement) {
      using UG ::UINT;
      using UG_NAMESPACE ::CONTROL_ENTRY;
      using UG_NAMESPACE ::control_entries;

      return LEAFELEM(theElement);
    }

    //! Return true if the node is a leaf node
    static bool isLeaf(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->isLeaf;
    }

    //! Return true if the edge is a leaf edge
    static bool isLeaf(const UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->leafIndex > -1;
    }

    //! Return true if the side vector is a leaf side vector
    static bool isLeaf(const UG_NS< UG_DIM >::Vector* theVector) {
      using UG_NAMESPACE ::VECTOR;
      using UG::UINT;
      // Since the vector cannot be asked directly,
      // the corresponding element is asked.
      return isLeaf((UG_NS< UG_DIM >::Element*)VOBJECT(theVector));
    }
    // /////////////////////////////////////////////
    //   Level indices
    // /////////////////////////////////////////////

    //! Gets the level index of a UG element
    static int& levelIndex(UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.levelIndex;
    }

    //! Gets the level index of a UG element
    static const int& levelIndex(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.levelIndex;
    }

    //! Gets the level index of a UG sidevector
    static UG::UINT& levelIndex(Vector* theVector) {
#if UG_DIM == 2
      DUNE_THROW(GridError, "levelIndex in side vector only in 3D!");
#endif
      return theVector->index;
    }

    //! Gets the level index of a UG sidevector
    static const UG::UINT& levelIndex(const Vector* theVector) {
#if UG_DIM == 2
      DUNE_THROW(GridError, "levelIndex in side vector only in 3D!");
#endif
      return theVector->index;
    }

    //! Gets the level index of a UG edge
    static int& levelIndex(UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->levelIndex;
    }

    //! Gets the level index of a UG edge
    static const int& levelIndex(const UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->levelIndex;
    }

    //! Gets the level index of a UG node
    static int& levelIndex(UG_NS< UG_DIM >::Node* theNode) {
      return theNode->levelIndex;
    }

    //! Gets the level index of a UG node
    static const int& levelIndex(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->levelIndex;
    }

    // /////////////////////////////////////////////
    //   Leaf indices
    // /////////////////////////////////////////////

    //! Gets the leaf index of a UG element
    static int& leafIndex(UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.leafIndex;
    }

    //! Gets the leaf index of a UG element
    static const int& leafIndex(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.leafIndex;
    }

    //! Gets the leaf index of a UG sidevector
    static UG::UINT& leafIndex(Vector* theVector) {
      return theVector->leafIndex;
    }

    //! Gets the leaf index of a UG sidevector
    static const UG::UINT& leafIndex(const Vector* theVector) {
      return theVector->leafIndex;
    }

    //! Gets the leaf index of a UG edge
    static int& leafIndex(UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->leafIndex;
    }

    //! Gets the leaf index of a UG edge
    static const int& leafIndex(const UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->leafIndex;
    }

    //! Gets the leaf index of a UG node
    static int& leafIndex(UG_NS< UG_DIM >::Node* theNode) {
      return theNode->myvertex->iv.leafIndex;
    }

    //! Gets the leaf index of a UG node
    static const int& leafIndex(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->myvertex->iv.leafIndex;
    }

    // /////////////////////////////////////////////
    //   IDs
    // /////////////////////////////////////////////

    //! Gets the index of a UG element
    static auto id(const UG_NS< UG_DIM >::Element* theElement) {
#if defined ModelP
      return theElement->ge.ddd.gid;
#else
      return theElement->ge.id;
#endif
    }

    //! Gets the id of a UG facet
    static auto id(const UG_NS< UG_DIM >::Vector* theVector) {
#ifdef ModelP
      return theVector->ddd.gid;
#else
      auto id = theVector->id;

      // In sequential UG, two entities of different codimension can have the same id,
      // so let's encode the codimension in the id to make them differ.
      constexpr unsigned int codim = 1;  //
      return id | (codim << 30);
#endif
    }

    //! Gets the id of a UG edge
    static auto id(const UG_NS< UG_DIM >::Edge* theEdge) {
#ifdef ModelP
      return theEdge->ddd.gid;
#else
      auto id = theEdge->id;

      // In sequential UG, two entities of different codimension can have the same id,
      // so let's encode the codimension in the id to make them differ.
      constexpr unsigned int codim = UG_DIM-1;  //
      return id | (codim << 30);
#endif
    }

    //! Gets the index of a UG node
    static auto id(const UG_NS< UG_DIM >::Node* theNode) {
#ifdef ModelP
      return theNode->myvertex->iv.ddd.gid;
#else
      auto id = theNode->myvertex->iv.id;

      // In sequential UG, two entities of different codimension can have the same id,
      // so let's encode the codimension in the id to make them differ.
      constexpr unsigned int codim = UG_DIM;
      return id | (codim << 30);
#endif
    }

    /** \brief Compute global coordinates of a point in local coordinates for an element */
    static void Local_To_Global(int n, double** y,
                                const FieldVector<double, UG_DIM>& local,
                                FieldVector<double, UG_DIM>& global) {
      using UG::DOUBLE;
      LOCAL_TO_GLOBAL(n,y,local,global);
    }

    /** \brief Compute global coordinates of a point in local coordinates for a vertex */
    static void Local_To_Global(int n, double** y,
                                const FieldVector<double,0>& local,
                                FieldVector<double, UG_DIM>& global) {
      for (int i=0; i<UG_DIM; i++)
        global[i] = y[0][i];
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
                              const FieldVector<double, UG_DIM>& local, FieldMatrix<double,UG_DIM,UG_DIM>& mat) {
      using UG_NAMESPACE ::DOUBLE_VECTOR;
      using UG::DOUBLE;
      double det;
#ifndef SMALL_D
      const double SMALL_D = DBL_EPSILON*10;
#endif
      INVERSE_TRANSFORMATION(n, x, local, mat, det);
      return 0;
    }

    /** \brief Dummy method for vertices */
    static int Transformation(int n, double** x,
                              const FieldVector<double, 0>& local, FieldMatrix<double,UG_DIM,0>& mat) {
      return 0;
    }

    /**
     * \param n Number of corners of the element
     * \param x Coordinates of the corners of the element
     * \param local Local evaluation point
     *
     * \return The return type is int because the macro TRANSFORMATION
     *  returns 1 on failure.
     */
    static int JacobianTransformation(int n, double** x,
                                      const FieldVector<double, UG_DIM>& local, FieldMatrix<double,UG_DIM,UG_DIM>& mat) {
      using UG_NAMESPACE ::DOUBLE_VECTOR;
      using UG::DOUBLE;
      TRANSFORMATION(n, x, local, mat);
      return 0;
    }

    /** \brief Dummy method for vertices */
    static int JacobianTransformation(int n, double** x,
                                      const FieldVector<double, 0>& local, FieldMatrix<double,0,UG_DIM>& mat) {
      return 0;
    }

    //! Returns the i-th corner of a UG element
    static UG_NS< UG_DIM >::Node* Corner(const UG_NS< UG_DIM >::Element* theElement, int i) {
      using UG_NAMESPACE ::NODE;
      using UG_NAMESPACE ::n_offset;
      using UG::UINT;
      return CORNER(theElement, i);
    }

    //! Returns the i-th edge of a UG element
    static UG_NS< UG_DIM >::Edge* ElementEdge(const UG_NS< UG_DIM >::Element* theElement, int i) {
      using UG_NAMESPACE ::NODE;
      using UG_NAMESPACE ::n_offset;
      using UG::UINT;
      using UG_NAMESPACE ::element_descriptors;
      return GetEdge(CORNER(theElement, CORNER_OF_EDGE(theElement, i, 0)),
                     CORNER(theElement, CORNER_OF_EDGE(theElement, i, 1)));
    }

    //! get edge from node i to node j (in UG's numbering !
    static UG_NS< UG_DIM >::Edge* GetEdge (UG_NS< UG_DIM >::Node* nodei, UG_NS< UG_DIM >::Node* nodej) {
      return UG_NAMESPACE ::GetEdge(nodei,nodej);
    }

    //! access side vector from element (this is just a dummy to compile code also in 2d)
    static Vector* SideVector (const UG_NS< UG_DIM >::Element* theElement, int i)
    {
#if UG_DIM == 2
      DUNE_THROW(GridError, "side vector only in 3D!");
#else
      using UG::D3::VECTOR;
      using UG::D3::svector_offset;
      using UG::UINT;
      return SVECTOR(theElement,i);
#endif
    }

    //! Access element from side vector
    static void GetElementAndSideFromSideVector(const UG_NS< UG_DIM >::Vector* theVector,
                                                UG_NS< UG_DIM >::Element*& theElement,
                                                unsigned int& side)
    {
      using UG_NAMESPACE ::VECTOR;
      using UG::UINT;

      theElement = (UG_NS< UG_DIM >::Element*)VOBJECT(theVector);
      side = VECTORSIDE(theVector);
    }

    /** \brief Return a pointer to the father of the given element */
    static UG_NS< UG_DIM >::Element* EFather(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::father_offset;
      using UG::UINT;
      return EFATHER(theElement);
    }

    //! get father element of vertex
    static UG_NS< UG_DIM >::Element* NFather(UG_NS< UG_DIM >::Node* theNode) {
      return theNode->myvertex->iv.father;
    }

    //! get father node of vertex
    static UG_NS< UG_DIM >::Node* NodeNodeFather(UG_NS< UG_DIM >::Node* theNode) {
      using UG_NAMESPACE ::NDOBJ;
      using UG_NAMESPACE ::GM_OBJECTS;
      using UG::UINT;
      if (theNode->father==0)
        return 0;         // no father at all
      if (OBJT(theNode->father)==NDOBJ)
        return (UG_NAMESPACE ::node*) theNode->father;
      else
        return 0;         // may be edge or element
    }

    static unsigned int ReadCW(void* obj, int ce) {
      return UG_NAMESPACE ::ReadCW(obj, ce);
    }

    static void WriteCW(void* obj, int ce, int n) {
      UG_NAMESPACE ::WriteCW(obj, ce, n);
    }

    //! \todo Please doc me!
    static int InitUg(int* argcp, char*** argvp) {
      return UG_NAMESPACE ::InitUg(argcp, argvp);
    }

    static void ExitUg() {
      UG_NAMESPACE ::ExitUg();
    }

    static int DisposeMultiGrid(UG_NAMESPACE ::multigrid* mg) {
      return UG_NAMESPACE ::DisposeMultiGrid(mg);
    }

    //! \todo Please doc me!
    static void* CreateBoundaryValueProblem(const char* BVPname,
                                            int numOfCoeffFunc,
                                            UG_NAMESPACE ::CoeffProcPtr coeffs[],
                                            int numOfUserFct,
                                            UG_NAMESPACE ::UserProcPtr userfct[]) {
      return UG_NAMESPACE ::CreateBoundaryValueProblem(BVPname, 0, numOfCoeffFunc, coeffs,
                                                       numOfUserFct, userfct);
    }

    //! Set the current boundary value problem
    static void Set_Current_BVP(void** thisBVP) {
      UG_NAMESPACE ::Set_Current_BVP(thisBVP);
    }

    //! Get UG boundary value problem from its name
    static void** BVP_GetByName(const char* name) {
      return UG_NAMESPACE ::BVP_GetByName(name);
    }

    //! Dispose of a boundary value problem
    static int BVP_Dispose(void** BVP) {
      return UG_NAMESPACE ::BVP_Dispose(BVP);
    }

    /** \brief Create new point on the grid boundary by giving local coordinates.

       Global coordinates are computed automatically using the domain boundary information.
       The point is not inserted as a new vertex in the grid!
     */
    static BNDP *BNDP_CreateBndP(UG::HEAP *Heap,
                                 BNDP *theBndP0,
                                 BNDP *theBndP1,
                                 UG::DOUBLE lcoord) {
      return UG_NAMESPACE::BNDP_CreateBndP(Heap, theBndP0, theBndP1, lcoord);
    }

    /** \brief Get global position of a point on the grid boundary */
    static UG::INT BNDP_Global(BNDP *theBndP,
                               UG::DOUBLE *global) {
      return UG_NAMESPACE::BNDP_Global(theBndP, global);
    }

    /** \brief Delete a grid boundary point */
    static UG::INT BNDP_Dispose(UG::HEAP *Heap,
                                BNDP *theBndP) {
      return UG_NAMESPACE::BNDP_Dispose(Heap, theBndP);
    }

    //! Get UG multigrid object from its name
    static UG_NS< UG_DIM >::MultiGrid* GetMultigrid(const char* name) {
      return UG_NAMESPACE ::GetMultigrid(name);
    }

    /** \brief Load-balance the grid by recursive coordinate bisection */
    static void lbs(const char *argv, UG_NAMESPACE ::multigrid *theMG) {
#ifdef ModelP
      return UG_NAMESPACE ::lbs(argv, theMG);
#endif
    }

    //! An UG-internal load balancing method
    static int TransferGridFromLevel(UG_NAMESPACE ::multigrid *theMG, int level) {
#ifdef ModelP
      return UG_NAMESPACE ::TransferGridFromLevel(theMG,level);
#else
      return 0;
#endif
    }

    static MultiGrid *CreateMultiGrid(const char *MultigridName, const char *BndValProblem,
                                      const char *format,
                                      int optimizedIE, int insertMesh,
                                      std::shared_ptr<PPIF::PPIFContext> ppifContext = nullptr) {
      return UG_NAMESPACE ::CreateMultiGrid(const_cast<char*>(MultigridName),
                                            const_cast<char*>(BndValProblem), format,
                                            optimizedIE, insertMesh, ppifContext);
    }

    static void* CreateDomain(const char* name, int segments, int corners) {
      return UG_NAMESPACE ::CreateDomain(name, segments, corners);
    }

    static void RemoveDomain(const char* name) {
      UG_NAMESPACE ::RemoveDomain(name);
    }

    static void* InsertInnerNode(UG_NAMESPACE ::grid* grid, const double* pos) {
      return UG_NAMESPACE ::InsertInnerNode(grid, pos);
    }

    static void* CreateBoundarySegment(const char *name, int left, int right,
                                       int index,
                                       UG::INT *point,
                                       const double *alpha, const double *beta,
                                       UG_NAMESPACE ::BndSegFuncPtr boundarySegmentFunction,
                                       void *userData) {
      return UG_NAMESPACE ::CreateBoundarySegment(name,            // internal name of the boundary segment
                                                  left,             //  id of left subdomain
                                                  right,             //  id of right subdomain
                                                  index,         // Index of the segment
                                                  point,
                                                  alpha,
                                                  beta,
                                                  boundarySegmentFunction,
                                                  userData);
    }

    static void* CreateLinearSegment(const char *name,
                                     int left, int right,
                                     int index, int numVertices,
                                     const UG::INT* cornerIndices,
                                     double cornerCoordinates[2][ UG_DIM ])
    {
      return UG_NAMESPACE ::CreateLinearSegment(name,            // internal name of the boundary segment
                                                left,            //  id of left subdomain
                                                right,           //  id of right subdomain
                                                index,           // Index of the segment
                                                numVertices,
                                                cornerIndices,
                                                cornerCoordinates);
    }

  };

  template <>
  class UG_NS< UG_DIM >::Entity<0>
  {
  public:
    typedef UG_NAMESPACE ::element T;
  };

  template <>
  class UG_NS< UG_DIM >::Entity< UG_DIM - 1 >
  {
  public:
    typedef UG_NAMESPACE ::edge T;
  };

#if UG_DIM == 3
  template <>
  class UG_NS< UG_DIM >::Entity<1>
  {
  public:
    typedef UG_NAMESPACE ::vector T;
  };
#endif

  template <>
  class UG_NS< UG_DIM >::Entity< UG_DIM > {
  public:
    typedef UG_NAMESPACE ::node T;
  };

#undef UG_NAMESPACE

} // namespace Dune
