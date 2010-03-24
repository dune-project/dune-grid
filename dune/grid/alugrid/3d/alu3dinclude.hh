// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DINCLUDE_HH
#define DUNE_ALU3DINCLUDE_HH

//////////////////////////////////////////////////////////////////////
// compile imp.cc into lib (1 yes, 0 no)
// if you change this, you'll get what you deserve
//////////////////////////////////////////////////////////////////////
#define COMPILE_ALUGRID_LIB 0

#if COMPILE_ALUGRID_LIB
  #define COMPILE_ALUGRID_INLINE 0
#else
  #define COMPILE_ALUGRID_INLINE 1
#endif

#if COMPILE_ALUGRID_INLINE
#define alu_inline inline
#else
#define alu_inline
#endif
/////////////////////////////////////////////////////////////////////

// all methods and classes of the ALUGrid are defined in the namespace
#define ALU3DSPACE ALUGridSpace ::

#include <alugrid_defineparallel.h>

#if HAVE_MPI
// if this variable is defined,
// // then parallel version of ALUGrid is compiled
  #if ALU3DGRID_BUILD_FOR_PARALLEL == 0
    #warning "The ALUGrid-library wasn't compiled for parallel usage. Reconfigure\
  using the MPI compiler script or compile Dune without the MPI support!\
  Defaulting to serial ALUGrid!"
    #define ALU3DGRID_PARALLEL 0
  #else
    #define ALU3DGRID_PARALLEL 1
  #endif
#else
  #define ALU3DGRID_PARALLEL 0
#endif

// if MPI was found include all headers
#if ALU3DGRID_PARALLEL
#include <alugrid_parallel.h>
#else
// if not, include only headers for serial version
#include <alugrid_serial.h>
#endif

namespace ALUGridSpace {

#if ALU3DGRID_PARALLEL
  typedef GitterDunePll GitterType;
  //typedef Gitter GitterType;
  typedef GitterDunePll GitterImplType;

  typedef Hbnd3PllInternal<GitterImplType::Objects::Hbnd3Default,
      BndsegPllBaseXClosure<GitterImplType::hbndseg3_GEO>,
      BndsegPllBaseXMacroClosure<GitterImplType::hbndseg3_GEO> > :: micro_t MicroType;
#else
  // the header
  typedef Gitter GitterType;
  typedef GitterDuneImpl GitterImplType;
#endif

  // value for boundary to other processes
  static const int ProcessorBoundary_t = GitterType::hbndseg_STI::closure;

  // general GatherScatter type
  typedef GatherScatter GatherScatterType;

  // typedefs of Element types
  typedef GitterType::helement_STI HElementType;               // Interface Element
  typedef GitterType::hface_STI HFaceType;                     // Interface Face
  typedef GitterType::hedge_STI HEdgeType;                     // Interface Edge
  typedef GitterType::vertex_STI VertexType;                   // Interface Vertex
  typedef GitterType::hbndseg_STI HBndSegType;
  typedef GitterType::ghostpair_STI GhostPairType;
  typedef GitterType::Geometric::hface3_GEO GEOFace3Type;     // Tetra Face
  typedef GitterType::Geometric::hface4_GEO GEOFace4Type; // Hexa Face
  typedef GitterType::Geometric::hedge1_GEO GEOEdgeT;     // * stays real Face
  typedef GitterType::Geometric::VertexGeo GEOVertexT;     // * stays real Face
  typedef GitterImplType::Objects::tetra_IMPL IMPLTetraElementType; //impl Element
  typedef GitterImplType::Objects::hexa_IMPL IMPLHexaElementType;
  typedef GitterType::Geometric::tetra_GEO GEOTetraElementType;  // real Element
  typedef GitterType::Geometric::hexa_GEO GEOHexaElementType;
  typedef GitterType::Geometric::hasFace3 HasFace3Type;    // has Face with 3 polygons
  typedef GitterType::Geometric::hasFace4 HasFace4Type;
  typedef GitterType::Geometric::Hface3Rule Hface3RuleType;
  typedef GitterType::Geometric::Hface4Rule Hface4RuleType;

  typedef GitterImplType::Objects::Hbnd3Default BNDFace3Type;    // boundary segment
  typedef GitterImplType::Objects::Hbnd4Default BNDFace4Type;
  typedef GitterImplType::Objects::hbndseg3_IMPL ImplBndFace3Type;    // boundary segment
  typedef GitterImplType::Objects::hbndseg4_IMPL ImplBndFace4Type;
} // end namespace ALUGridSpace

//- local includes
#include "topology.hh"

namespace Dune {

  // typedef of ALU3dGridElementType see topology.hh

  // i.e. double or float
  typedef double alu3d_ctype;

  template <ALU3dGridElementType elType>
  struct ALU3dImplTraits {};

  template <>
  struct ALU3dImplTraits<tetra> {
    typedef ALU3DSPACE GEOFace3Type GEOFaceType;
    typedef ALU3DSPACE GEOEdgeT GEOEdgeType;
    typedef ALU3DSPACE GEOVertexT GEOVertexType;
    typedef ALU3DSPACE IMPLTetraElementType IMPLElementType;
    typedef ALU3DSPACE GEOTetraElementType GEOElementType;
    typedef ALU3DSPACE HasFace3Type HasFaceType;
    typedef ALU3DSPACE Hface3RuleType HfaceRuleType;
    typedef ALU3DSPACE BNDFace3Type BNDFaceType;
    typedef ALU3DSPACE ImplBndFace3Type ImplBndFaceType;
    typedef ALU3DSPACE BNDFace3Type PLLBndFaceType;
    typedef ALU3DSPACE HElementType HElementType;    // Interface Element
    typedef ALU3DSPACE HFaceType HFaceType;          // Interface Face
    typedef ALU3DSPACE HEdgeType HEdgeType;          // Interface Edge
    typedef ALU3DSPACE VertexType VertexType;        // Interface Vertex
    typedef ALU3DSPACE HBndSegType HBndSegType;      // Interface Boundaries

    // refinement and coarsening enum for tetrahedons
    enum { refine_element_t =
             ALU3DSPACE GitterType::Geometric::TetraRule::iso8 };
    enum { coarse_element_t =
             ALU3DSPACE GitterType::Geometric::TetraRule::crs  };
    enum { nosplit_element_t = ALU3DSPACE GitterType::Geometric::TetraRule::nosplit };

    typedef ALU3DSPACE GitterType::Geometric::TetraRule MarkRuleType;

    typedef std::pair<GEOFaceType*, int> NeighbourFaceType;
    typedef std::pair<HasFaceType*, int> NeighbourPairType;
    typedef ALU3DSPACE GhostPairType GhostPairType;

    template <int cdim>
    struct Codim;

  };

  template <>
  struct ALU3dImplTraits<tetra>::Codim<0> {
    typedef ALU3DSPACE GitterType::helement_STI InterfaceType;
    typedef IMPLElementType ImplementationType;
    typedef ALU3DSPACE HBndSegType GhostInterfaceType;
    typedef PLLBndFaceType GhostImplementationType;
  };

  template <>
  struct ALU3dImplTraits<tetra>::Codim<1> {
    typedef ALU3DSPACE GitterType::hface_STI InterfaceType;
    typedef GEOFaceType ImplementationType;
  };

  template <>
  struct ALU3dImplTraits<tetra>::Codim<2> {
    typedef ALU3DSPACE GitterType::hedge_STI InterfaceType;
    typedef GEOEdgeType ImplementationType;
  };

  template <>
  struct ALU3dImplTraits<tetra>::Codim<3> {
    typedef ALU3DSPACE GitterType::vertex_STI InterfaceType;
    typedef ALU3DSPACE GitterType::Geometric::VertexGeo ImplementationType;
  };



  template <>
  struct ALU3dImplTraits<hexa> {
    typedef ALU3DSPACE GEOFace4Type GEOFaceType;
    typedef ALU3DSPACE GEOEdgeT GEOEdgeType;
    typedef ALU3DSPACE GEOVertexT GEOVertexType;
    typedef ALU3DSPACE IMPLHexaElementType IMPLElementType;
    typedef ALU3DSPACE GEOHexaElementType GEOElementType;
    typedef ALU3DSPACE HasFace4Type HasFaceType;
    typedef ALU3DSPACE Hface4RuleType HfaceRuleType;
    typedef ALU3DSPACE BNDFace4Type BNDFaceType;
    typedef ALU3DSPACE ImplBndFace4Type ImplBndFaceType;
    typedef ALU3DSPACE BNDFace4Type PLLBndFaceType;
    typedef ALU3DSPACE HElementType HElementType;    // Interface Element
    typedef ALU3DSPACE HFaceType HFaceType;          // Interface Face
    typedef ALU3DSPACE HEdgeType HEdgeType;          // Interface Edge
    typedef ALU3DSPACE VertexType VertexType;        // Interface Vertex
    typedef ALU3DSPACE HBndSegType HBndSegType;      // Interface Boundaries

    // refinement and coarsening enum for hexahedrons
    enum { refine_element_t  = ALU3DSPACE GitterType::Geometric::HexaRule::iso8 };
    enum { coarse_element_t  = ALU3DSPACE GitterType::Geometric::HexaRule::crs  };
    enum { nosplit_element_t = ALU3DSPACE GitterType::Geometric::HexaRule::nosplit };

    typedef ALU3DSPACE GitterType::Geometric::HexaRule MarkRuleType;
    typedef std::pair<GEOFaceType*, int> NeighbourFaceType;
    typedef std::pair<HasFaceType*, int> NeighbourPairType;
    typedef ALU3DSPACE GhostPairType GhostPairType;

    template <int cdim>
    struct Codim;
  };

  template <>
  struct ALU3dImplTraits<hexa>::Codim<0> {
    typedef ALU3DSPACE GitterType::helement_STI InterfaceType;
    typedef IMPLElementType ImplementationType;
    typedef ALU3DSPACE HBndSegType GhostInterfaceType;
    typedef PLLBndFaceType GhostImplementationType;
  };

  template <>
  struct ALU3dImplTraits<hexa>::Codim<1> {
    typedef ALU3DSPACE GitterType::hface_STI InterfaceType;
    typedef GEOFaceType ImplementationType;
  };

  template <>
  struct ALU3dImplTraits<hexa>::Codim<2> {
    typedef ALU3DSPACE GitterType::hedge_STI InterfaceType;
    typedef GEOEdgeType ImplementationType;
  };

  template <>
  struct ALU3dImplTraits<hexa>::Codim<3> {
    typedef ALU3DSPACE GitterType::vertex_STI InterfaceType;
    typedef ALU3DSPACE GitterType::Geometric::VertexGeo ImplementationType;
  };

  //! contains list of vertices of one level
  //! needed for VertexLevelIterator
  class ALU3dGridVertexList
  {
  public:
    // level vertex iterator list
    typedef std::vector < ALU3DSPACE VertexType * > VertexListType;
    typedef VertexListType :: iterator IteratorType;

    ALU3dGridVertexList () : up2Date_(false) {}

    size_t size () const { return vertexList_.size(); }

    bool up2Date () const { return up2Date_;  }
    void unsetUp2Date ()  { up2Date_ = false; }

    // make grid walkthrough and calc global size
    template <class GridType>
    void setupVxList (const GridType & grid, int level);

    IteratorType begin () { return vertexList_.begin(); }
    IteratorType end   () { return vertexList_.end(); }

    VertexListType & getItemList() { return vertexList_; }
  private:
    bool up2Date_;
    VertexListType vertexList_;
  };

  //! contains list of vertices of one level
  //! needed for VertexLevelIterator
  class ALU3dGridLeafVertexList
  {
  public:
    // level vertex iterator list
    typedef std::pair < ALU3DSPACE VertexType * , int > ItemType;
    typedef std::vector < ItemType > VertexListType;
    typedef VertexListType :: iterator IteratorType;

    ALU3dGridLeafVertexList () : up2Date_(false) {}

    size_t size () const { return vertexList_.size(); }

    bool up2Date () const { return up2Date_;  }
    void unsetUp2Date ()  { up2Date_ = false; }

    // make grid walkthrough and calc global size
    template <class GridType>
    void setupVxList (const GridType & grid);

    IteratorType begin () { return vertexList_.begin(); }
    IteratorType end   () { return vertexList_.end(); }

    VertexListType & getItemList() { return vertexList_; }
    int getLevel(const ALU3DSPACE VertexType & vertex) const
    {
      const int idx = vertex.getIndex();
      assert( idx >= 0 );
      assert( idx < (int)size());
      const ItemType & p = vertexList_[idx];
      if( p.first == 0 )
        return vertex.level();
      else
        return p.second;
    }
  private:
    bool up2Date_;
    VertexListType vertexList_;
  };
  typedef ALU3dGridLeafVertexList LeafVertexListType;

  class ALU3dGridItemList
  {
  public:
    // level vertex iterator list
    typedef std::vector < void * > ItemListType;
    typedef ItemListType :: iterator IteratorType;

    ALU3dGridItemList () : up2Date_(false) {}

    size_t size () const { return itemList_.size(); }

    bool up2Date () const { return up2Date_;  }
    void unsetUp2Date ()  { up2Date_ = false; }

    void markAsUp2Date() { up2Date_ = true; }

    IteratorType begin () { return itemList_.begin(); }
    IteratorType end   () { return itemList_.end(); }

    ItemListType & getItemList() { return itemList_; }

  private:
    bool up2Date_;
    ItemListType itemList_;
  };

  typedef ALU3dGridItemList ALU3dGridItemListType;

  /////////////////////////////////////////////////////////////////////////
  //  some helper functions
  /////////////////////////////////////////////////////////////////////////

  inline const ALU3dImplTraits<tetra>::GEOFaceType*
  getFace(const ALU3DSPACE GEOTetraElementType& elem, int index) {
    assert(index >= 0 && index < 4);
    return elem.myhface3(ElementTopologyMapping<tetra>::dune2aluFace(index));
  }

  inline const ALU3dImplTraits<hexa>::GEOFaceType*
  getFace(const ALU3DSPACE GEOHexaElementType& elem, int index) {
    assert(index >= 0 && index < 6);
    return elem.myhface4(ElementTopologyMapping<hexa>::dune2aluFace(index));
  }

} // end namespace Dune

#endif
