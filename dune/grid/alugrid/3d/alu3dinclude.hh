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

#include <dune/common/mpicollectivecommunication.hh>
#include <dune/grid/alugrid/checkparallel.hh>

// if MPI was found include all headers
#if ALU3DGRID_PARALLEL
#include <alugrid_parallel.h>
#else
// if not, include only headers for serial version
#include <alugrid_serial.h>
#endif

//- local includes
#include <dune/grid/alugrid/3d/topology.hh>

namespace ALUGridSpace
{

  static const int ProcessorBoundary_t = Gitter::hbndseg_STI::closure;

  // general GatherScatter type
  typedef GatherScatter GatherScatterType;

}

namespace Dune
{

  // typedef of ALU3dGridElementType see topology.hh

  // i.e. double or float
  typedef double alu3d_ctype;



  // ALU3dBasicImplTraits
  // --------------------

  template< class Comm >
  struct ALU3dBasicImplTraits;

  template<>
  struct ALU3dBasicImplTraits< No_Comm >
  {
    typedef ALU3DSPACE Gitter GitterType;
    typedef ALU3DSPACE GitterDuneImpl GitterImplType;

    typedef GitterType::helement_STI HElementType;    // Interface Element
    typedef GitterType::hface_STI HFaceType;          // Interface Face
    typedef GitterType::hedge_STI HEdgeType;          // Interface Edge
    typedef GitterType::vertex_STI VertexType;        // Interface Vertex
    typedef GitterType::hbndseg_STI HBndSegType;
    typedef GitterType::ghostpair_STI GhostPairType;

    typedef HElementType PllElementType;

    typedef GitterType::Geometric::hedge1_GEO GEOEdgeType;

    //! method for ghost check
    template <class BndFaceType>
    static bool isGhost( const BndFaceType* ghost )
    {
      return false ;
    }
  };

#if ALU3DGRID_PARALLEL
  template<>
  struct ALU3dBasicImplTraits< MPI_Comm >
  {
    typedef ALU3DSPACE GitterDunePll GitterType;
    typedef ALU3DSPACE GitterDunePll GitterImplType;

    typedef GitterType::helement_STI HElementType;    // Interface Element
    typedef GitterType::hface_STI HFaceType;          // Interface Face
    typedef GitterType::hedge_STI HEdgeType;          // Interface Edge
    typedef GitterType::vertex_STI VertexType;        // Interface Vertex
    typedef GitterType::hbndseg_STI HBndSegType;
    typedef GitterType::ghostpair_STI GhostPairType;

    typedef ALU3DSPACE ElementPllXIF_t PllElementType;

    typedef GitterType::Geometric::hedge1_GEO GEOEdgeType;

    // method for ghost check
    template <class BndFaceType>
    static bool isGhost( const BndFaceType* ghost )
    {
      return ( ghost != 0 );
    }
  };
#endif // #if ALU3DGRID_PARALLEL



  // ALU3dCodimImplTraits
  // --------------------

  template< ALU3dGridElementType elType, class Comm, int codim >
  struct ALU3dCodimImplTraits;

  template< class Comm >
  struct ALU3dCodimImplTraits< tetra, Comm, 0 >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;
    typedef typename ALU3dBasicImplTraits< Comm >::GitterImplType GitterImplType;

    typedef typename GitterType::helement_STI InterfaceType;
    typedef typename GitterImplType::Objects::tetra_IMPL ImplementationType;
    typedef typename GitterType::hbndseg_STI GhostInterfaceType;
    typedef typename GitterImplType::Objects::Hbnd3Default GhostImplementationType;
  };

  template< class Comm >
  struct ALU3dCodimImplTraits< hexa, Comm, 0 >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;
    typedef typename ALU3dBasicImplTraits< Comm >::GitterImplType GitterImplType;

    typedef typename GitterType::helement_STI InterfaceType;
    typedef typename GitterImplType::Objects::hexa_IMPL ImplementationType;
    typedef typename GitterType::hbndseg_STI GhostInterfaceType;
    typedef typename GitterImplType::Objects::Hbnd4Default GhostImplementationType;
  };

  template< class Comm >
  struct ALU3dCodimImplTraits< tetra, Comm, 1 >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;

    typedef typename GitterType::hface_STI InterfaceType;
    typedef typename GitterType::Geometric::hface3_GEO ImplementationType;
  };

  template< class Comm >
  struct ALU3dCodimImplTraits< hexa, Comm, 1 >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;

    typedef typename GitterType::hface_STI InterfaceType;
    typedef typename GitterType::Geometric::hface4_GEO ImplementationType;
  };

  template< ALU3dGridElementType elType, class Comm >
  struct ALU3dCodimImplTraits< elType, Comm, 2 >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;

    typedef typename GitterType::hedge_STI InterfaceType;
    typedef typename GitterType::Geometric::hedge1_GEO ImplementationType;
  };

  template< ALU3dGridElementType elType, class Comm >
  struct ALU3dCodimImplTraits< elType, Comm, 3 >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;

    typedef typename GitterType::vertex_STI InterfaceType;
    typedef typename GitterType::Geometric::VertexGeo ImplementationType;
  };



  // ALU3dImplTraits
  // ---------------

  template< ALU3dGridElementType elType, class Comm >
  struct ALU3dImplTraits;

  template< class Comm >
  struct ALU3dImplTraits< tetra, Comm >
    : public ALU3dBasicImplTraits< Comm >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;
    typedef typename ALU3dBasicImplTraits< Comm >::GitterImplType GitterImplType;

    typedef typename GitterType::Geometric::hface3_GEO GEOFaceType;
    typedef typename GitterType::Geometric::VertexGeo GEOVertexType;
    typedef typename GitterImplType::Objects::tetra_IMPL IMPLElementType;
    typedef typename GitterType::Geometric::tetra_GEO GEOElementType;
    typedef typename GitterType::Geometric::periodic3_GEO GEOPeriodicType;
    typedef typename GitterType::Geometric::hasFace3 HasFaceType;
    typedef typename GitterType::Geometric::Hface3Rule HfaceRuleType;
    typedef typename GitterImplType::Objects::Hbnd3Default BNDFaceType;
    typedef typename GitterImplType::Objects::hbndseg3_IMPL ImplBndFaceType;

    typedef typename GitterType::Geometric::TetraRule MarkRuleType;

    // refinement and coarsening enum
    enum { refine_element_t = MarkRuleType::iso8 };
    enum { coarse_element_t = MarkRuleType::crs };
    enum { nosplit_element_t = MarkRuleType::nosplit };

    typedef std::pair< GEOFaceType *, int > NeighbourFaceType;
    typedef std::pair< HasFaceType *, int > NeighbourPairType;

    template< int codim >
    struct Codim
      : public ALU3dCodimImplTraits< tetra, Comm, codim >
    {};

    // access of faces
    template <class Elem>
    static const GEOFaceType* getFace( const Elem& elem, const int aluFace )
    {
      return elem.myhface3( aluFace );
    }
  };

  template< class Comm >
  struct ALU3dImplTraits< hexa, Comm >
    : public ALU3dBasicImplTraits< Comm >
  {
    typedef typename ALU3dBasicImplTraits< Comm >::GitterType GitterType;
    typedef typename ALU3dBasicImplTraits< Comm >::GitterImplType GitterImplType;

    typedef typename GitterType::Geometric::hface4_GEO GEOFaceType;
    typedef typename GitterType::Geometric::VertexGeo GEOVertexType;
    typedef typename GitterImplType::Objects::hexa_IMPL IMPLElementType;
    typedef typename GitterType::Geometric::hexa_GEO GEOElementType;
    typedef typename GitterType::Geometric::periodic4_GEO GEOPeriodicType;
    typedef typename GitterType::Geometric::hasFace4 HasFaceType;
    typedef typename GitterType::Geometric::Hface4Rule HfaceRuleType;
    typedef typename GitterImplType::Objects::Hbnd4Default BNDFaceType;
    typedef typename GitterImplType::Objects::hbndseg4_IMPL ImplBndFaceType;

    typedef typename GitterType::Geometric::HexaRule MarkRuleType;

    // refinement and coarsening enum
    enum { refine_element_t = MarkRuleType::iso8 };
    enum { coarse_element_t = MarkRuleType::crs };
    enum { nosplit_element_t = MarkRuleType::nosplit };

    typedef std::pair< GEOFaceType *, int > NeighbourFaceType;
    typedef std::pair< HasFaceType *, int > NeighbourPairType;

    template< int codim >
    struct Codim
      : public ALU3dCodimImplTraits< hexa, Comm, codim >
    {};

    // access of faces
    template <class Elem>
    static const GEOFaceType* getFace( const Elem& elem, const int aluFace )
    {
      return elem.myhface4( aluFace );
    }
  };



  //! contains list of vertices of one level
  //! needed for VertexLevelIterator
  template< class Comm >
  struct ALU3dGridVertexList
  {
    // level vertex iterator list
    typedef typename ALU3dBasicImplTraits< Comm >::VertexType VertexType;
    typedef std::vector< VertexType * > VertexListType;
    typedef typename VertexListType::iterator IteratorType;

    ALU3dGridVertexList ()
      : up2Date_( false )
    {}

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
  template< class Comm >
  struct ALU3dGridLeafVertexList
  {
    // level vertex iterator list
    typedef typename ALU3dBasicImplTraits< Comm >::VertexType VertexType;
    typedef std::pair< VertexType *, int > ItemType;
    typedef std::vector< ItemType > VertexListType;
    typedef typename VertexListType::iterator IteratorType;

    ALU3dGridLeafVertexList ()
      : up2Date_( false )
    {}

    size_t size () const { return vertexList_.size(); }

    bool up2Date () const { return up2Date_;  }
    void unsetUp2Date ()  { up2Date_ = false; }

    // make grid walkthrough and calc global size
    template <class GridType>
    void setupVxList (const GridType & grid);

    IteratorType begin () { return vertexList_.begin(); }
    IteratorType end   () { return vertexList_.end(); }

    VertexListType & getItemList() { return vertexList_; }

    int getLevel ( const VertexType &vertex ) const
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

  template< class Comm >
  struct ALU3dGridFaceGetter
  {
    static const typename ALU3dImplTraits< tetra, Comm >::GEOFaceType *
    getFace( const typename ALU3dImplTraits< tetra, Comm >::GEOElementType& elem, int index)
    {
      assert(index >= 0 && index < 4);
      return elem.myhface3( ElementTopologyMapping< tetra >::dune2aluFace(index) );
    }

    static const typename ALU3dImplTraits< hexa, Comm >::GEOFaceType*
    getFace( const typename ALU3dImplTraits< hexa, Comm >::GEOElementType &elem, int index )
    {
      assert(index >= 0 && index < 6);
      return elem.myhface4( ElementTopologyMapping< hexa >::dune2aluFace(index) );
    }
  };

} // end namespace Dune

#endif // #ifndef DUNE_ALU3DINCLUDE_HH
