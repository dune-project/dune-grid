// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INCLUDE_HH
#define DUNE_ALU2DGRID_INCLUDE_HH

#include <alugrid_2d.h>
#define ALU2DSPACE ALUGridSpace ::
#define ALU2DSPACENAME ALUGridSpace

#ifdef ALUGRID_SURFACE_2D
#define ALU2DDIMWORLD(dimw) < dimw >
#else
#define ALU2DDIMWORLD(dimw)
#endif

// use the ALU3dGrid Parallel detection
//#define ALU2DGRID_PARALLEL ALU3DGRID_PARALLEL
#define ALU2DGRID_PARALLEL 0

#include <dune/common/collectivecommunication.hh>

#if ALU2DGRID_PARALLEL
//#include "communicator.hh"
#warning "Using ALU2dGrid in parallel"
#endif

namespace Dune {

  template< int dim, int dimw >
  struct ALU2dImplInterface;

  template< int dimw >
  struct ALU2dImplInterface< 0, dimw >
  {
#ifdef ALUGRID_SURFACE_2D
    typedef ALU2DSPACE Hmesh_basic::vertex_t Type;
#else
    typedef ALU2DSPACE Vertex Type;
#endif
  };

  template< int dimw >
  struct ALU2dImplInterface< 1, dimw >
  {
    typedef ALU2DSPACE Hmesh_basic::helement_t Type;
  };

  template< int dimw >
  struct ALU2dImplInterface< 2, dimw >
  {
    typedef ALU2DSPACE Hmesh_basic::helement_t Type;
  };


  template< int dimw >
  struct ALU2dImplTraits
  {
    template< int cdim >
    struct Codim
    {
      typedef typename ALU2dImplInterface< 2-cdim, dimw >::Type InterfaceType;
    };
  };


  class ALU2dGridMarkerVector
  {
    typedef std::vector< int > VectorType;
  public:
    ALU2dGridMarkerVector() : up2Date_(false) {}

    bool up2Date() const { return up2Date_; }

    void unsetUp2Date() { up2Date_ = false; }

    bool isOnElement(int elementIndex, int idx, int codim) const
    {
      return marker_[codim-1][idx] == elementIndex;
    }

    template <class GridType>
    void update (const GridType & grid, int level )
    {
      static const int dimworld = GridType::dimensionworld;
      typedef typename Dune::ALU2dImplTraits< dimworld >::template Codim<0>::InterfaceType ElementType;
      typedef typename Dune::ALU2dImplTraits< dimworld >::template Codim<2>::InterfaceType VertexType;
      typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

      // resize
      for(int i=0; i<2; ++i)
      {
        int s = grid.hierSetSize(i+1);
        if((int) marker_[i].size() < s ) marker_[i].resize(s);

        size_t markerSize = marker_[i].size();
        // reset marker vector to default value
        for(size_t k=0; k<markerSize; ++k) marker_[i][k] = -1;
      }

      enum { dim = GridType::dimension };
      IteratorType iter(grid.myGrid(), level);

      for(iter->first(); !iter->done(); iter->next())
      {
        ElementType & elem = iter->getitem();
        int elIdx = elem.getIndex();

        // if element is not valid, go to next
#if ALU2DGRID_PARALLEL
        if( ! grid.rankManager().isValid( elIdx , All_Partition ) ) continue;
#endif
        for(int i=0; i<dim+1; ++i)
        {
          enum { vxCodim = 1 };
          int vxIdx = elem.getVertex(i)->getIndex();
          if( marker_[vxCodim][vxIdx] < 0) marker_[vxCodim][vxIdx] = elIdx;

          enum { edgeCodim = 0 };
          int edgeIdx = elem.edge_idx(i);
          if( marker_[edgeCodim][edgeIdx] < 0) marker_[edgeCodim][edgeIdx] = elIdx;
        }
      }
      up2Date_ = true;
    }

  private:
    VectorType marker_[2];

    bool up2Date_;
  };

  class ALU2dGridLeafMarkerVector
  {
    typedef std::vector< int > VectorType;
  public:
    ALU2dGridLeafMarkerVector() : up2Date_(false) {}

    bool up2Date() const { return up2Date_; }

    void unsetUp2Date() { up2Date_ = false; }

    // return true, if edge is visited on given element
    bool isOnElement(int elementIndex, int idx, int codim) const
    {
      assert( up2Date_ );
      // this marker only works for codim 1, i.e. edges
      assert( codim == 1 );
      return marker_[idx] == elementIndex;
    }

    // this is for the LeafIterator
    template <class GridType>
    void update (const GridType & grid)
    {
      static const int dimworld = GridType::dimensionworld;
      typedef typename Dune::ALU2dImplTraits< dimworld >::template Codim<0>::InterfaceType ElementType;
      typedef typename Dune::ALU2dImplTraits< dimworld >::template Codim<2>::InterfaceType VertexType;
      typedef ALU2DSPACE Listwalkptr< ElementType > IteratorType;

      // resize edge marker
      {
        int s = grid.hierSetSize(1);
        if((int) marker_.size() < s ) marker_.resize(s);

        size_t markerSize = marker_.size();
        // reset marker vector to default value
        for(size_t k=0; k<markerSize; ++k) marker_[k] = -1;
      }

      // resize vertex levels
      {
        int s = grid.hierSetSize(2);
        if((int) vertexLevels_.size() < s ) vertexLevels_.resize(s);

        // initialize with -1
        size_t vxSize = vertexLevels_.size();
        for(size_t k=0; k<vxSize; ++k) vertexLevels_[k] = -1;
      }

      enum { dim = GridType::dimension };
      IteratorType iter(grid.myGrid());

      for(iter->first(); !iter->done(); iter->next())
      {
        ElementType & elem = iter->getitem();
        int elIdx = elem.getIndex();

#if ALU2DGRID_PARALLEL
        // is element is not valid, go to next
        if( ! grid.rankManager().isValid( elIdx , All_Partition ) ) continue;
#endif
        int level = elem.level();

        for(int i=0; i<dim+1; ++i)
        {
          int vxIdx = elem.getVertex(i)->getIndex();

          // set max level to vertices, see Grid docu paper
          if(level > vertexLevels_[vxIdx]) vertexLevels_[vxIdx] = level;

          int edgeIdx = elem.edge_idx(i);
          if( marker_[edgeIdx] < 0) marker_[edgeIdx] = elIdx;
        }
      }
      up2Date_ = true;
    }

    //! return level of vertex
    int levelOfVertex(const int vxIdx) const
    {
      assert( up2Date_ );
      assert( vxIdx >= 0 && vxIdx < (int) vertexLevels_.size());
      // if this assertion is thrown, the level has not been initialized
      assert( vertexLevels_[vxIdx] >= 0 );
      return vertexLevels_[vxIdx];
    }

    //! return level of vertex
    bool isValidVertex(const int vxIdx) const
    {
      assert( up2Date_ );
      assert( vxIdx >= 0 && vxIdx < (int) vertexLevels_.size());
      return (vertexLevels_[vxIdx] >= 0);
    }

  private:
    VectorType marker_;
    VectorType vertexLevels_;

    bool up2Date_;
  };

  // dummy object stream class
  class ALU2dGridObjectStream
  {
  public:
    class EOFException {} ;
    template <class T>
    void readObject (T &) {}
    void readObject (int) {}
    void readObject (double) {}
    template <class T>
    void writeObject (T &) {}
    void writeObject (int) {}
    void writeObject (double) {}

    template <class T>
    void read (T &) const {}
    template <class T>
    void write (const T &) {}
  };

} //end namespace Dune
#endif
