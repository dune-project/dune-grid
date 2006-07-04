// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INCLUDE_HH
#define DUNE_ALU2DGRID_INCLUDE_HH

#include <alugrid_2d.h>
#define ALU2DSPACE ALU2dGridSpace ::

namespace Dune {

  struct ALU2dImplTraits {
    template <int cdim>
    struct Codim;
  };

  template<>
  struct ALU2dImplTraits::Codim<0> {
    typedef ALU2DSPACE Hmesh_basic::helement_t InterfaceType;
  };

  template<>
  struct ALU2dImplTraits::Codim<1> {
    typedef ALU2DSPACE Hmesh_basic::helement_t InterfaceType;
  };

  template <>
  struct ALU2dImplTraits::Codim<2> {
    typedef ALU2DSPACE Vertex InterfaceType;
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
      typedef typename Dune::ALU2dImplTraits::template Codim<0>::InterfaceType ElementType;
      typedef typename Dune::ALU2dImplTraits::template Codim<2>::InterfaceType VertexType;
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
        for(int i=0; i<dim+1; ++i)
        {
          enum { vxCodim = 1 };
          int vxIdx = elem.vertex(i)->getIndex();
          if( marker_[vxCodim][vxIdx] < 0) marker_[vxCodim][vxIdx] = elIdx;

          enum { edgeCodim = 0 };
          int edgeIdx = elem.edge_idx(i);
          if( marker_[edgeCodim][edgeIdx] < 0) marker_[edgeCodim][edgeIdx] = elIdx;
        }
      }

      //for(size_t k=0; k<marker_[0].size(); ++k)
      //  std::cout << "Edge[" << k <<"] is watched on el " << marker_[0][k] << "\n";
      up2Date_ = true;
    }

  private:
    VectorType marker_[2];

    bool up2Date_;
  };

} //end namespace Dune
#endif
