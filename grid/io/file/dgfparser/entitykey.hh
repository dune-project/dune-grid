// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef MACROGRIDENTITYKEY_HH
#define MACROGRIDENTITYKEY_HH

#include <vector>
#include <algorithm>
#include <dune/grid/alugrid/3d/topology.hh>
namespace Dune {

  template < class A> class DGFEntityKey {
    inline DGFEntityKey () ;
    std::vector<A> key_,origKey_;
    bool origKeySet_;
  public:
    inline DGFEntityKey (const DGFEntityKey < A > &k) : key_(k.key_.size()), origKey_(k.key_.size()), origKeySet_(k. origKeySet_) {
      for (size_t i=0; i<key_.size(); i++) {
        key_[i]=k.key_[i];
        origKey_[i]=k.origKey_[i];
      }
    }
    inline DGFEntityKey& operator=(const DGFEntityKey < A > &k)
    {
      assert(key_.size()==k.key_.size());
      for (size_t i=0; i<key_.size(); i++) {
        key_[i]=k.key_[i];
        origKey_[i]=k.origKey_[i];
      }
      origKeySet_ = k.origKeySet_;
      return *this;
    }
    inline DGFEntityKey (std::vector<A>& key,bool setOrigKey=true) : key_(key.size()), origKey_(key.size()), origKeySet_(setOrigKey) {
      for (size_t i=0; i<key_.size(); i++) {
        key_[i]=key[i];
        origKey_[i]=key_[i];
      }
      std::sort(key_.begin(),key_.end());
    }
    inline DGFEntityKey (std::vector<A>& key,int N,int offset,bool setOrigKey=true) : key_(N), origKey_(N), origKeySet_(setOrigKey) {
      for (size_t i=0; i<key_.size(); i++) {
        key_[i]=key[(i+offset)%key.size()];
        origKey_[i]=key[(i+offset)%key.size()];
      }
      std::sort(key_.begin(),key_.end());
    }
    inline void orientation (int base,std::vector<std::vector<double> >& vtx) {
      if (key_.size()==3)  {
        assert( (size_t) origKey_[0] < vtx.size() );
        std::vector<double>& p0 = vtx[origKey_[0]];
        assert( (size_t) origKey_[1] < vtx.size() );
        std::vector<double>& p1 = vtx[origKey_[1]];
        assert( (size_t) origKey_[2] < vtx.size() );
        std::vector<double>& p2 = vtx[origKey_[2]];
        assert( (size_t) base < vtx.size() );
        std::vector<double>& q  = vtx[base];
        double n[3];
        n[0] = (p1[1]-p0[1])*(p2[2]-p0[2])-(p2[1]-p0[1])*(p1[2]-p0[2]);
        n[1] = (p1[2]-p0[2])*(p2[0]-p0[0])-(p2[2]-p0[2])*(p1[0]-p0[0]);
        n[2] = (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]);
        double test = n[0]*(q[0]-p0[0])+n[1]*(q[1]-p0[1])+n[2]*(q[2]-p0[2]);
        bool reorient = (test>0);
        if (reorient) {
          A key1=origKey_[1];
          origKey_[1]=origKey_[2];
          origKey_[2]=key1;
        }
      }
    }
    inline bool operator < (const DGFEntityKey <A> &k) const {
      // assert(k.key_.size()==key_.size());
      return key_<k.key_;
    }
    inline void print() const {
      for (size_t i=0; i<key_.size(); i++) {
        std::cerr << key_[i] << " ";
      }
      std::cerr << std::endl;
    }
    const A& operator[](int i) const {
      return key_[i];
    }
    bool origKeySet() const {
      return origKeySet_;
    }
    const A& origKey(int i) const {
      return origKey_[i];
    }
    int size() const {
      return key_.size();
    }
  } ;
  class ElementFaceUtil {
  public:
    inline static int nofFaces(int dimw,std::vector<int>& element) {
      if (dimw==1)
        return 2;
      else if (dimw==2)
        switch (element.size()) {
        case 3 : return 3; break;
        case 4 : return 4; break;
        }
      else if (dimw==3)
        switch (element.size()) {
        case 4 : return 4; break;
        case 8 : return 6; break;
        }
      return -1;
    }
    inline static int faceSize(int dimw,bool simpl) {
      if (dimw==1)
        return 1;
      else if (dimw==2)
        return 2;
      else if (dimw==3)
        return ((simpl) ? 3 : 4);
      return -1;
    }
    template <int dimworld>
    inline static DGFEntityKey<int>
    generateCubeFace(std::vector<int>& element,int f) {
      ReferenceCube<double,dimworld> ref;
      int size=ref.size(f,1,dimworld);
      std::vector<int> k(size);
      /*
         for (int i=0;i<size;i++) {
         k[i] = element[ref.subEntity(f,1,i,dimworld)];
         }
         if (dimworld==3) {
         if (f==2 || f==1 || f==5) {
          int ktmp=k[0];
          k[0]=k[1];
          k[1]=ktmp;
         }
         else {
          int ktmp=k[2];
          k[2]=k[3];
          k[3]=ktmp;
         }
         }
       */
      int face=ElementTopologyMapping<hexa>::dune2aluFace(f);
      for (int i=0; i<size; i++) {
        // int idxdune = ref.subEntity(f,1,i,dimworld);
        int idx = ElementTopologyMapping<hexa>::alu2duneFaceVertex(face,i);
        int idxdune = ref.subEntity(f,1,idx,dimworld);
        k[size-1-i] = element[idxdune];
      }
      return DGFEntityKey<int> (k);
    }
    template <int dimworld>
    inline static DGFEntityKey<int>
    generateSimplexFace(std::vector<int>& element,int f) {
      ReferenceSimplex<double,dimworld> ref;
      int size=ref.size(f,1,dimworld);
      std::vector<int> k(size);
      for (int i=0; i<size; i++) {
        k[i] = element[ref.subEntity(f,1,i,dimworld)];
      }
      return DGFEntityKey<int> (k);
    }
    inline static DGFEntityKey<int>
    generateFace(int dimw,std::vector<int>& element,int f) {
      if (element.size()==size_t(dimw+1)) { // Simplex element
        if (dimw==3)
          return generateSimplexFace<3>(element,f);
        else if (dimw==2)
          return generateSimplexFace<2>(element,f);
        else if (dimw==1)
          return generateSimplexFace<1>(element,f);
      }
      else { // Cube element
        if (dimw==3)
          return generateCubeFace<3>(element,f);
        else if (dimw==2)
          return generateCubeFace<2>(element,f);
        else if (dimw==1)
          return generateCubeFace<1>(element,f);
      }
      DUNE_THROW(DGFException,"WRONG DIMENSION");
      return generateCubeFace<1>(element,f);
    }
  };

} //end namespace Dune
#endif
