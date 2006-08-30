// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef MACROGRIDENTITYKEY_HH
#define MACROGRIDENTITYKEY_HH

#include <vector>

namespace Dune {
  namespace {
    template < class A> class EntityKey {
      inline EntityKey () ;
      std::vector<A> key_,origKey_;
    public:
      inline EntityKey (const EntityKey < A > &k) : key_(k.key_.size()), origKey_(k.key_.size()) {
        for (size_t i=0; i<key_.size(); i++) {
          key_[i]=k.key_[i];
          origKey_[i]=k.origKey_[i];
        }
      }
      inline EntityKey (std::vector<A>& key) : key_(key.size()), origKey_(key.size()) {
        for (size_t i=0; i<key_.size(); i++) {
          key_[i]=key[i];
          origKey_[i]=key_[i];
        }
        std::sort(key_.begin(),key_.end());
      }
      inline EntityKey (std::vector<A>& key,int N,int offset) : key_(N), origKey_(N) {
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
      inline bool operator < (const EntityKey <A> &k) const {
        assert(k.key_.size()==key_.size());
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
      const A& origKey(int i) const {
        return origKey_[i];
      }
      int size() const {
        return key_.size();
      }
    } ;
  }
}
#endif
