// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_BOUNDARYDOMBLOCK_HH
#define DUNE_DGF_BOUNDARYDOMBLOCK_HH

#include <iostream>
#include <vector>

#include <dune/grid/io/file/dgfparser/dgfbasicblock.hh>


namespace Dune
{

  namespace dgf
  {
    class BoundaryDomBlock
      : public BasicBlock
    {
      int dimworld;        // the dimension of the vertices (is given  from user)
      bool goodline;       // active line describes a vertex
      std::vector<double> p1,p2;    // active vertex
      int bndid;
      bool withdefault;
      int defaultvalue;
    public:
      // initialize vertex block and get first vertex
      BoundaryDomBlock(std::istream& in,int cdimworld );
      bool next ();
      bool inside(const std::vector<double>& v) const;
      int id() const {
        return bndid;
      }
      bool defaultValueGiven() {
        return withdefault;
      }
      int defaultValue() {
        return defaultvalue;
      }
      // some information
      bool ok() {
        return goodline;
      }
      int nofdombound() {
        return noflines();
      }
    };

  } // end namespace dgf

} // end namespace Dune

#endif
