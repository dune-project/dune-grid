// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  template <int dim,int dimworld>
  inline OneDGrid<dim,dimworld>*
  MacroGrid ::
  Impl<OneDGrid<dim,dimworld> >::generate(MacroGrid& mg,
                                          const char* filename, MPICommunicatorType ) {
    mg.element=Cube;
    std::ifstream gridin(filename);
    if(mg.readDuneGrid(gridin))
    {
      typedef std::map<size_t,double> VtxMapType;
      typedef typename VtxMapType :: iterator iterator;
      VtxMapType vtxlist;
      size_t size = mg.vtx.size();
      for(size_t i=0; i<size; ++i)
      {
        iterator vtxend = vtxlist.end();
        bool found = false;
        for(iterator it=vtxlist.begin(); it!= vtxend; ++it)
        {
          if( std::abs((*it).second - mg.vtx[i][0]) < 1e-10)
          {
            found = true;
            break;
          }
        }

        if(!found) vtxlist[i] = mg.vtx[i][0];
      }

      std::vector<double> vtx;
      vtx.reserve(vtxlist.size());
      iterator vtxend = vtxlist.end();
      for(iterator it=vtxlist.begin(); it!= vtxend; ++it)
      {
        vtx.push_back( (*it).second );
      }

      // sort vector, otherwise OneDGrid not satisfied
      std::sort(vtx.begin(), vtx.end());

      return new OneDGrid<dim,dimworld>(vtx);
    }
    DUNE_THROW(DGFException, "Unrecoverable Error in dgfpaser<OneDGrid>");
    return 0;
  }
}
