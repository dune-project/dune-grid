// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/io/file/dgfparser/dgfparserblocks.hh>

namespace Dune
{

  namespace dgf
  {

    // BasicBlock
    // ----------

    void BasicBlock :: getblock ( std :: istream &in )
    {
      std::string id;
      getline(in,id);
      while (in.good()) {
        std::stringstream idstream(id);
        std::string upcaseid;
        idstream >> upcaseid;
        makeupcase(upcaseid);
        if (upcaseid==identifier)
        {
          active=true;
          break;
        }
        getline(in,id);
      }
      if (active) {
        bool blockend=false;
        while (in.good()) {
          getline(in,oneline);
          if (oneline.size()==0)
            continue;
          std::stringstream onelinestream(oneline);
          std::string test;
          onelinestream >> test;
          if (test[0] == '#') {
            blockend=true;
            break;
          }
          empty=false;
          block << oneline << "\n";
        }
        if (!blockend) {
          DUNE_THROW(DGFException,
                     "Error: block must end with a #-line");
        }
      }
      else {
        // std::cerr << "Warning: Block " << identifier << " not found" << std::endl;
      }
    }


    int BasicBlock :: countlines ()
    {
      if (empty)
        return 0;
      int ret=0;
      while (1) {
        getnextline();
        if (oneline.size()==0)
          break;
        ret++;
      }
      return ret;
    }


    void BasicBlock :: reset ()
    {
      pos=-1;
      block.clear();
      block.seekg(0);
      linecount=countlines();
      pos=-1;
      block.clear();
      block.seekg(0);
    }


    // get next line and store in string stream
    void BasicBlock :: getnextline ()
    {
      line.clear();
      oneline.clear();
      getline(block,oneline);
      if (oneline.size()>0) {
        std::size_t comment=oneline.find("%");
        if (comment!=std::string::npos) {
          oneline.erase(comment);
          if (oneline.size()==0) {
            getnextline();
            return;
          }
        }
      }
      line.str(oneline);
      pos++;
    }


    bool BasicBlock :: gettokenparam ( std :: string token, std :: string &entry )
    {
      makeupcase(token);
      std::string ltoken;
      reset();
      do {
        getnextline();
        if (oneline.size()==0)
          return false;
        line >> ltoken;
        makeupcase(ltoken);
      } while (ltoken!=token);
      getline(line,entry);
      return true;
    }


    bool BasicBlock :: findtoken ( std :: string token )
    {
      makeupcase(token);
      std::string ltoken;
      reset();
      do {
        getnextline();
        if (oneline.size()==0)
          return false;
        line >> ltoken;
        makeupcase(ltoken);
      } while (ltoken!=token);
      return true;
    }


    BasicBlock :: BasicBlock ( std :: istream &in, const char* id )
      : pos(-1),
        active(false),
        empty(true),
        identifier(id),
        linecount(0)
    {
      makeupcase(identifier);
      in.clear();
      in.seekg(0);
      if (!in) {
        DUNE_THROW(DGFException,
                   "file not found in BasicBlock::BasicBlock");
      }
      getblock(in);
      if (active && !empty) {
        linecount=countlines();
        reset();
      }
      in.clear();
      in.seekg(0);
    }



    // VertexBlock
    // -----------

    const char *VertexBlock :: ID = "Vertex";


    VertexBlock :: VertexBlock ( std :: istream &in, int &pdimworld )
      : BasicBlock(in,ID),
        dimworld(pdimworld),
        goodline(true),
        p(0),
        vtxoffset(0),
        nofParam(0),
        vtxparam(0)
    {
      if (!isactive())
        return;
      if (dimworld<0)
        dimworld=0;
      {
        int x;
        if (findtoken("firstindex")) {
          if (getnextentry(x)) {
            vtxoffset=x;
          }
        }
      }
      {
        int x;
        if (findtoken("parameters")) {
          if (getnextentry(x)) {
            nofParam=x;
            vtxparam.resize(nofParam);
          }
        }
      }
      dimworld=getDimW();
      if (dimworld>0) {
        p.resize(dimworld);
      } else {
        DUNE_THROW(DGFException,
                   "ERROR in " << *this
                               << "      no line with enough entries found");
      }
      reset();
      next();
      pdimworld=dimworld;
    }


    int VertexBlock :: get ( std :: vector< std :: vector< double > > &vtx,
                             std :: vector< std :: vector< double > > &param,
                             int &nofp )
    {
      nofp=nofParam;
      size_t nofvtx;
      size_t old_size = vtx.size();
      // vtx.resize(old_size+nofvertex());
      for (nofvtx=old_size; ok(); next(),nofvtx++) {
        vtx.push_back(p);
        if (nofParam>0)
          param.push_back(vtxparam);
      }
      return nofvtx;
    }


    int VertexBlock :: getDimW ()
    {
      reset();
      int dimworld=0;
      getnextline();
      while (dimworld<1 && linenumber()<noflines()) {
        dimworld = 0;
        double x;
        while (getnextentry(x)) {
          dimworld++;
        }
        dimworld-=nofParam;
        getnextline();
      }
      return dimworld;
    }


    bool VertexBlock :: next ()
    {
      assert(ok());
      int n=0;
      getnextline();
      if (linenumber()==noflines()) {
        goodline=false;
        return goodline;
      }
      double x;
      while (getnextentry(x)) {
        if (n<dimworld)
          p[n]=x;
        else if (n-dimworld<nofParam) {
          vtxparam[n-dimworld]=x;
        }
        n++;
      }
      if (n>0 && n!=dimworld+nofParam) {
        return next();
        /*
           DUNE_THROW(DGFException,
                   "ERROR in " << *this
                   << "      wrong number of coordinates and parameters: "
                   << n << " read but expected " << dimworld+nofParam);
         */
      }
      if (n==0) {
        return next();
      }
      goodline=true;
      if (!goodline) {
        DUNE_THROW(DGFException,
                   "ERROR in " << *this
                               << "      wrong number of coordinates: "
                               << n << " read but expected " << dimworld);
      }
      return goodline;
    }



    // SimplexGenerationBlock
    // ----------------------

    const char *SimplexGenerationBlock :: ID = "Simplexgenerator";


    SimplexGenerationBlock :: SimplexGenerationBlock ( std :: istream &in )
      : BasicBlock(in,ID),
        area_(-1),
        angle_(-1),
        display_(false),
        haspath_(false),
        filetype_(),
        parameter_(),
        hasfile_(false),
        dimension_(-1)
    {
      double x;
      bool b;
      int i;
      std::string p;
      if (findtoken("max-area"))
        if (getnextentry(x))
          area_=x;
      if (findtoken("min-angle"))
        if (getnextentry(x))
          angle_=x;
      if (findtoken("display"))
        if (getnextentry(b))
          display_=b;
      if (findtoken("path"))
        if (getnextentry(p)) {
          path_=p;
          haspath_=true;
        }
      if (findtoken("file")) {
        if (getnextentry(p)) {
          filename_=p;
          hasfile_=true;
        }
        if (getnextentry(p)) {
          filetype_=p;
        }
        if (findtoken("dimension"))
          if (getnextentry(i)) {
            dimension_=i;
          }
        gettokenparam("parameter",parameter_);
      }
    }


    // SimplexBlock
    // ------------

    const char *SimplexBlock :: ID = "Simplex";


    SimplexBlock
    :: SimplexBlock ( std :: istream &in, int pnofvtx, int pvtxoffset, int adimworld )
      : BasicBlock(in,ID),
        nofvtx(pnofvtx),
        vtxoffset(pvtxoffset),
        dimworld(adimworld),
        goodline(true),
        p(adimworld+1),
        nofparams(0),
        psimpl(0)
    {
      if (!isactive()) return;
      assert((dimworld+1)>0);
      if (findtoken("parameters")) {
        int x=0;
        if (getnextentry(x)) {
          if (x>0) {
            nofparams = x;
            psimpl.resize(nofparams);
          }
        }
        if (x<=0) {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      parameter key found with no or non-positive value "
                                 << x);
        }
      }
      reset();
      next();
    }


    int SimplexBlock :: get ( std :: vector< std :: vector< unsigned int > > &simplex,
                              std :: vector< std :: vector< double > > &params,
                              int &nofp)
    {
      nofp=nofparams;
      int nofsimpl;
      for (nofsimpl=0; ok(); next(), nofsimpl++) {
        simplex.push_back(p);
        for (size_t j=0; j<p.size(); j++) {
          simplex[nofsimpl][j] = p[j];
        }
        if (nofparams>0) {
          params.push_back(psimpl);
        }
      }
      /*
         // make numbering starting from zero
         // not matter whether offset is positive or negative
         offset = vtxoffset;
         if(offset != 0) {
         for (int i=0; i<nofsimpl; ++i) {
          for (size_t j=0;j<simplex[i].size();++j) {
            simplex[i][j] -= offset;
          }
         }
         }
       */
      return nofsimpl;
    }


    int SimplexBlock
    :: cube2simplex ( std :: vector< std :: vector< double > > &vtx,
                      std :: vector< std :: vector< unsigned int > > &elements,
                      std :: vector< std :: vector< double > > &params )
    {
      static int offset3[6][4][3] = {{{0,0,0},{1,1,1},{1,0,0},{1,1,0}},
                                     {{0,0,0},{1,1,1},{1,0,1},{1,0,0}},
                                     {{0,0,0},{1,1,1},{0,0,1},{1,0,1}},
                                     {{0,0,0},{1,1,1},{1,1,0},{0,1,0}},
                                     {{0,0,0},{1,1,1},{0,1,0},{0,1,1}},
                                     {{0,0,0},{1,1,1},{0,1,1},{0,0,1}} };
      static int offset2[2][3][2] = {{{0,0},{1,0},{0,1}},
                                     {{1,1},{0,1},{1,0}}};
      if( vtx.size()==0 )
        DUNE_THROW( DGFException, "Cannot convert emty cube grid to simplex grid.");
      const int dimworld = vtx[ 0 ].size();

      dverb << "generating simplices...";
      dverb.flush();

      if( dimworld == 1 )
        return elements.size();

      std :: vector< std :: vector< unsigned int > > cubes;
      std :: vector< std :: vector< double > > cubeparams;
      elements.swap( cubes );
      params.swap( cubeparams );

      if(dimworld == 3)
      {
        elements.resize(6*cubes.size());
        if (cubeparams.size()>0)
          params.resize(6*cubes.size());
        for(size_t countsimpl=0; countsimpl < elements.size(); countsimpl++)
          elements[countsimpl].resize(4);
        for (size_t c=0; c<cubes.size(); c++)  {
          for(int tetra=0; tetra < 6 ; tetra++) {
            for (int v=0; v<4; v++) {
              elements[c*6+tetra][v]=
                cubes[c][offset3[tetra][v][0]+
                         offset3[tetra][v][1]*2+
                         offset3[tetra][v][2]*4];
            }
            if (cubeparams.size()>0)
              params[c*6+tetra] = cubeparams[c];
          }
        }
      }
      else if (dimworld==2) {
        elements.resize(2*cubes.size() );
        if (cubeparams.size()>0)
          params.resize(2*cubes.size());
        for(size_t countsimpl=0; countsimpl < elements.size(); countsimpl++)
          elements[countsimpl].resize(3);
        for (size_t c=0; c<cubes.size(); c++) {
          int diag = 0;
          double mind = 0;
          for (int d=0; d<2; d++) {
            double diaglen =
              pow(vtx[cubes[c][d]][0]-vtx[cubes[c][2+((d+1)%2)]][0],2) +
              pow(vtx[cubes[c][d]][1]-vtx[cubes[c][2+((d+1)%2)]][1],2);
            if (diaglen<mind) {
              mind=diaglen;
              diag = d;
            }
          }

          if (diag == 0) {
            int tmp0 = cubes[c][0];
            cubes[c][0] = cubes[c][1];
            cubes[c][1] = cubes[c][3];
            cubes[c][3] = cubes[c][2];
            cubes[c][2] = tmp0;
          }
          for(int tetra=0; tetra < 2 ; tetra++) {
            for (int v=0; v<3; v++) {
              elements[c*2+tetra][v]=
                cubes[c][offset2[tetra][v][0]+
                         offset2[tetra][v][1]*2];
            }
            if (cubeparams.size()>0)
              params[c*2+tetra] = cubeparams[c];
          }
        }
      }
      return elements.size();
    }


    bool SimplexBlock :: next ()
    {
      assert(ok());
      int n=0;
      getnextline();
      if (linenumber()==noflines()) {
        goodline=false;
        return goodline;
      }
      double x;
      while (getnextentry(x)) {
        if (n<(int)p.size()) {
          p[n]=int(x)-vtxoffset;
          if (p[n]<0 || p[n]>=nofvtx) {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      wrong index of vertices: "
                                   << x << " read but expected value between "
                                   << vtxoffset << " and "
                                   << nofvtx+vtxoffset);
          }
        }
        else if (n-int(p.size())<nofparams) {
          psimpl[n-p.size()]=x;
        }
        n++;
      }
      // tests if the written block is ok in its size
      if (n!=(int)p.size()+nofparams) {
        return next();
      }
      // tests if the written block is ok in its size
      goodline=(n==(int)p.size()+nofparams);
      return goodline;
    }



    // CubeBlock
    // ---------

    const char *CubeBlock :: ID = "Cube";

    CubeBlock :: CubeBlock(std::istream& in,int pnofvtx, int pvtxoffset, int adimworld)
      : BasicBlock(in,ID),
        nofvtx(pnofvtx),
        dimworld(adimworld),
        goodline(true),
        p(0),
        map(0),
        nofparams(0),
        vtxoffset(pvtxoffset),
        psimpl(0)
    {
      if (!isactive()) return;
      assert((dimworld+1)>0);
      p.resize(1<<dimworld);
      map.resize(1<<dimworld);
      int x;
      if (findtoken("map")) {
        for (size_t i=0; i<map.size(); i++) {
          if (getnextentry(x)) {
            map[i]=x;
          } else {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      reference maping not complete "
                                   << i
                                   << " entries read but expected "
                                   << map.size());
          }
        }
      } else {
        for (size_t i=0; i<map.size(); i++) {
          map[i]=i;
        }
      }
      if (findtoken("parameters")) {
        int x=0;
        if (getnextentry(x)) {
          if (x>0) {
            nofparams = x;
            psimpl.resize(nofparams);
          }
        }
        if (x<=0) {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      parameter key found with no or non-positive value "
                                 << x);
        }
      }
      reset();
      next();
    }

    int CubeBlock :: get ( std :: vector< std :: vector< unsigned int> > &simplex,
                           std :: vector< std :: vector< double > > &params,
                           int &nofp )
    {
      nofp=nofparams;
      int nofsimpl;
      // simplex.resize(nofsimplex());
      for (nofsimpl=0; ok(); next(), nofsimpl++) {
        simplex.push_back(p);
        for (size_t j=0; j<p.size(); j++) {
          simplex[nofsimpl][map[j]] = p[j];
        }
        if (nofparams>0) {
          params.push_back(psimpl);
        }
      }
      /*
         int offset = vtxoffset;
         if(offset != 0)
         {
         for (int i=0; i<nofsimpl; ++i)
         {
          for (size_t j=0;j<simplex[i].size();++j)
          {
            simplex[i][j] -= offset;
          }
         }
         }
       */
      return nofsimpl;
    }


    bool CubeBlock :: next ()
    {
      assert(ok());
      int n=0;
      getnextline();
      if (linenumber()==noflines()) {
        goodline=false;
        return goodline;
      }
      double x;
      while (getnextentry(x))
      {
        if (n<(int)p.size())
        {
          p[n]=int(x)-vtxoffset;
          if (p[n]<0 || p[n]>=nofvtx)
          {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      wrong index of vertices: "
                                   << x << " read but expected value between "
                                   << vtxoffset << " and "
                                   << nofvtx+vtxoffset);
          }
        } else if (n-int(p.size())<nofparams) {
          psimpl[n-p.size()]=x;
        }
        n++;
      }
      // tests if the written block is ok in its size
      if (n!=(int)p.size()+nofparams) {
        return next();
      }
      goodline=(n==(int)p.size()+nofparams);
      return goodline;
    }



    // BoundaryDomBlock
    // ----------------

    const char *BoundaryDomBlock :: ID = "boundarydomain";

    BoundaryDomBlock :: BoundaryDomBlock(std::istream& in,int cdimworld )
      : BasicBlock(in,ID),
        dimworld(cdimworld),
        goodline(true),
        p1(cdimworld),
        p2(cdimworld),
        bndid(0),
        withdefault(false),
        defaultvalue(0)
    {
      if (!isactive())
        return;
      assert(cdimworld>0);
      {
        int x;
        if (findtoken("default"))
        {
          if (getnextentry(x))
          {
            if( x <= 0 )
            {
              DUNE_THROW(DGFException,
                         "ERROR in " << *this
                                     << "      non-positive boundary id (" << x << ") read!");
            }
            defaultvalue=x;
            withdefault = true;
          }
        }
      }
      reset();
      next();
    }


    bool BoundaryDomBlock :: next ()
    {
      assert(ok());
      getnextline();
      if (linenumber()==noflines()) {
        goodline=false;
        return goodline;
      }
      int id;
      if (getnextentry(id))
      {
        if( id <= 0 )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      non-positive boundary id (" << id << ") read!");
        }
        bndid = id;
        double x;
        int n=0;
        while (getnextentry(x))
        {
          if (0<=n && n<dimworld)
            p1.at(n)=x;
          else if (dimworld<=n && n<2*dimworld)
          {
            p2.at(n-dimworld)=x;
            if (p2.at(n-dimworld)<p1.at(n-dimworld))
            {
              DUNE_THROW(DGFException,
                         "ERROR in " << *this
                                     << "      second coordinate smaller than first coordinate: "
                                     << p2.at(n-dimworld)
                                     << " read but expected value larger or equal to "
                                     << p1.at(n-dimworld));
            }
          }
          n++;
        }

        goodline=(n==dimworld*2);
        if (!goodline)
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      wrong number of coordinates: "
                                 << n << " read but expected " << dimworld);
        }
      }
      else
        next();
      return goodline;
    }


    bool BoundaryDomBlock :: inside ( const std :: vector< double > &v ) const
    {
      assert( v.size() == (size_t)dimworld );
      for( int i = 0; i < dimworld; ++i )
      {
        if( (v[ i ] < p1[ i ]) || (v[ i ] > p2[ i ]) )
          return false;
      }
      return true;
    }


    // BoundarySegBlock
    // ----------------

    const char *BoundarySegBlock :: ID = "boundarysegments";


    BoundarySegBlock :: BoundarySegBlock ( std :: istream &in, int pnofvtx,
                                           int pdimworld,bool psimplexgrid )
      : BasicBlock(in,ID),
        dimworld(pdimworld),
        goodline(true),
        p(),
        bndid(-1),
        simplexgrid(psimplexgrid)
    {
      if (!isactive())
        return;
      assert(dimworld>0);
      next();
    }


    int BoundarySegBlock
    :: get( std :: map< DGFEntityKey< unsigned int>, int > &facemap,
            bool fixedsize,
            int vtxoffset )
    {
      static int cube2simplex[3][3] = {
        {0,1,3},
        {0,2,3},
        {1,2,3}
      };

      int lnofbound;
      int face = ElementFaceUtil::faceSize(dimworld,simplexgrid);
      for (lnofbound=0; ok(); next())
      {
        for (size_t i=0; i<p.size(); i++)
        {
          p[i] -= vtxoffset;
        }
        if (fixedsize)
        {
          if ((dimworld==2 && size()< 2) ||
              (dimworld==3 && simplexgrid && size()!=3 && size()!=4) ||
              (dimworld==3 && !simplexgrid && size()!=4))
            continue;
          std :: vector< unsigned int > bound( face );
          for (int j=0; j<face; j++) {
            bound[j] = p[j];
          }

          DGFEntityKey< unsigned int > key( bound, false );

          facemap[key] = bndid;
          ++lnofbound;

          if (size()>face)
          {
            assert(dimworld==2 || face==3);
            if (dimworld==3) {
              for (int i=0; i<3; i++)
              {
                for (int j=0; j<face; j++)
                {
                  bound[j] = p[cube2simplex[i][j]];
                }

                DGFEntityKey< unsigned int > key( bound, false );
                facemap[key] = bndid;
                ++lnofbound;
              }
            }
            else
            {
              for (int i=2; i<=size(); i++)
              {
                bound[0] = p[i-1];
                bound[1] = p[i%size()];
                DGFEntityKey< unsigned int > key( bound, false );
                facemap[key] = bndid;
                ++lnofbound;
              }
            }
          }
        }
        else {
          if (dimworld==3) {
            DGFEntityKey< unsigned int > key( p, false );
            facemap[key] = bndid;
            ++lnofbound;
          } else {
            std :: vector< unsigned int > k( 2 );
            for (size_t i=0; i<p.size()-1; i++) {
              k[0]=p[i];
              k[1]=p[(i+1)%p.size()];
              DGFEntityKey< unsigned int > key( k, false );
              facemap[key] = bndid;
              ++lnofbound;
            }
          }
        }
      }
      return lnofbound;
    }


    bool BoundarySegBlock :: next ()
    {
      assert(ok());
      int n=0;
      getnextline();
      if (linenumber()==noflines()) {
        goodline=false;
        return goodline;
      }
      p.clear();
      int x;
      if (getnextentry(x)) {
        bndid = x;
        if (bndid<=0)
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      non-positive boundary id (" << bndid << ") read!");
        }
        while (getnextentry(x))
        {
          p.push_back(x);
          n++;
        }
        // goodline=(n==dimworld+1);
        goodline=true;
        return goodline;
      } else return next();
    }


    // DimBlock
    // --------

    const char *DimBlock :: ID = "Dimensions";


    DimBlock :: DimBlock ( std :: istream &in )
      : BasicBlock ( in, ID )
    {
      if (isempty()) {
        DUNE_THROW(DGFException,
                   "no dimension of world specified!");
      } else {
        getnextline();
        line >> _dim;
        if (_dim<1) {
          DUNE_THROW(DGFException,
                     "negative dimension of world specified!");
        }
        else {
          if (noflines()==1)
            _dimworld=_dim;
          else {
            getnextline();
            line >> _dimworld;
            if (_dimworld < _dim) {
              DUNE_THROW(DGFException,
                         "negative dimension of world smaller than dim!");
            }
          }
        }
      }
    }



    // GridParameterBlock
    // ------------------

    const char *GridParameterBlock :: ID = "GridParameter";


    GridParameterBlock :: GridParameterBlock ( std :: istream &in,
                                               const bool readOverlapAndBnd )
      : BasicBlock(in,ID),
        _periodic(),
        _overlap(0), // default value
        _noClosure(false) // default value
    {
      if (! isempty() )
      {
        if( readOverlapAndBnd  )
        {
          // check overlap
          if (findtoken("overlap"))
          {
            int x;
            if( getnextentry(x) ) _overlap = x;
            else
            {
              dwarn << "GridParameterBlock: found keyword `overlap' but no value, defaulting to `" <<  _overlap  <<"' !\n";
            }

            if (_overlap < 0)
            {
              DUNE_THROW(DGFException,"Negative overlap specified!");
            }
          }
          else
          {
            dwarn << "GridParameterBlock: could not find keyword `overlap' in DGF file, defaulting to `"<<_overlap<<"' !\n";
          }

          // check periodic grid
          if (findtoken("periodic"))
          {
            int x;
            while (getnextentry(x))
            {
              _periodic.insert(x);
            }
          }
          else
          {
            dwarn << "GridParameterBlock: could not find keyword `periodic' in DGF file, defaulting to no periodic boundary! \n";
          }
        }
        else
        {
          // check closure
          if (findtoken("closure"))
          {
            std::string clo;
            if(getnextentry(clo))
            {
              makeupcase(clo);
              if(clo == "NONE")
              {
                _noClosure = true;
              }
            }
          }
          else
          {
            dwarn << "GridParameterBlock: could not find keyword `closure' in DGF file, defaulting to `GREEN' !\n";
          }
        }
      }
    }


    // IntervalBlock
    // -------------

    const char *IntervalBlock :: ID = "Interval";


    IntervalBlock :: IntervalBlock ( std :: istream &in )
      : BasicBlock(in,ID),
        p0_(0),
        p1_(0),
        h_(0),
        nofcells_(0),
        good_(false),
        dimw_(0)
    {
      if(isactive()) {
        getnextline();
        double x;
        while (getnextentry(x)) {
          dimw_++;
        }
        if (dimw_==0) {
          DUNE_THROW(DGFException,
                     "Too few coordinates for point p0 in IntervalBlock");
        } else if (dimw_>3) {
          DUNE_THROW(DGFException,
                     "Interval block only implemented for dimension 1,2, and 3");
        }
        p0_.resize(dimw_);
        p1_.resize(dimw_);
        h_.resize(dimw_);
        nofcells_.resize(dimw_);
        reset();
        next();
      }
    }


    int IntervalBlock :: getVtx ( std::vector<std::vector<double> >& vtx )
    {
      size_t countvtx;
      size_t old_size = vtx.size();
      //fill vtx
      dverb << "reading vertices...";
      vtx.resize(vtx.size()+nofvtx());
      for (countvtx=old_size; countvtx < vtx.size(); countvtx++)
        vtx[countvtx].resize(dimw_);
      int m = old_size;
      if(dimw_ == 3) {
        for(int k=0; k < nofcells_[2]+1; k++)        // z-dir
          for(int j=0; j < nofcells_[1]+1; j++)       // y-dir
            for(int i =0; i < nofcells_[0]+1; i++) // x-dir
            {
              vtx[m][0] = p0_[0] + i*h_[0];
              vtx[m][1] = p0_[1] + j*h_[1];
              vtx[m][2] = p0_[2] + k*h_[2];
              m++;
            }
      }
      else if (dimw_==2)
      {
        for(int j=0; j < nofcells_[1]+1; j++)
          for(int i =0; i < nofcells_[0]+1; i++)
          {
            vtx[m][0] = p0_[0] + i*h_[0];
            vtx[m][1] = p0_[1] + j*h_[1];
            m++;
          }
      } else {
        for(int j=0; j < nofcells_[0]+1; j++) {
          vtx[m][0] = p0_[0] + j*h_[0];
          m++;
        }
      }
      dverb << "done" << std::endl;
      return nofvtx();
    }


    int IntervalBlock
    :: getHexa ( std :: vector< std :: vector< unsigned int > > &simplex, int offset )
    {
      int oldsize=simplex.size();
      //fill simplex with Hexaeder
      size_t counthexa;
      int verticesPerCube = -1;
      int m=oldsize;
      if(dimw_ == 3)
        verticesPerCube = 8;
      else if (dimw_ == 2)
        verticesPerCube = 4;
      else if (dimw_ == 1)
        verticesPerCube = 2;
      else
        DUNE_THROW(DGFException,
                   "Invalid dimension world "<< dimw_ << " in IntervalBlock::getHexa!");

      dverb << "generating hexaeder...";
      simplex.resize(oldsize+nofhexa());
      for (counthexa=m; counthexa < simplex.size(); counthexa++)
        simplex[counthexa].resize(verticesPerCube);
      if(dimw_ == 3) {
        for(int k=0; k < nofcells_[2]; k++)
          for(int j=0; j < nofcells_[1]; j++)
            for(int i =0; i < nofcells_[0]; i++)
            {
              simplex[m][0] = offset+getIndex(i,j,k);
              simplex[m][1] = offset+getIndex(i+1,j,k);
              simplex[m][2] = offset+getIndex(i,j+1,k);
              simplex[m][3] = offset+getIndex(i+1,j+1,k);
              simplex[m][4] = offset+getIndex(i,j,k+1);
              simplex[m][5] = offset+getIndex(i+1,j,k+1);
              simplex[m][6] = offset+getIndex(i,j+1,k+1);
              simplex[m][7] = offset+getIndex(i+1,j+1,k+1);
              m++;
            }
      }
      else if (dimw_==2)
      {
        for(int j=0; j < nofcells_[1]; j++)
        {
          for(int i =0; i < nofcells_[0]; i++)
          {
            simplex[m][0] = offset+getIndex(i,j);
            simplex[m][1] = offset+getIndex(i+1,j);
            simplex[m][2] = offset+getIndex(i,j+1);
            simplex[m][3] = offset+getIndex(i+1,j+1);
            m++;
          }
        }
      }
      else
      {
        for(int i =0; i < nofcells_[0]; i++)
        {
          simplex[m][0] = offset+getIndex(i);
          simplex[m][1] = offset+getIndex(i+1);
          m++;
        }
      }
      dverb << "done" << std::endl;
      assert((size_t)m==simplex.size());
      return nofhexa();
    }


    bool IntervalBlock :: next ()
    {
      if (linenumber()==noflines()-1) {
        good_=false;
        return good_;
      }
      //read p0_
      getnextline();
      double x;
      for(int i = 0; i<dimw_; i++)
        if(getnextentry(x))
          p0_[i] = x;
        else  {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "Too few coordinates for point p0");
        }
      //read p1_
      getnextline();
      for(int i = 0; i<dimw_; i++)
        if(getnextentry(x))
          p1_[i] = x;
        else  {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "Too few coordinates for point p1");
        }

      //find real upper and lower edge
      std::vector<double> p0h(dimw_),p1h(dimw_);  //help variables
      for(int i = 0; i<dimw_; i++) {
        p0h[i] = p0_[i] < p1_[i] ? p0_[i] : p1_[i];
        p1h[i] = p0_[i] > p1_[i] ? p0_[i] : p1_[i];
      }
      p0_ = p0h;
      p1_ = p1h;
      //get numbers of cells for every direction
      getnextline();
      int number;
      for(int i = 0; i<dimw_; i++)
        if(getnextentry(number))
          nofcells_[i] = number;
        else {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "Couldn't detect a number of cells for every direction");
        }
      good_ = true;
      for(int i =0; i < dimw_; i++)
        h_[i] = (p1_[i] - p0_[i])/double(nofcells_[i]);
      //Printing Data on screen
      dverb << "p0 = (";
      for(int i = 0; i<dimw_-1; i++)
        dverb << p0_[i]  <<",";
      dverb << p0_[dimw_-1] <<") \n";
      dverb << "p1 = (";
      for(int i = 0; i<dimw_-1; i++)
        dverb << p1_[i]  <<",";
      dverb << p1_[dimw_-1] <<") \n";
      dverb << "n = (";
      for(int i = 0; i<dimw_-1; i++)
        dverb << nofcells_[i]  <<",";
      dverb << nofcells_[dimw_-1] <<") \n";
      dverb << std::endl;
      return good_;
    }

  } // end namespace dgf

} // end namespace Dune
