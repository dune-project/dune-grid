// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSERBLOCKS_HH
#define DUNE_MACROGRIDPARSERBLOCKS_HH

namespace Dune {

  // *************************************************************
  // Read one block with given identifier from disk
  // and allows the line-wise extraction from this block
  // *************************************************************
  namespace {
    void makeupcase(std::string &s) {
      for (size_t i=0; i<s.size(); i++)
        s[i]=toupper(s[i]);
    }
    class BasicBlock {
      int pos;               // line number
      bool active;           // block was found
      bool empty;            // block was found but was empty
      std::string identifier; // identifier of this block
      int linecount;         // total number of lines in the block
      std::stringstream block; // the block itself
      std::string oneline;   // the active line in the block
      // get the block (if it exists)
      void getblock(std::istream &in) {
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
      // count the number of lines in the block
      int countlines() {
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
    protected:
      std::stringstream line; // the active line as string buffer
                              // for use in the derived classes
      // go back to beginning of block
      void reset() {
        pos=-1;
        block.clear();
        block.seekg(0);
        linecount=countlines();
        pos=-1;
        block.clear();
        block.seekg(0);
      }
      // get next line and store in string stream
      void getnextline() {
        line.clear();
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
      // get next entry in line
      template <class ENTRY>
      bool getnextentry(ENTRY &entry) {
        line >> entry;
        return line;
      }
      bool gettokenparam(std::string token,std::string& entry) {
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
      bool findtoken(std::string token) {
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
    public:
      // search for block in file and store in buffer
      BasicBlock(std::istream& in, const char* id) :
        pos(-1),
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
      // some information on this block
      bool isactive() {
        return active;
      }
      bool isempty() {
        return empty;
      }
      int& noflines() {
        return linecount;
      }
      int linenumber() {
        return pos;
      }
      // for error messages
      friend std::ostream& operator<<(std::ostream& os, const BasicBlock &b) {
        return os
               << "block " << b.identifier
               << " on line " << b.pos << std::endl;
      }
    };
    // *************************************************************
    // derived classes for each block in grid file
    // *************************************************************
    class VertexBlock : public BasicBlock {
      int dimworld;      // the dimesnsion of the verticies (is given from user)
      bool goodline;     // active line describes a vertex
      std::vector<double> p; // active vertex
      int vtxoffset;
    public:
      static const char* ID;
      // initialize vertex block and get first vertex
      VertexBlock(std::istream& in,int &pdimworld) :
        BasicBlock(in,ID),
        dimworld(pdimworld),
        goodline(true),
        p(0),
        vtxoffset(0)
      {
        /*
           if (dimworld==-1) {
           assert(ok());
           getnextline();
           dimworld=0;
           double x;
           while (getnextentry(x))
            dimworld++;
           pdimworld=dimworld;
           reset();
           }
           p.resize(dimworld);
         */
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
        reset();
        next();
        pdimworld=dimworld;
      }
      ~VertexBlock() {}
      int offset() {
        return vtxoffset;
      }
      int get(std::vector<std::vector<double> >& vtx) {
        size_t nofvtx;
        size_t old_size = vtx.size();
        // vtx.resize(old_size+nofvertex());
        for (nofvtx=old_size; ok(); next(),nofvtx++) {
          vtx.push_back(p);
          /*
             vtx[nofvtx].resize(dimworld);
             for (int j=0;j<dimworld;j++) {
             vtx[nofvtx][j] = p[j];
             }
           */
        }
        /*
           if (nofvtx!=vtx.size()) {
           DUNE_THROW(DGFException, "Wrong number of verticies read!");
           }
         */
        return nofvtx;
      }
      // some information
      bool ok() {
        return goodline;
      }
      /*
         int nofvertex() {
         return noflines();
         }
       */
    private:
      // get next vertex
      bool next() {
        assert(ok());
        int n=0;
        getnextline();
        if (linenumber()==noflines()) {
          goodline=false;
          return goodline;
        }
        int olddimworld=dimworld;
        double x;
        while (getnextentry(x)) {
          if (n<dimworld)
            p[n]=x;
          else {
            p.push_back(x);
            dimworld++;
          }
          n++;
        }
        if (olddimworld>0 && dimworld!=olddimworld) {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      wrong number of coordinates: "
                                 << n << " read but expected " << dimworld);
        }
        if (n!=dimworld || dimworld<1) {
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
      // get coordinates of active vertex
      double operator[](int i) {
        assert(ok());
        assert(linenumber()>=0);
        assert(0<=i && i<dimworld);
        return p[i];
      }
    };
    const char *VertexBlock::ID = "Vertex";
    // *************************************************************
    class SimplexGenerationBlock : public BasicBlock {
      double area_;
      double angle_;
      bool display_;
      std::string path_;
      bool haspath_;
      std::string filename_;
      std::string filetype_;
      std::string parameter_;
      bool hasfile_;
      int dimension_;
    public:
      const static char* ID;
      SimplexGenerationBlock(std::istream& in) :
        BasicBlock(in,ID),
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
      double maxArea() {
        return area_;
      }
      double minAngle() {
        return angle_;
      }
      bool display() {
        return display_;
      }
      bool haspath() {
        return haspath_;
      }
      std::string path() {
        return path_;
      }
      bool hasfile() {
        return hasfile_;
      }
      std::string filename() {
        return filename_;
      }
      std::string filetype() {
        return filetype_;
      }
      int dimension() {
        return dimension_;
      }
      std::string parameter() {
        return parameter_;
      }
    };
    const char* SimplexGenerationBlock::ID = "Simplexgenerator";
    // *************************************************************
    class SimplexBlock : public BasicBlock {
      int nofvtx;
      int dimworld;
      bool goodline;  // active line describes a vertex
      std::vector<int> p; // active vertex
    public:
      const static char* ID;
      SimplexBlock(std::istream& in,int pnofvtx, int adimworld) :
        BasicBlock(in,ID),
        nofvtx(pnofvtx),
        dimworld(adimworld),
        goodline(true),
        p(adimworld+1)
      {
        assert((dimworld+1)>0);
        next();
      }
      ~SimplexBlock() {}
      int get(std::vector<std::vector<int> >& simplex,int vtxoffset)
      {
        int nofsimpl;
        simplex.resize(nofsimplex());
        int offset = 0;
        bool first = true;
        for (nofsimpl=0; ok(); next(), nofsimpl++)
        {
          assert(nofsimpl<nofsimplex());
          simplex[nofsimpl].resize(dimworld+1);
          if(first)
          {
            offset = p[0];
            first = false;
          }
          for (int j=0; j<dimworld+1; j++)
          {
            simplex[nofsimpl][j] = p[j];
            offset = std::min(p[j],offset);
          }
        }

        if (nofsimpl!=nofsimplex()) {
          DUNE_THROW(DGFException,
                     "Wrong number of simplex elements read!" << std::endl
                                                              << "        Read " << nofsimpl
                                                              << " elements, expected " << nofsimplex());
        }
        // make numbering starting from zero
        // not matter whether offset is positive or negative
        offset = vtxoffset;
        if(offset != 0)
        {
          for (int i=0; i<nofsimpl; ++i)
          {
            for (size_t j=0; j<simplex[i].size(); ++j)
            {
              simplex[i][j] -= offset;
            }
          }
        }
        return nofsimpl;
      }
      // cubes -> simplex
      static int cube2simplex(std::vector<std::vector<double> >& vtx,
                              std::vector<std::vector<int> >& elements) {
        static int offset3[6][4][3] = {{{0,0,0},{1,1,1},{1,0,0},{1,1,0}},
                                       {{0,0,0},{1,1,1},{1,0,1},{1,0,0}},
                                       {{0,0,0},{1,1,1},{0,0,1},{1,0,1}},
                                       {{0,0,0},{1,1,1},{1,1,0},{0,1,0}},
                                       {{0,0,0},{1,1,1},{0,1,0},{0,1,1}},
                                       {{0,0,0},{1,1,1},{0,1,1},{0,0,1}} };
        static int offset2[2][3][2] = {{{0,0},{1,0},{0,1}},
                                       {{1,1},{0,1},{1,0}}};
        if (vtx.size()==0)
          DUNE_THROW(DGFException, "Converting Cune- to Simplexgrid with no verticies given");
        int dimworld = vtx[0].size();
        dverb << "generating simplices...";
        dverb.flush();
        std::vector<std::vector<int> > cubes = elements;
        if(dimworld == 3) {
          elements.resize(6*cubes.size());
          for(size_t countsimpl=0; countsimpl < elements.size(); countsimpl++)
            elements[countsimpl].resize(4);
          for (size_t c=0; c<cubes.size(); c++)
          {
            for(int tetra=0; tetra < 6 ; tetra++)
            {
              for (int v=0; v<4; v++) {
                elements[c*6+tetra][v]=
                  cubes[c][offset3[tetra][v][0]+
                           offset3[tetra][v][1]*2+
                           offset3[tetra][v][2]*4];
              }
            }
          }
        }
        else {
          elements.resize(2*cubes.size() );
          for(size_t countsimpl=0; countsimpl < elements.size(); countsimpl++)
            elements[countsimpl].resize(3);
          for (size_t c=0; c<cubes.size(); c++)
          {
            int diag = 0;
            double mind = 0;
            for (int d=0; d<2; d++)
            {
              double diaglen =
                pow(vtx[cubes[c][d]][0]-vtx[cubes[c][2+((d+1)%2)]][0],2) +
                pow(vtx[cubes[c][d]][1]-vtx[cubes[c][2+((d+1)%2)]][1],2);
              if (diaglen<mind)
              {
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
            }
          }
        }
        return elements.size();
      }
      // some information
      bool ok() {
        return goodline;
      }
      int nofsimplex() {
        return noflines();
      }
    private:
      // get next simplex
      bool next() {
        assert(ok());
        int n=0;
        getnextline();
        if (linenumber()==noflines()) {
          goodline=false;
          return goodline;
        }
        int x;
        while (getnextentry(x)) {
          if (n<(dimworld+1)) {
            p[n]=x;
            if (n< 0 && x>= nofvtx) {
              DUNE_THROW(DGFException,
                         "ERROR in " << *this
                                     << "      wrong index of vertices: "
                                     << x << " read but expected value between 0 and"
                                     << nofvtx-1);
            }
          }
          n++;
        }
        // tests if the written block is ok in its size
        int dimnew = dimworld +1 ;
        goodline=(n==(dimworld+1));
        if (!goodline) {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      wrong number of vertices: "
                                 << n << " read but expected " << dimnew);
        }
        return goodline;
      }
      // get coordinates of active simplex
      int operator[](int i) {
        assert(ok());
        assert(linenumber()>=0);
        assert(0<=i && i<(dimworld+1));
        return p[i];
      }
    };
    const char* SimplexBlock::ID = "Simplex";
    /// *************************************************************
    class CubeBlock : public BasicBlock {
      int nofvtx;
      int dimworld;
      bool goodline;    // active line describes a vertex
      std::vector<int> p; // active vertex
      std::vector<int> map; // active vertex
    public:
      static const char* ID;
      CubeBlock(std::istream& in,int pnofvtx, int adimworld) :
        BasicBlock(in,ID),
        nofvtx(pnofvtx),
        dimworld(adimworld),
        goodline(true),
        p(1<<adimworld),
        map(1<<adimworld)
      {
        assert((dimworld+1)>0);
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
        reset();
        next();
      }
      ~CubeBlock() {}
      int get(std::vector<std::vector<int> >& simplex,int vtxoffset) {
        int nofsimpl;
        // simplex.resize(nofsimplex());
        for (nofsimpl=0; ok(); next(), nofsimpl++) {
          simplex.push_back(p);
          for (size_t j=0; j<p.size(); j++) {
            simplex[nofsimpl][map[j]] = p[j];
          }
        }
        int offset = vtxoffset;
        if(offset != 0)
        {
          for (int i=0; i<nofsimpl; ++i)
          {
            for (size_t j=0; j<simplex[i].size(); ++j)
            {
              simplex[i][j] -= offset;
            }
          }
        }
        /*
           if (nofsimpl!=nofsimplex()){
           DUNE_THROW(DGFException,
                     "Error occured while reading element information!"
                     << std::endl
                     << "        expected " << nofsimplex()
                     << " element information "
                     << "        but only read " << nofsimpl << " elements.");
           }
         */
        return nofsimpl;
      }
      // some information
      bool ok() {
        return goodline;
      }
      int nofsimplex() {
        return noflines();
      }
    private:
      // get next simplex
      bool next() {
        assert(ok());
        int n=0;
        getnextline();
        if (linenumber()==noflines()) {
          goodline=false;
          return goodline;
        }
        int x;
        while (getnextentry(x))
        {
          if (n<(int)p.size())
          {
            p[n]=x;
            if (x< 0 && x>= nofvtx)
            {
              DUNE_THROW(DGFException,
                         "ERROR in " << *this
                                     << "      wrong index of vertices: "
                                     << x << " read but expected value between 0 and "
                                     << nofvtx-1);
            }
          }
          n++;
        }
        // tests if the written block is ok in its size
        if (n!=(int)p.size()) {
          return next();
          /*
             DUNE_THROW(DGFException,
                     "ERROR in " << *this
                     << "      wrong number of vertices: "
                     << "      read " << n << " but expected " << p.size());
           */
        }
        goodline=(n==(int)p.size());
        return goodline;
      }
      // get coordinates of active simplex
      int operator[](int i) {
        assert(ok());
        assert(linenumber()>=0);
        assert(0<=i && i< (int)p.size());
        return p[i];
      }
    };
    const char* CubeBlock::ID = "Cube";
    // *************************************************************
    // the block BoundaryDomBlock looks for a domain which is characterized by two points in R^dimworld
    class BoundaryDomBlock : public BasicBlock {
      int dimworld;    // the dimesnsion of the verticies (is given  from user)
      bool goodline;   // active line describes a vertex
      std::vector<double> p1,p2; // active vertex
      int bndid;
      bool withdefault;
      int defaultvalue;
    public:
      static const char* ID;
      // initialize vertex block and get first vertex
      BoundaryDomBlock(std::istream& in,int cdimworld ) :
        BasicBlock(in,ID),
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
            if (getnextentry(x)) {
              defaultvalue=x;
              withdefault = true;
            }
        }
        reset();
        next();
      }
      ~BoundaryDomBlock() {}
      bool next() {
        assert(ok());
        getnextline();
        if (linenumber()==noflines()) {
          goodline=false;
          return goodline;
        }
        int id;
        if (getnextentry(id)) {
          bndid = id;
          double x;
          int n=0;
          while (getnextentry(x)) {
            if (0<=n && n<dimworld)
              p1.at(n)=x;
            else if (dimworld<=n && n<2*dimworld) {
              p2.at(n-dimworld)=x;
              if (p2.at(n-dimworld)<p1.at(n-dimworld)) {
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
          if (!goodline) {
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
      bool inside(const std::vector<double>& v) const {
        assert(v.size()==(size_t)dimworld);
        for (int i=0; i<dimworld; i++)
          if (v[i]<p1[i] || v[i]>p2[i])
            return false;
        return true;
      }
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
    const char *BoundaryDomBlock::ID = "boundarydomain";
    // *************************************************************
    // the BoundarySegBlock looks for given boundary values unless they aren't given they got the value zero
    class BoundarySegBlock : public BasicBlock {
      int dimworld;      // the dimesnsion of the verticies (is given  from user)
      bool goodline;     // active line describes a vertex
      std::vector<int> p; // active vertex
      int bndid;
      bool simplexgrid;
    public:
      static const char* ID;
      // initialize vertex block and get first vertex
      BoundarySegBlock(std::istream& in,int pnofvtx,
                       int pdimworld,bool psimplexgrid ) :
        BasicBlock(in,ID),
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
      ~BoundarySegBlock() {}
      // some information
      int get(std::map<EntityKey<int>,int>& facemap,bool fixedsize,int vtxoffset) {
        static int cube2simplex[3][3] = {
          {0,1,3},
          {0,2,3},
          {1,2,3} };
        int lnofbound;
        int face=ElementFaceUtil::faceSize(dimworld,simplexgrid);
        for (lnofbound=0; ok(); next()) {
          for (size_t i=0; i<p.size(); i++) {
            p[i] -= vtxoffset;
          }
          if (fixedsize) {
            if ((dimworld==2 && size()<=2) ||
                (dimworld==3 && simplexgrid && size()!=3 && size()!=4) ||
                (dimworld==3 && !simplexgrid && size()!=4))
              continue;
            std::vector<int> bound(face);
            for (int j=0; j<face; j++) {
              bound[j] = p[j];
            }
            EntityKey<int> key(bound,false);
            facemap[key] = bndid;
            ++lnofbound;
            if (size()>face) {
              assert(dimworld==2 || face==3);
              if (dimworld==3) {
                for (int i=0; i<3; i++) {
                  for (int j=0; j<face; j++) {
                    bound[j] = p[cube2simplex[i][j]];
                  }
                  EntityKey<int> key(bound,false);
                  facemap[key] = bndid;
                  ++lnofbound;
                }
              } else {
                for (int i=2; i<=size(); i++) {
                  bound[0] = p[i-1];
                  bound[1] = p[i%size()];
                  EntityKey<int> key(bound,false);
                  facemap[key] = bndid;
                  ++lnofbound;
                }
              }
            }
          }
          else {
            if (dimworld==3) {
              EntityKey<int> key(p,false);
              facemap[key] = bndid;
              ++lnofbound;
            } else {
              std::vector<int> k(2);
              for (size_t i=0; i<p.size()-1; i++) {
                k[0]=p[i];
                k[1]=p[(i+1)%p.size()];
                EntityKey<int> key(k,false);
                facemap[key] = bndid;
                ++lnofbound;
              }
            }
          }
        }
        return lnofbound;
      }
      bool ok() {
        return goodline;
      }
      int nofbound() {
        return noflines();
      }
    private:
      bool next() {
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
          if (bndid<=0) {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      non-positive boundary id read!");
          }
          while (getnextentry(x)) {
            p.push_back(x);
            n++;
          }
          // goodline=(n==dimworld+1);
          goodline=true;
          return goodline;
        } else return next();
      }
      // get coordinates of active vertex
      int operator[](int i) {
        assert(ok());
        assert(linenumber()>=0);
        assert(0<=i && i<dimworld+1);
        return p[i];
      }
      int size() {
        return p.size();
      }
    };
    const char *BoundarySegBlock::ID = "boundarysegments";
    // *************************************************************
    class DimBlock : public BasicBlock {
      int _dimworld; // dimension of world
      int _dim;      // dimension of grid
    public:
      const static char* ID;
      // initialize block and get dimension of world
      DimBlock(std::istream& in) :
        BasicBlock(in,ID)
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
      // get dimension of world found in block
      int dim() {
        return _dim;
      }
      int dimworld() {
        return _dimworld;
      }
      // some information
      bool ok() {
        return true;
      }
    };
    const char* DimBlock::ID = "Dimensions";
    // *************************************************************
    class IntervalBlock : public BasicBlock {
      std::vector<double> p0_,p1_; //lower and upper boundary points
      std::vector<double> h_;      // width of the cells in every direction
      std::vector<int> nofcells_;  // number of cells in every direction
      bool good_;                  //data read correctly
      int dimw_;                   //dimension of world
    public:
      const static char* ID;
      IntervalBlock(std::istream& in) :
        BasicBlock(in,ID),
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
          }
          p0_.resize(dimw_);
          p1_.resize(dimw_);
          h_.resize(dimw_);
          nofcells_.resize(dimw_);
          reset();
          next();
        }
      }
      void get(std::vector<std::vector<double> >& vtx,int& nofvtx,
               std::vector<std::vector<int> >& simplex,int& nofsimpl) {
        do {
          int oldvtx = nofvtx;
          nofvtx  +=getVtx(vtx);
          nofsimpl+=getHexa(simplex,oldvtx);
        } while (next());
      }
      void get(std::vector<std::vector<double> >& vtx,int& nofvtx) {
        do {
          // int oldvtx = nofvtx;
          nofvtx  +=getVtx(vtx);
        } while (next());
      }
      int getVtx(std::vector<std::vector<double> >& vtx) {
        size_t countvtx;
        size_t old_size = vtx.size();
        //fill vtx
        dverb << "reading vertices...";
        vtx.resize(vtx.size()+nofvtx());
        for (countvtx=old_size; countvtx < vtx.size(); countvtx++)
          vtx[countvtx].resize(dimw_);
        int m = old_size;
        if(dimw_ == 3) {
          for(int i =0; i < nofcells_[0]+1; i++)
            for(int j=0; j < nofcells_[1]+1; j++)
              for(int k=0; k < nofcells_[2]+1; k++) {
                vtx[m][0] = p0_[0] + i*h_[0];
                vtx[m][1] = p0_[1] + j*h_[1];
                vtx[m][2] = p0_[2] + k*h_[2];
                m++;
              }
        }
        else {
          for(int i =0; i < nofcells_[0]+1; i++)
            for(int j=0; j < nofcells_[1]+1; j++) {
              vtx[m][0] = p0_[0] + i*h_[0];
              vtx[m][1] = p0_[1] + j*h_[1];
              m++;
            }
        }
        dverb << "done" << std::endl;
        return nofvtx();
      }
      int getHexa(std::vector<std::vector<int> >& simplex,
                  int offset=0) {
        int oldsize=simplex.size();
        //fill simplex with Hexaeder
        size_t counthexa;
        int verticesPerCube;
        int m=oldsize;
        if(dimw_ == 3)
          verticesPerCube = 8;
        else
          verticesPerCube = 4;
        dverb << "generating hexaeder...";
        simplex.resize(oldsize+nofhexa());
        for (counthexa=m; counthexa < simplex.size(); counthexa++)
          simplex[counthexa].resize(verticesPerCube);
        if(dimw_ == 3) {
          for(int i =0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++)
              for(int k=0; k < nofcells_[2]; k++) {
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
        else {
          for(int i =0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++) {
              simplex[m][0] = offset+getIndex(i,j);
              simplex[m][1] = offset+getIndex(i+1,j);
              simplex[m][2] = offset+getIndex(i,j+1);
              simplex[m][3] = offset+getIndex(i+1,j+1);
              m++;
            }
        }
        dverb << "done" << std::endl;
        assert((size_t)m==simplex.size());
        return nofhexa();
      }

      int nofvtx() {
        if(dimw_ == 3)
          return (nofcells_[0]+1)*(nofcells_[1]+1)*(nofcells_[2]+1);
        else
          return (nofcells_[0]+1)*(nofcells_[1]+1);
      }

      int nofhexa() {
        if(dimw_ == 3)
          return (nofcells_[0])*(nofcells_[1])*(nofcells_[2]);
        else
          return (nofcells_[0])*(nofcells_[1]);

      }

      int segments(int i) {
        return nofcells_[i];
      }

      double length(int i) {
        return p1_[i]-p0_[i];
      }

      double start(int i) {
        return p0_[i];
      }

      double end(int i) {
        return p1_[i];
      }
      /*
         bool ok() {
         return good_;
         }
       */
      int dimw() {
        return dimw_;
      }

      int getIndex(int i,int j, int k = 0) {
        if(dimw_ == 3)
          return i*(nofcells_[1]+1)*(nofcells_[2]+1) + j*(nofcells_[2]+1) + k;
        else
          return i*(nofcells_[1]+1) + j;
      }
    private:
      bool next() {
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
        std::vector<double> p0h(dimw_),p1h(dimw_); //help variables
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
    };
    const char* IntervalBlock::ID = "Interval";
  }
}
#endif
