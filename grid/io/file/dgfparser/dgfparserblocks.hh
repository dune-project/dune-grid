// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSERBLOCKS_HH
#define DUNE_MACROGRIDPARSERBLOCKS_HH

#include "entitykey.hh"

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
            std::cerr << "Error: block must end with a #-line" << std::endl;
            abort();
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
          std::cerr << "ERROR: file not found in BasicBlock::BasicBlock" << std::endl;
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
      static const char* ID;
      int dimworld;      // the dimesnsion of the verticies (is given from user)
      bool goodline;     // active line describes a vertex
      std::vector<double> p; // active vertex
    public:
      // initialize vertex block and get first vertex
      VertexBlock(std::istream& in,int &pdimworld) :
        BasicBlock(in,ID),
        dimworld(pdimworld),
        goodline(true),
        p(0)
      {
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
        next();
      }
      ~VertexBlock() {}
      int get(std::vector<std::vector<double> >& vtx) {
        size_t nofvtx;
        size_t old_size = vtx.size();
        vtx.resize(old_size+nofvertex());
        for (nofvtx=old_size; ok(); next(),nofvtx++) {
          vtx[nofvtx].resize(dimworld);
          for (int j=0; j<dimworld; j++) {
            vtx[nofvtx][j] = p[j];
          }
        }
        if (nofvtx!=vtx.size()) {
          std::cerr << "ERROR: Wrong number of verticies read!" << std::endl;
          return 0;
        }
        return nofvertex();
      }
      // some information
      bool ok() {
        return goodline;
      }
      int nofvertex() {
        return noflines();
      }
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
        double x;
        while (getnextentry(x)) {
          if (n<dimworld)
            p[n]=x;
          n++;
        }
        goodline=(n==dimworld);
        if (!goodline) {
          std::cerr << "ERROR in " << *this
                    << "      wrong number of coordinates: "
                    << n << " read but expected " << dimworld << std::endl;
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
      const static char* ID;
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
      const static char* ID;
      int nofvtx;
      int dimworld;
      bool goodline;  // active line describes a vertex
      std::vector<int> p; // active vertex
    public:
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
      int get(std::vector<std::vector<int> >& simplex) {
        int nofsimpl;
        simplex.resize(nofsimplex());
        for (nofsimpl=0; ok(); next(), nofsimpl++) {
          assert(nofsimpl<nofsimplex());
          simplex[nofsimpl].resize(dimworld+1);
          for (int j=0; j<dimworld+1; j++) {
            simplex[nofsimpl][j] = p[j];
          }
        }
        if (nofsimpl!=nofsimplex()) {
          std::cerr << "ERROR:  Wrong number of simplex elements read!" << std::endl;
          std::cerr << "        Read " << nofsimpl << " elements, expected " << nofsimplex() << std::endl;
          return 0;
        }
        return nofsimpl;
      }
      // cubes -> simplex
      int cube2simplex(std::vector<std::vector<double> >& vtx,
                       std::vector<std::vector<int> >& elements) {
        static int offset3[6][4][3] = {{{0,0,0},{1,1,1},{1,0,0},{1,1,0}},
                                       {{0,0,0},{1,1,1},{1,0,1},{1,0,0}},
                                       {{0,0,0},{1,1,1},{0,0,1},{1,0,1}},
                                       {{0,0,0},{1,1,1},{1,1,0},{0,1,0}},
                                       {{0,0,0},{1,1,1},{0,1,0},{0,1,1}},
                                       {{0,0,0},{1,1,1},{0,1,1},{0,0,1}} };
        static int offset2[2][3][2] = {{{0,0},{1,0},{0,1}},
                                       {{1,1},{0,1},{1,0}}};
        std::cout << "generating simplices...";
        std::cout.flush();
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
            if (n< 0 && x>= nofvtx) { std::cerr << "ERROR in " << *this
                                                << "      wrong index of vertices: "
                                                << x << " read but expected value between 0 and" << nofvtx-1 << std::endl;}
          }
          n++;
        }
        // tests if the written block is ok in its size
        int dimnew = dimworld +1 ;
        goodline=(n==(dimworld+1));
        if (!goodline) {
          std::cerr << "ERROR in " << *this
                    << "      wrong number of vertices: "
                    << n << " read but expected " << dimnew << std::endl;
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
      const static char* ID;
      int nofvtx;
      int dimworld;
      bool goodline;  // active line describes a vertex
      std::vector<int> p; // active vertex
    public:
      CubeBlock(std::istream& in,int pnofvtx, int adimworld) :
        BasicBlock(in,ID),
        nofvtx(pnofvtx),
        dimworld(adimworld),
        goodline(true),
        p(1<<adimworld)
      {
        assert((dimworld+1)>0);
        next();
      }
      ~CubeBlock() {}
      int get(std::vector<std::vector<int> >& simplex) {
        int nofsimpl;
        simplex.resize(nofsimplex());
        for (nofsimpl=0; ok(); next(), nofsimpl++) {
          simplex[nofsimpl].resize(p.size());
          for (size_t j=0; j<p.size(); j++) {
            simplex[nofsimpl][j] = p[j];
          }
        }
        if (nofsimpl!=nofsimplex()) {
          std::cerr << "ERROR:  Error occured while reading element information!" << std::endl;
          std::cerr << "        expected " << nofsimplex() << " element information " << std::endl;
          std::cerr << "        but only read " << nofsimpl << " elements." << std::endl;
          return 0;
        }
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
            { std::cerr << "ERROR in " << *this
                        << "      wrong index of vertices: "
                        << x << " read but expected value between 0 and" << nofvtx-1 << std::endl;}
          }
          n++;
        }
        // tests if the written block is ok in its size
        goodline=(n==(int)p.size());
        if (!goodline) {
          std::cerr << "ERROR in " << *this
                    << "      wrong number of vertices: "
                    << "      read " << n << " but expected " << p.size() << std::endl;
        }
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
      const static char* ID;
      int dimworld;    // the dimesnsion of the verticies (is given  from user)
      bool goodline;   // active line describes a vertex
      std::vector<double> p1,p2; // active vertex
      int bndid;
      bool withdefault;
      int defaultvalue;
    public:
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
                std::cerr << "ERROR in " << *this
                          << "      second coordinate smaller than first coordinate: "
                          << p2.at(n-dimworld)
                          << " read but expected value larger or equal to "
                          << p1.at(n-dimworld) << std::endl;
                goodline=false;
                return goodline;
              }
            }
            n++;
          }
          goodline=(n==dimworld*2);
          if (!goodline) {
            std::cerr << "ERROR in " << *this
                      << "      wrong number of coordinates: "
                      << n << " read but expected " << dimworld << std::endl;
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
      const static char* ID;
      int dimworld;      // the dimesnsion of the verticies (is given  from user)
      bool goodline;     // active line describes a vertex
      std::vector<int> p; // active vertex
      int bndid;
    public:
      // initialize vertex block and get first vertex
      BoundarySegBlock(std::istream& in,int pnofvtx, int cdimworld ) :
        BasicBlock(in,ID),
        dimworld(cdimworld),
        goodline(true),
        p(cdimworld+1)
      {
        assert(dimworld>0);
        next();
      }
      ~BoundarySegBlock() {}
      // some information
      bool ok() {
        return goodline;
      }
      int nofbound() {
        return noflines();
      }
      // private:
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
          if (0<=n && n<dimworld+1)
            p[n]=x;
          n++;
        }
        goodline=(n==dimworld+1);
        if (!goodline) {
          std::cerr << "ERROR in " << *this
                    << "      wrong number of coordinates: "
                    << n << " read but expected " << dimworld << std::endl;
        }
        return goodline;
      }



      // get coordinates of active vertex
      int operator[](int i) {
        assert(ok());
        assert(linenumber()>=0);
        assert(0<=i && i<dimworld+1);
        return p[i];
      }
    };
    const char *BoundarySegBlock::ID = "boundarysegments";
    // *************************************************************
    class DimBlock : public BasicBlock {
      const static char* ID;
      int _dimworld; // dimension of world
      int _dim;      // dimension of grid
      bool good;     // dimension of world > 0, dim <= dimworld
    public:
      // initialize block and get dimension of world
      DimBlock(std::istream& in) :
        BasicBlock(in,ID)
      {
        if (isempty()) {
          // std::cerr << "ERROR: no dimension of world specified!" << std::endl;
          good=false;
        } else {
          getnextline();
          line >> _dim;
          if (_dim<1) {
            std::cerr << "ERROR: negative dimension of world specified!" << std::endl;
            good=false;
          }
          else {
            good=true;
            if (noflines()==1)
              _dimworld=_dim;
            else {
              getnextline();
              line >> _dimworld;
              if (_dimworld < _dim) {
                std::cerr << "ERROR: dimension of world smaller than dim!" << std::endl;
                good=false;
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
        return good;
      }
    };
    const char* DimBlock::ID = "Dimensions";
    // *************************************************************
    class IntervalBlock : public BasicBlock {
      const static char* ID;
      std::vector<double> p0_,p1_; //lower and upper boundary points
      std::vector<double> h_;    // width of the cells in every direction
      int nofcells_[3];                // number of cells in every direction
      bool good_;                  //data read correctly
      int dimw_;                   //dimension of world
    public:
      IntervalBlock(std::istream& in) :
        BasicBlock(in,ID),
        p0_(0),
        p1_(0),
        h_(0),
        good_(false),
        dimw_(0)
      {
        if(isactive()) {
          //read p0_
          getnextline();
          double x;
          while (getnextentry(x)) {
            p0_.push_back(x);
            dimw_++;
          }
          if (dimw_==0) {
            std::cerr << "ERROR: Too few coordinates for lower point p0" << std::endl;
            abort();
          }
          p1_.resize(dimw_);
          h_.resize(dimw_);
          //read p1_
          getnextline();
          for(int i = 0; i<dimw_; i++)
            if(getnextentry(x))
              p1_[i] = x;
            else  {
              std::cerr << "ERROR :Too few coordinates for upper point p1" << std::endl;
              abort();
            }
          //assert((p0_[0] < p1_[0]) && (p0_[1] < p1_[1]) && (p0_[2] < p1_[2]));

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
              std::cerr << "ERROR: Couldn't detect a number of cells for every direction"
                        << std::endl;
              abort();
            }
          good_ = true;
          for(int i =0; i < dimw_; i++)
            h_[i] = (p1_[i] - p0_[i])/nofcells_[i];
          //Printing Data on screen
          std::cout << "p0 = (";
          for(int i = 0; i<dimw_-1; i++)
            std::cout << p0_[i]  <<",";
          std::cout << p0_[dimw_-1] <<") \n";
          std::cout << "p1 = (";
          for(int i = 0; i<dimw_-1; i++)
            std::cout << p1_[i]  <<",";
          std::cout << p1_[dimw_-1] <<") \n";
          std::cout << "n = (";
          for(int i = 0; i<dimw_-1; i++)
            std::cout << nofcells_[i]  <<",";
          std::cout << nofcells_[dimw_-1] <<") \n";
          std::cout << std::endl;
        }
      }
      int getVtx(std::vector<std::vector<double> >& vtx) {
        size_t countvtx;
        size_t old_size = vtx.size();
        //fill vtx
        std::cout << "reading vertices...";
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
        std::cout << "done" << std::endl;
        return nofvtx();
      }

      int getSimplex(std::vector<std::vector<int> >& simplex) {
        //fill simplex
        size_t countsimpl;
        int m=0;
        //static int offset[6][4][3] = {{{0,0,0},{1,1,1},{1,0,0},{1,1,0}},
        //				  {{0,0,0},{1,1,1},{1,0,0},{1,0,1}},
        //			  {{0,0,0},{1,1,1},{0,0,1},{1,0,1}},
        //			  {{0,0,0},{1,1,1},{0,1,0},{1,1,0}},
        //			  {{0,0,0},{1,1,1},{0,1,0},{0,1,1}},
        //			  {{0,0,0},{1,1,1},{0,0,1},{0,1,1}}};
        static int offset[6][4][3] = {{{0,0,0},{1,1,1},{1,0,0},{1,1,0}},
                                      {{0,0,0},{1,1,1},{1,0,1},{1,0,0}},
                                      {{0,0,0},{1,1,1},{0,0,1},{1,0,1}},
                                      {{0,0,0},{1,1,1},{1,1,0},{0,1,0}},
                                      {{0,0,0},{1,1,1},{0,1,0},{0,1,1}},
                                      {{0,0,0},{1,1,1},{0,1,1},{0,0,1}} };
        std::cout << "generating simplices...";
        std::cout.flush();
        if(dimw_ == 3) {
          simplex.resize(6*nofhexa() );
          for(countsimpl=0; countsimpl < simplex.size(); countsimpl++)
            simplex[countsimpl].resize(4);
          for(int i =0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++)
              for(int k=0; k < nofcells_[2]; k++) {
                for(int tetra=0; tetra < 6 ; tetra++)
                  for(int vert = 0 ; vert < 4 ; vert++)
                    simplex[m+tetra][vert] =
                      getIndex(i+offset[tetra][vert][0],
                               j+offset[tetra][vert][1],
                               k+offset[tetra][vert][2]);
                m+=6;
              }
        }
        else {
          simplex.resize(2*nofhexa() );
          for(countsimpl=0; countsimpl < simplex.size(); countsimpl++)
            simplex[countsimpl].resize(3);
          for(int i =0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++) {
              simplex[m][0] = getIndex(i+1,j);
              simplex[m][1] = getIndex(i+1,j+1);
              simplex[m][2] = getIndex(i,j);
              simplex[m+1][0] = getIndex(i,j+1);
              simplex[m+1][1] = getIndex(i,j);
              simplex[m+1][2] = getIndex(i+1,j+1);
              m+=2;
            }
        }

        std::cout << "done" << std::endl;
        std::cout.flush();
        return simplex.size();
      }
      int getHexa(std::vector<std::vector<int> >& simplex) {
        //fill simplex with Hexaeder
        size_t counthexa;
        int verticesPerCube;
        int m=0;
        if(dimw_ == 3)
          verticesPerCube = 8;
        else
          verticesPerCube = 4;
        std::cout << "generating hexaeder...";
        simplex.resize(nofhexa());
        for (counthexa=0; counthexa < simplex.size(); counthexa++)
          simplex[counthexa].resize(verticesPerCube);
        if(dimw_ == 3) {
          for(int i =0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++)
              for(int k=0; k < nofcells_[2]; k++) {
                simplex[m][0] = getIndex(i,j,k);
                simplex[m][1] = getIndex(i+1,j,k);
                simplex[m][2] = getIndex(i,j+1,k);
                simplex[m][3] = getIndex(i+1,j+1,k);
                simplex[m][4] = getIndex(i,j,k+1);
                simplex[m][5] = getIndex(i+1,j,k+1);
                simplex[m][6] = getIndex(i,j+1,k+1);
                simplex[m][7] = getIndex(i+1,j+1,k+1);
                m++;
              }
        }
        else {
          for(int i =0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++) {
              simplex[m][0] = getIndex(i,j);
              simplex[m][1] = getIndex(i+1,j);
              simplex[m][2] = getIndex(i,j+1);
              simplex[m][3] = getIndex(i+1,j+1);
              m++;
            }
        }
        std::cout << "done" << std::endl;
        return simplex.size();
      }

      void getCubeBoundary(std::map<EntityKey<int>,int>& facemap) {

        //fill facemap
        std::cout << "Filling boundary facemap of Hexagrid...";
        if(dimw_ == 3) {
          for(int j=0; j < nofcells_[1]; j++)
            for(int k=0; k < nofcells_[2]; k++) {
              std::vector<int> entity(4);
              entity[0] = getIndex(0,j,k);
              entity[1] = getIndex(0,j+1,k);
              entity[2] = getIndex(0,j+1,k+1);
              entity[3] = getIndex(0,j,k+1);
              EntityKey<int> key3(entity);
              facemap[key3] = 0;
              entity[0] = getIndex(nofcells_[0],j,k);
              entity[1] = getIndex(nofcells_[0],j,k+1);
              entity[2] = getIndex(nofcells_[0],j+1,k+1);
              entity[3] = getIndex(nofcells_[0],j+1,k);
              EntityKey<int> key4(entity);
              facemap[key4] = 0;
            }
          for(int i=0; i < nofcells_[0]; i++)
            for(int k=0; k < nofcells_[2]; k++) {
              std::vector<int> entity(4);
              entity[0] = getIndex(i,0,k);
              entity[1] = getIndex(i,0,k+1);
              entity[2] = getIndex(i+1,0,k+1);
              entity[3] = getIndex(i+1,0,k);
              EntityKey<int> key3(entity);
              facemap[key3] = 0;
              entity[0] = getIndex(i,nofcells_[1],k);
              entity[1] = getIndex(i+1,nofcells_[1],k);
              entity[2] = getIndex(i+1,nofcells_[1],k+1);
              entity[3] = getIndex(i,nofcells_[1],k+1);
              EntityKey<int> key4(entity);
              facemap[key4] = 0;
            }
          for(int i=0; i < nofcells_[0]; i++)
            for(int j=0; j < nofcells_[1]; j++) {
              std::vector<int> entity(4);
              entity[0] = getIndex(i,j,0);
              entity[1] = getIndex(i+1,j,0);
              entity[2] = getIndex(i+1,j+1,0);
              entity[3] = getIndex(i,j+1,0);
              EntityKey<int> key3(entity);
              facemap[key3] = 0;
              entity[0] = getIndex(i,j,nofcells_[2]);
              entity[1] = getIndex(i,j+1,nofcells_[2]);
              entity[2] = getIndex(i+1,j+1,nofcells_[2]);
              entity[3] = getIndex(i+1,j,nofcells_[2]);
              EntityKey<int> key4(entity);
              facemap[key4] = 0;
            }
        }
        else {
          for(int i=0; i < nofcells_[0]; i++) {
            std::vector<int> entity(2);
            entity[0] = getIndex(i,0);
            entity[1] = getIndex(i+1,0);
            EntityKey<int> key3(entity);
            facemap[key3] = 0;
            entity[0] = getIndex(i+1,nofcells_[1]);
            entity[1] = getIndex(i,nofcells_[1]);
            EntityKey<int> key4(entity);
            facemap[key4] = 0;
          }
          for(int j=0; j < nofcells_[1]; j++) {
            std::vector<int> entity(2);
            entity[0] = getIndex(0,j);
            entity[1] = getIndex(0,j+1);
            EntityKey<int> key3(entity);
            facemap[key3] = 0;
            entity[0] = getIndex(nofcells_[0],j+1);
            entity[1] = getIndex(nofcells_[0],j);
            EntityKey<int> key4(entity);
            facemap[key4] = 0;
          }
        }
        std::cout << "done" << std::endl;
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

      bool ok() {
        return good_;
      }

      int dimw() {
        return dimw_;
      }

      int getIndex(int i,int j, int k = 0) {
        if(dimw_ == 3)
          return i*(nofcells_[1]+1)*(nofcells_[2]+1) + j*(nofcells_[2]+1) + k;
        else
          return i*(nofcells_[1]+1) + j;
      }
    };
    const char* IntervalBlock::ID = "Interval";
  }
}
#endif
