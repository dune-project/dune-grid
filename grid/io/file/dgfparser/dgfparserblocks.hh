// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSERBLOCKS_HH
#define DUNE_MACROGRIDPARSERBLOCKS_HH

#include <cassert>
#include <map>
#include <set>
#include <vector>

#include <dune/common/stdstreams.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>

namespace Dune
{

  // *************************************************************
  // Read one block with given identifier from disk
  // and allows the line-wise extraction from this block
  // *************************************************************
  namespace dgf {

    inline void makeupcase( std :: string &s )
    {
      for (size_t i=0; i<s.size(); i++)
        s[i]=toupper(s[i]);
    }


    class BasicBlock
    {
      int pos;               // line number
      bool active;           // block was found
      bool empty;            // block was found but was empty
      std::string identifier; // identifier of this block
      int linecount;         // total number of lines in the block
      std::stringstream block; // the block itself
      std::string oneline;   // the active line in the block

      // get the block (if it exists)
      void getblock ( std::istream &in );
      // count the number of lines in the block
      // int countlines ();

    protected:
      std::stringstream line; // the active line as string buffer
                              // for use in the derived classes

      // go back to beginning of block
      void reset ()
      {
        pos = -1;
        block.clear();
        block.seekg( 0 );
      }

      // get next line and store in string stream
      bool getnextline ();

      // get next entry in line
      template< class ENTRY >
      bool getnextentry( ENTRY &entry )
      {
        line >> entry;
        return line;
      }

      bool gettokenparam ( std :: string token, std :: string &entry );
      bool findtoken( std :: string token );

    public:
      // search for block in file and store in buffer
      BasicBlock ( std::istream &in, const char* id );

      // some information on this block
      bool isactive ()
      {
        return active;
      }

      bool isempty ()
      {
        return empty;
      }

      int &noflines ()
      {
        return linecount;
      }

      int linenumber ()
      {
        return pos;
      }

      // for error messages
      friend std :: ostream &operator<< ( std :: ostream &os, const BasicBlock &b )
      {
        return os << "block " << b.identifier << " (line " << b.pos << ")";
      }
    };


    // *************************************************************
    // derived classes for each block in grid file
    // *************************************************************
    class VertexBlock
      : public BasicBlock
    {
      int dimworld;        // the dimesnsion of the vertices (is given from user)
      bool goodline;       // active line describes a vertex
      int vtxoffset;
      int nofParam;

    public:
      static const char* ID;

      // initialize vertex block
      VertexBlock ( std :: istream &in, int &pdimworld );

      int get ( std :: vector< std :: vector< double > > &vtx,
                std :: vector< std :: vector< double > > &param,
                int &nofp );

      // some information
      bool ok () const
      {
        return goodline;
      }

      int offset () const
      {
        return vtxoffset;
      }

    private:
      // get dimworld
      int getDimWorld ();

      // get next vertex
      bool next ( std :: vector< double > &point, std :: vector< double > &param );
    };


    // *************************************************************
    class SimplexGenerationBlock
      : public BasicBlock
    {
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
      SimplexGenerationBlock ( std :: istream &in );

      double maxArea ()
      {
        return area_;
      }

      double minAngle ()
      {
        return angle_;
      }

      bool display ()
      {
        return display_;
      }

      bool haspath ()
      {
        return haspath_;
      }

      std :: string path ()
      {
        return path_;
      }

      bool hasfile ()
      {
        return hasfile_;
      }

      std :: string filename ()
      {
        return filename_;
      }

      std :: string filetype ()
      {
        return filetype_;
      }

      int dimension ()
      {
        return dimension_;
      }

      std :: string parameter ()
      {
        return parameter_;
      }
    };


    // SimplexBlock
    // ------------

    class SimplexBlock
      : public BasicBlock
    {
      unsigned int nofvtx;
      int vtxoffset;
      int dimgrid;
      bool goodline;               // active line describes a vertex
      int nofparams;               // nof parameters

    public:
      const static char* ID;

      SimplexBlock ( std :: istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid );

      int get ( std :: vector< std :: vector< unsigned int > > &simplex,
                std :: vector< std :: vector< double > > &params,
                int &nofp );

      // cubes -> simplex
      static int
      cube2simplex ( std :: vector< std :: vector< double > > &vtx,
                     std :: vector< std :: vector< unsigned int > > &elements,
                     std :: vector< std :: vector< double > > &params );

      // some information
      bool ok ()
      {
        return goodline;
      }

      int nofsimplex ()
      {
        return noflines();
      }

    private:
      // get the dimension of the grid
      int getDimGrid ();
      // get next simplex
      bool next ( std :: vector< unsigned int > &simplex,
                  std :: vector< double > &param );
    };



    // CubeBlock
    // ---------

    class CubeBlock
      : public BasicBlock
    {
      unsigned int nofvtx;
      int dimgrid;
      bool goodline;      // active line describes a vertex
      std :: vector< unsigned int > map; // active vertex
      int nofparams;
      int vtxoffset;

    public:
      static const char* ID;

      CubeBlock ( std :: istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid );

      int get ( std :: vector< std :: vector< unsigned int> > &simplex,
                std :: vector< std :: vector< double > > &params,
                int &nofp );

      // some information
      bool ok ()
      {
        return goodline;
      }

      int nofsimplex ()
      {
        return noflines();
      }

    private:
      // get the dimension of the grid
      int getDimGrid ();
      // get next simplex
      bool next ( std :: vector< unsigned int > &simplex,
                  std :: vector< double > &param );
    };


    // *************************************************************
    // the block BoundaryDomBlock looks for a domain which is characterized by two points in R^dimworld
    class BoundaryDomBlock : public BasicBlock {
      int dimworld;    // the dimesnsion of the vertices (is given  from user)
      bool goodline;   // active line describes a vertex
      std::vector<double> p1,p2; // active vertex
      int bndid;
      bool withdefault;
      int defaultvalue;
    public:
      static const char* ID;
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


    // *************************************************************
    // the BoundarySegBlock looks for given boundary values unless they aren't given they got the value zero
    class BoundarySegBlock : public BasicBlock {
      int dimworld;                // the dimesnsion of the vertices (is given  from user)
      bool goodline;               // active line describes a vertex
      std :: vector< unsigned int > p; // active vertex
      int bndid;
      bool simplexgrid;
    public:
      static const char* ID;
      // initialize vertex block and get first vertex
      BoundarySegBlock ( std :: istream &in, int pnofvtx,
                         int pdimworld, bool psimplexgrid );

      // some information
      int get( std :: map< DGFEntityKey< unsigned int>, int > &facemap,
               bool fixedsize,
               int vtxoffset );
      bool ok() {
        return goodline;
      }
      int nofbound() {
        return noflines();
      }
    private:
      bool next();
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


    // *************************************************************
    class DimBlock : public BasicBlock {
      int _dimworld; // dimension of world
      int _dim;      // dimension of grid
    public:
      const static char* ID;
      // initialize block and get dimension of world
      DimBlock ( std :: istream &in );
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


    // *************************************************************
    class GridParameterBlock
      : public BasicBlock
    {
    public:
      typedef unsigned int Flags;

      static const Flags foundPeriodic = 1 << 1;
      static const Flags foundOverlap = 1 << 2;
      static const Flags foundClosure = 1 << 3;
      static const Flags foundCopies = 1 << 4;
      static const Flags foundName = 1 << 5;

    protected:
      Flags foundFlags_; // supportFlags, this block was created with
      std::set<int> _periodic; // periodic grid
      int _overlap; // overlap for YaspGrid
      bool _noClosure; // no closure for UGGrid
      bool _noCopy; // no copies for UGGrid
      std::string name_; // name of the grid

    private:
      // copy not implemented
      GridParameterBlock(const GridParameterBlock&);

    public:
      const static char* ID;
      // initialize block and get dimension of world
      GridParameterBlock ( std::istream &in, const bool readOverlapAndBnd = true );

      // get dimension of world found in block
      int overlap () const
      {
        if( (foundFlags_ & foundOverlap) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'overlap' not specified, "
                << "defaulting to '" << _overlap << "'." << std::endl;
        }
        return _overlap;
      }

      std::string name ( const std::string &defaultValue ) const
      {
        if( (foundFlags_ & foundName) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'name' not specified, "
                << "defaulting to '" << defaultValue << "'." << std::endl;
          return defaultValue;
        }
        else
          return name_;
      }

      // returns true if no closure should be used for UGGrid
      bool noClosure () const
      {
        if( (foundFlags_ & foundClosure) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'closure' not specified, "
                << "defaulting to 'GREEN'." << std::endl;
        }
        return _noClosure;
      }

      // returns true if no closure should be used for UGGrid
      bool noCopy () const
      {
        if( (foundFlags_ & foundCopies) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'copies' not specified, "
                << "no copies will be generated." << std::endl;
        }
        return _noCopy;
      }

      // returns true if dimension is periodic
      bool isPeriodic ( const int dim ) const
      {
        if( (foundFlags_ & foundPeriodic) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'copies' not specified, "
                << "defaulting to no periodic boundary." << std::endl;
        }
        return (_periodic.find(dim) != _periodic.end());
      }

      // some information
      bool ok()
      {
        return true;
      }
    };


    // *************************************************************
    class IntervalBlock : public BasicBlock {
      std::vector<double> p0_,p1_; //lower and upper boundary points
      std::vector<double> h_;      // width of the cells in every direction
      std::vector<int> nofcells_;  // number of cells in every direction
      bool good_;                  //data read correctly
      int dimw_;                   //dimension of world
    public:
      const static char* ID;
      IntervalBlock ( std :: istream &in );

      void get ( std::vector<std::vector<double> >& vtx,int& nofvtx,
                 std::vector<std::vector<unsigned int> >& simplex,int& nofsimpl )
      {
        do {
          int oldvtx = nofvtx;
          nofvtx  +=getVtx(vtx);
          nofsimpl+=getHexa(simplex,oldvtx);
        } while (next());
      }
      void get ( std::vector<std::vector<double> >& vtx,int& nofvtx )
      {
        do {
          // int oldvtx = nofvtx;
          nofvtx  +=getVtx(vtx);
        } while (next());
      }
      int getVtx(std::vector<std::vector<double> >& vtx);
      int getHexa ( std :: vector< std :: vector< unsigned int > > &simplex,
                    int offset = 0 );

      int nofvtx() {
        if(dimw_ == 3)
          return (nofcells_[0]+1)*(nofcells_[1]+1)*(nofcells_[2]+1);
        else if (dimw_ == 2)
          return (nofcells_[0]+1)*(nofcells_[1]+1);
        else
          return nofcells_[0]+1;
      }

      int nofhexa() {
        if(dimw_ == 3)
          return (nofcells_[0])*(nofcells_[1])*(nofcells_[2]);
        else if (dimw_ == 2)
          return (nofcells_[0])*(nofcells_[1]);
        else
          return nofcells_[0];
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

      int getIndex(int i,int j = 0, int k = 0)
      {
        if(dimw_ == 3)
          return k*(nofcells_[1]+1)*(nofcells_[0]+1) + j*(nofcells_[0]+1) + i;
        else if (dimw_ == 2)
          return j * (nofcells_[0]+1) + i;
        else
          return i;
      }
    private:
      bool next ();
    };

  } // end namespace dgf

} // end namespace Dune

#endif
