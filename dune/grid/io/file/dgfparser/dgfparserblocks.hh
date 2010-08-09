// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSERBLOCKS_HH
#define DUNE_MACROGRIDPARSERBLOCKS_HH

#include <cassert>
#include <iostream>
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
      int dimvertex;       // the dimension of the vertices (determined from DGF file)
      int dimworld;        // the dimension of the world (either dimvertex or given by user)
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
      std::string dumpfilename_;
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

      const std::string dumpFileName ( ) const
      {
        return dumpfilename_;
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


    //-  *************************************************************

    /** \brief Common Grid parameters
        \ingroup DGFGridParameter

         For each grid implementation there is a set of parameters
         that can be passed via the GridParameter block to the momment of
         grid construction.
         Currently implemented common parameters are: \n\n
         1. \b name: The name of the grid ( later returned by the method grid.name() ). \n
         2. \b refinementedge: parameter to specify the refinement edge in simplices.
            Valid values are \b arbitrary (which is the default value) and \b longest
            which marks the longest edge/face of each simplex to be the refinement edge. \n
        See also the \b examplegrid5.dgf file for examples.
     */
    class GridParameterBlock
      : public BasicBlock
    {
    public:
      typedef unsigned int Flags;

      static const Flags foundName = 1 << 0;
      static const Flags foundDumpFileName = 1 << 1;
      static const Flags foundLongestEdge = 1 << 5;

    protected:
      Flags foundFlags_; // supportFlags, this block was created with
      std::string name_; // name of the grid
      std::string dumpFileName_; // name of the grid
      bool markLongestEdge_; // Mark longest edge for AlbertaGrid

    private:
      // copy not implemented
      GridParameterBlock(const GridParameterBlock&);

    public:
      const static char* ID;
      //! constructor: read commmon parameters
      GridParameterBlock ( std::istream &in );

      //! return the name of the grid
      const std::string &name ( const std::string &defaultValue ) const
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

      const std::string &dumpFileName ( ) const
      {
        if( (foundFlags_ & foundDumpFileName) != 0 )
        {
          dwarn << "GridParameterBlock: found Parameter 'dumpfilename', "
                << "dumping file to `" << dumpFileName_ << "'" << std::endl;
        }
        return dumpFileName_;
      }

      //! returns true if longest edge should be marked for AlbertaGrid
      bool markLongestEdge () const
      {
        if( (foundFlags_ & foundLongestEdge) == 0 )
        {
          dwarn << "GridParameterBlock: Parameter 'refinementedge' not specified, "
                << "defaulting to 'ARBITRARY'." << std::endl;
        }
        return markLongestEdge_;
      }

      // some information
      bool ok()
      {
        return true;
      }
    };



    // IntervalBlock
    // -------------

    struct IntervalBlock
      : public BasicBlock
    {
      struct Interval
      {
        std::vector< double > p[ 2 ]; // lower and upper boundary points
        std::vector< double > h;      // width of the cells in each direction
        std::vector< int > n;         // number of cells in each direction
      };

    private:
      std::vector< Interval > intervals_;
      bool good_;                      //data read correctly
      int dimw_;                       //dimension of world

    public:
      const static char* ID;
      explicit IntervalBlock ( std::istream &in );

      void get ( std::vector< std::vector< double > > &vtx, int &nofvtx,
                 std::vector< std::vector< unsigned int > > &simplex, int &nofsimpl )
      {
        for( size_t i = 0; i < intervals_.size(); ++i )
        {
          int oldvtx = nofvtx;
          nofvtx += getVtx( i, vtx );
          nofsimpl += getHexa( i, simplex, oldvtx );
        }
      }

      void get ( std::vector< std::vector< double > > &vtx, int &nofvtx )
      {
        for( size_t i = 0; i < intervals_.size(); ++i )
          nofvtx += getVtx( i, vtx );
      }

      const Interval &get ( int block ) const
      {
        return intervals_[ block ];
      }

      int numIntervals () const
      {
        return intervals_.size();
      }

      int dimw () const
      {
        return dimw_;
      }

      int getVtx ( int block, std::vector< std::vector< double > > &vtx ) const;
      int getHexa ( int block, std::vector< std::vector< unsigned int > > &cubes,
                    int offset = 0 ) const;

      int nofvtx ( int block ) const
      {
        const Interval &interval = get( block );
        int n = 1;
        for( int i = 0; i < dimw_; ++i )
          n *= (interval.n[ i ] + 1);
        return n;
      }

      int nofhexa ( int block ) const
      {
        const Interval &interval = get( block );
        int n = 1;
        for( int i = 0; i < dimw_; ++i )
          n *= interval.n[ i ];
        return n;
      }

    private:
      template< class T >
      void parseLine ( std::vector< T > &v );

      bool next ();
    };


    inline std::ostream &
    operator<< ( std::ostream &out, const IntervalBlock::Interval &interval )
    {
      if( interval.p[ 0 ].empty() || interval.p[ 1 ].empty() || interval.n.empty() )
        return out << "Interval {}";

      out << "Interval { p0 = (" << interval.p[ 0 ][ 0 ];
      for( size_t i = 1; i < interval.p[ 0 ].size(); ++i )
        out << ", " << interval.p[ 0 ][ i ];
      out << "), p1 = (" << interval.p[ 1 ][ 0 ];
      for( size_t i = 1; i < interval.p[ 1 ].size(); ++i )
        out << ", " << interval.p[ 1 ][ i ];
      out << "), n = (" << interval.n[ 0 ];
      for( size_t i = 1; i < interval.n.size(); ++i )
        out << ", " << interval.n[ i ];
      return out << ") }";
    }



    // PeriodicFaceTransformationBlock
    // -------------------------------

    struct PeriodicFaceTransformationBlock
      : public BasicBlock
    {
      template< class T >
      class Matrix;

      struct AffineTransformation;

    private:
      std::vector< AffineTransformation > transformations_;

      // copy not implemented
      PeriodicFaceTransformationBlock ( const PeriodicFaceTransformationBlock & );

    public:
      static const char *ID;

      // initialize block and get dimension of world
      PeriodicFaceTransformationBlock ( std::istream &in, int dimworld );

      const AffineTransformation &transformation ( int i ) const
      {
        assert( i < numTransformations() );
        return transformations_[ i ];
      }

      int numTransformations () const
      {
        return transformations_.size();
      }

    private:
      void match ( char what );
    };



    // PeriodicFaceTransformationBlock::Matrix
    // ---------------------------------------

    template< class T >
    class PeriodicFaceTransformationBlock::Matrix
    {
      int rows_;
      int cols_;
      std::vector< T > fields_;

    public:
      Matrix ( int rows, int cols )
        : rows_( rows ),
          cols_( cols ),
          fields_( rows * cols )
      {}

      const T &operator() ( int i, int j ) const
      {
        return fields_[ i * cols_ + j ];
      }

      T &operator() ( int i, int j )
      {
        return fields_[ i * cols_ + j ];
      }

      int rows () const
      {
        return rows_;
      }

      int cols () const
      {
        return cols_;
      }
    };


    // PeriodicFaceTransformationBlock::AffineTransformation
    // -----------------------------------------------------

    struct PeriodicFaceTransformationBlock::AffineTransformation
    {
      Matrix< double > matrix;
      std::vector< double > shift;

      explicit AffineTransformation ( int dimworld )
        : matrix( dimworld, dimworld ),
          shift( dimworld )
      {}
    };



    inline std::ostream &
    operator<< ( std::ostream &out, const PeriodicFaceTransformationBlock::AffineTransformation &trafo )
    {
      for( int i = 0; i < trafo.matrix.rows(); ++i )
      {
        out << (i > 0 ? ", " : "");
        for( int j = 0; j < trafo.matrix.cols(); ++j )
          out << (j > 0 ? " " : "") << trafo.matrix( i, j );
      }
      out << " +";
      for( unsigned int i = 0; i < trafo.shift.size(); ++i )
        out << " " << trafo.shift[ i ];
      return out;
    }


  } // end namespace dgf

} // end namespace Dune

#endif
