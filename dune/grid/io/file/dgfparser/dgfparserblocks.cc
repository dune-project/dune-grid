// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/dgfparserblocks.hh>

namespace Dune
{

  namespace dgf
  {

    // BasicBlock
    // ----------

    BasicBlock :: BasicBlock ( std :: istream &in, const char *id )
      : pos(-1),
        active(false),
        empty(true),
        identifier(id),
        linecount(0)
    {
      makeupcase( identifier );
      in.clear();
      in.seekg(0);
      if (!in) {
        DUNE_THROW(DGFException,
                   "file not found in BasicBlock::BasicBlock");
      }
      getblock( in );
      empty = (linecount == 0);
      if (active && !empty) {
        //linecount=countlines();
        reset();
      }
      in.clear();
      in.seekg(0);
    }


    void BasicBlock :: getblock ( std :: istream &in )
    {
      linecount = 0;
      while( in.good() )
      {
        std :: string line;
        getline( in, line );

        std :: istringstream linestream( line );
        std :: string id;
        linestream >> id;

        makeupcase( id );
        if( id == identifier )
          break;
      }
      if( in.eof() )
        return;

      active = true;
      while( in.good() )
      {
        std :: string line;
        getline( in, line );

        // strip comments
        if( !line.empty() )
        {
          std::size_t comment = line.find( '%' );
          if( comment != std::string::npos )
            line.erase( comment );
        }
        if( line.empty() )
          continue;

        std :: istringstream linestream( line );
        char test = 0;
        linestream >> test;
        if( test == '#' )
          return;

        ++linecount;
        block << line << std :: endl;
      }
      DUNE_THROW( DGFException, "Error reading from stream." );
    }


    // get next line and store in string stream
    bool BasicBlock :: getnextline ()
    {
      getline( block, oneline );
      line.clear();
      line.str( oneline );
      ++pos;
      return !oneline.empty();
    }


    bool BasicBlock :: gettokenparam ( std :: string token, std :: string &entry )
    {
      reset();
      makeupcase( token );
      while( getnextline() )
      {
        std :: string ltoken;
        line >> ltoken;
        makeupcase( ltoken );
        if( ltoken == token )
        {
          getline( line, entry );
          return true;
        }
      }
      return false;
    }


    bool BasicBlock :: findtoken ( std :: string token )
    {
      reset();
      makeupcase( token );
      while( getnextline() )
      {
        std :: string ltoken;
        line >> ltoken;
        makeupcase( ltoken );
        if( ltoken == token )
          return true;
      }
      return false;
    }



    // VertexBlock
    // -----------

    const char *VertexBlock :: ID = "Vertex";


    VertexBlock :: VertexBlock ( std :: istream &in, int &pdimworld )
      : BasicBlock( in, ID ),
        dimvertex( -1 ),
        dimworld( pdimworld ),
        goodline( true ),
        vtxoffset( 0 ),
        nofParam( 0 )
    {
      if (!isactive())
        return;

      if( findtoken( "firstindex" ) )
      {
        int x;
        if( getnextentry( x ) )
          vtxoffset=x;
      }

      if( findtoken( "parameters" ) )
      {
        int x;
        if( getnextentry( x ) )
          nofParam=x;
      }

      dimvertex = getDimWorld();
      if( pdimworld < 0 )
        pdimworld = dimvertex;
      dimworld = pdimworld;

      if( dimworld < dimvertex )
      {
        DUNE_THROW( DGFException,
                    "Error in " << *this << ": "
                                << "Vertex dimension greater than world dimension." );
      }
      if( dimworld > dimvertex )
      {
        dwarn << ID << " block: Embedding "
              << dimvertex << "-dimensional vertices into "
              << dimworld << "-dimensional space." << std::endl;
      }
    }


    int VertexBlock :: get ( std :: vector< std :: vector< double > > &points,
                             std :: vector< std :: vector< double > > &params,
                             int &nofp )
    {
      nofp = nofParam;
      reset();

      std::vector< double > point( dimworld );
      std::vector< double > param( nofParam );
      while( next( point, param ) )
      {
        points.push_back( point );
        if( nofParam > 0 )
          params.push_back( param );
      }
      return points.size();
    }


    int VertexBlock :: getDimWorld ()
    {
      if( findtoken( "dimension" ) )
      {
        int dimworld;
        if( !getnextentry( dimworld ) || (dimworld <= 0) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid value given for 'dimension'." );
        }
        return dimworld;
      }

      reset();
      while( getnextline() )
      {
        int dimworld = -nofParam;
        double x;
        while( getnextentry( x ) )
          ++dimworld;
        if( dimworld > 0 )
          return dimworld;
      }

      DUNE_THROW( DGFException,
                  "Error in " << *this << ": "
                              << "Unable to determine dimension of vertices." );
    }


    bool VertexBlock :: next ( std :: vector< double > &point,
                               std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      int n = 0;
      double x;
      for( ; getnextentry( x ); ++n )
      {
        if( n < dimvertex )
          point[ n ] = x;
        else if( n-dimvertex < nofParam )
          param[ n-dimvertex ] = x;
      }

      if( n == 0 )
        return next( point, param );
      else if( n != dimvertex + nofParam )
      {
        DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                               << "Wrong number of coordinates and parameters "
                                               << "(got " << n
                                               << ", expected " << (dimvertex + nofParam) << ")" );
      }

      for( int i = dimvertex; i < dimworld; ++i )
        point[ i ] = double( 0 );
      return (goodline = true);
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
        dumpfilename_(),
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
      if (findtoken("dumpfilename"))
        if (getnextentry(p)) {
          dumpfilename_=p;
        }
    }


    // SimplexBlock
    // ------------

    const char *SimplexBlock :: ID = "Simplex";


    SimplexBlock :: SimplexBlock
      ( std :: istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid )
      : BasicBlock( in, ID ),
        nofvtx( pnofvtx ),
        vtxoffset( pvtxoffset ),
        dimgrid( pdimgrid ),
        goodline( true ),
        nofparams( 0 )
    {
      if( !isactive() )
        return;

      if( findtoken( "parameters" ) )
      {
        int x = 0;
        if( getnextentry( x ) )
        {
          if( x > 0 )
            nofparams = x;
        }
        if( x <= 0 )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Key 'parameters' found with no or non-positive value." );
        }
      }

      if( dimgrid < 0 )
        dimgrid = getDimGrid();
      pdimgrid = dimgrid;
    }


    int SimplexBlock :: getDimGrid ()
    {
      reset();
      while( getnextline() )
      {
        int count = 0;
        double x;
        while( getnextentry( x ) )
          ++count;
        if( count > nofparams )
          return (count - nofparams) - 1;
      }
      return 0;
    }


    int SimplexBlock
    :: get ( std :: vector< std :: vector< unsigned int > > &simplices,
             std :: vector< std :: vector< double > > &params,
             int &nofp)
    {
      nofp = nofparams;
      reset();

      std :: vector< unsigned int > simplex( dimgrid+1 );
      std :: vector< double > param( nofparams );
      int nofsimpl = 0;
      for( ; next( simplex, param ); ++nofsimpl )
      {
        simplices.push_back( simplex );
        /*
           for( size_t j = 0; j < simplex.size(); ++j )
           simplices[ nofsimpl ][ j ] = simplex[ j ];
         */
        if( nofparams > 0 )
          params.push_back( param );
      }
      return nofsimpl;
    }


    bool SimplexBlock :: next ( std :: vector< unsigned int > &simplex,
                                std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      for( std :: size_t n = 0; n < simplex.size(); ++n )
      {
        int idx;
        if( !getnextentry( idx ) )
        {
          if( n > 0 )
          {
            DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                                   << "Wrong number of vertex indices "
                                                   << "(got " << idx
                                                   << ", expected " << simplex.size() << ")" );
          }
          else
            return next( simplex, param );
        }
        if( (vtxoffset > idx) || (idx >= int(nofvtx + vtxoffset)) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid vertex index "
                                  << "(" << idx << " not in ["
                                  << vtxoffset << ", " << (nofvtx + vtxoffset) << "[)" );
        }
        simplex[ n ] = idx - vtxoffset;
      }

      std :: size_t np = 0;
      double x;
      for( ; getnextentry( x ); ++np )
      {
        if( np < param.size() )
          param[ np ] = x;
      }

      if( np != param.size() )
      {
        DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                               << "Wrong number of simplex parameters "
                                               << "(got " << np
                                               << ", expected " << param.size() << ")" );
      }
      return (goodline = true);
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
      if (elements.empty())
        return 0;

      const int dimworld = vtx[ 0 ].size();

      int dimgrid = 0;
      for( size_t n = elements[ 0 ].size(); n > 1; ++dimgrid, n /= 2 ) ;
      if( size_t( 1 << dimgrid ) != elements[ 0 ].size() )
        DUNE_THROW( DGFException, "cube2simplex: all elements must be cubes." );

      dverb << "generating simplices...";
      dverb.flush();

      if( dimgrid == 1 )
        return elements.size();

      std::vector< std::vector< unsigned int > > cubes;
      std::vector< std::vector< double > > cubeparams;
      elements.swap( cubes );
      params.swap( cubeparams );

      if( dimgrid == 3 )
      {
        elements.resize( 6*cubes.size() );
        if( cubeparams.size() > 0 )
          params.resize( 6*cubes.size() );
        for( size_t countsimpl = 0; countsimpl < elements.size(); ++countsimpl )
          elements[ countsimpl ].resize( 4 );
        for( size_t c = 0; c < cubes.size(); ++c )
        {
          for( int tetra = 0; tetra < 6; ++tetra )
          {
            for( int v = 0; v < 4; ++v )
            {
              elements[ c*6+tetra ][ v ]
                = cubes[ c ][ offset3[ tetra ][ v ][ 0 ] +2*offset3[ tetra ][ v ][ 1 ] +4*offset3[ tetra ][ v ][ 2 ] ];
            }
            if( cubeparams.size() > 0 )
              params[ c*6+tetra ] = cubeparams[ c ];
          }
        }
      }
      else if( dimgrid == 2 )
      {
        elements.resize( 2*cubes.size() );
        if( cubeparams.size() > 0 )
          params.resize( 2*cubes.size() );
        for( size_t countsimpl = 0; countsimpl < elements.size(); ++countsimpl )
          elements[ countsimpl ].resize( 3 );
        for( size_t c = 0; c < cubes.size(); ++c )
        {
          int diag = 0;
          double mind = 0;
          for( int d = 0; d < 2; ++d )
          {
            // let's cut the longer diagonal
            double diaglen = 0;
            for( int i = 0; i < dimworld; ++i )
            {
              const double dist = vtx[ cubes[ c ][ d ] ][ i ] - vtx[ cubes[ c ][ 3-d ] ][ i ];
              diaglen += dist*dist;
            }
            if( diaglen < mind )
            {
              mind = diaglen;
              diag = d;
            }
          }
          if( diag == 0 )
          {
            int tmp0 = cubes[ c ][ 0 ];
            cubes[ c ][ 0 ] = cubes[ c ][ 1 ];
            cubes[ c ][ 1 ] = cubes[ c ][ 3 ];
            cubes[ c ][ 3 ] = cubes[ c ][ 2 ];
            cubes[ c ][ 2 ] = tmp0;
          }

          for( int triangle = 0; triangle < 2; ++triangle )
          {
            for( int v = 0; v < 3; ++v )
            {
              elements[ c*2+triangle ][ v ]
                = cubes[ c ][ offset2[ triangle ][ v ][ 0 ] + 2*offset2[ triangle ][ v ][ 1 ] ];
            }
            if( cubeparams.size() > 0 )
              params[ c*2+triangle ] = cubeparams[ c ];
          }
        }
      }
      return elements.size();
    }



    // CubeBlock
    // ---------

    const char *CubeBlock :: ID = "Cube";

    CubeBlock :: CubeBlock
      ( std::istream &in, int pnofvtx, int pvtxoffset, int &pdimgrid )
      : BasicBlock( in, ID ),
        nofvtx(pnofvtx),
        dimgrid( pdimgrid ),
        goodline(true),
        map(0),
        nofparams(0),
        vtxoffset(pvtxoffset)
    {
      if( !isactive() )
        return;

      if( findtoken( "parameters" ) )
      {
        int x = 0;
        if( getnextentry( x ) )
        {
          if( x > 0 )
            nofparams = x;
        }
        if( x <= 0 )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Key 'parameters' found with no or non-positive value." );
        }
      }

      if( dimgrid < 0 )
        dimgrid = getDimGrid();
      pdimgrid = dimgrid;

      map.resize( 1 << dimgrid );
      for( size_t i = 0; i < map.size(); ++i )
        map[ i ] = i;

      if( findtoken( "map" ) )
      {
        for( size_t i = 0; i < map.size(); ++i )
        {
          int x;
          if( !getnextentry( x ) )
          {
            DUNE_THROW( DGFException,
                        "Error in " << *this << ": "
                                    << "Incomplete reference mapping "
                                    << "(got " << i << " entries, "
                                    << "expected " << map.size() << " entries." );
          }
          map[ i ] = x;
        }
      }
    }


    int CubeBlock :: getDimGrid ()
    {
      reset();
      while( getnextline() )
      {
        int count = 0;
        double x;
        while( getnextentry( x ) )
          ++count;
        if( count > nofparams )
        {
          count -= nofparams;
          // int dim = (int)(log( count ) / M_LN2);
          int dim = 1;
          while (1<<dim < count)
            dim++;
          if( (dim < 0) || ((1 << dim) != count) )
          {
            DUNE_THROW( DGFException,
                        "Error in " << *this << ": Number of vertex indices ("
                                    << count << ") is not a power of 2." );
          }
          return dim;
        }
      }
      return 0;
    }


    int CubeBlock :: get ( std :: vector< std :: vector< unsigned int> > &cubes,
                           std :: vector< std :: vector< double > > &params,
                           int &nofp )
    {
      nofp = nofparams;
      reset();

      std :: vector< unsigned int > cube( 1 << dimgrid );
      std :: vector< double > param( nofparams );
      int nofcubes = 0;
      for( ; next( cube, param ); ++nofcubes )
      {
        cubes.push_back( cube );
        /*
           for( size_t j = 0; j < cube.size(); ++j )
           cubes[ nofcubes ][ j ] = cubes[ j ];
         */
        if( nofparams > 0 )
          params.push_back( param );
      }
      return nofcubes;
    }


    bool CubeBlock :: next ( std :: vector< unsigned int > &cube,
                             std :: vector< double > &param )
    {
      assert( ok() );
      if( !getnextline() )
        return (goodline = false);

      for( std :: size_t n = 0; n < cube.size(); ++n )
      {
        int idx;
        if( !getnextentry( idx ) )
        {
          if( n > 0 )
          {
            DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                                   << "Wrong number of vertex indices "
                                                   << "(got " << idx
                                                   << ", expected " << cube.size() << ")" );
          }
          else
            return next( cube, param );
        }
        if( (vtxoffset > idx) || (idx >= int(nofvtx + vtxoffset)) )
        {
          DUNE_THROW( DGFException,
                      "Error in " << *this << ": "
                                  << "Invalid vertex index "
                                  << "(" << idx << " not in ["
                                  << vtxoffset << ", " << (nofvtx + vtxoffset) << "[)" );
        }
        cube[ map[ n ] ] = idx - vtxoffset;
      }

      std :: size_t np = 0;
      double x;
      for( ; getnextentry( x ); ++np )
      {
        if( np < param.size() )
          param[ np ] = x;
      }

      if( np != param.size() )
      {
        DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                               << "Wrong number of simplex parameters "
                                               << "(got " << np
                                               << ", expected " << param.size() << ")" );
      }
      return (goodline = true);
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


    GridParameterBlock::GridParameterBlock ( std :: istream &in )
      : BasicBlock( in, ID ),
        foundFlags_( 0 ),
        name_( "Unnamed Grid" ), // default value (used if name is empty)
        dumpFileName_( "" ),
        markLongestEdge_( false ) // default value
    {
      if( isempty() )
        return;

      // check name
      if( findtoken( "name" ) )
      {
        std::string entry;
        if( getnextentry( entry ) )
          name_ = entry;
        else
          dwarn << "GridParameterBlock: Found keyword 'name' without value." << std::endl;
        foundFlags_ |= foundName;
      }

      // get file name
      if ( findtoken( "dumpfilename" ) )
      {
        std::string entry;
        if( getnextentry( entry ) )
          dumpFileName_ = entry;
        else
          dwarn << "GridParameterBlock: Found keyword 'dumpFileName' without value." << std::endl;
        foundFlags_ |= foundDumpFileName;
      }

      // check for markLongestEdge
      if( findtoken( "refinementedge" ) )
      {
        std::string entry;
        if( getnextentry( entry ) )
        {
          makeupcase( entry );
          if( entry == "LONGEST" )
            markLongestEdge_ = true;
          else if( entry != "ARBITRARY" )
            dwarn << "GridParameterBlock: Invalid value for keyword 'refinementedge': " << entry << std::endl;
        }
        else
          dwarn << "GridParameterBlock: Found keyword 'refinementedge' without value." << std::endl;
        foundFlags_ |= foundLongestEdge;
      }

    }


    // IntervalBlock
    // -------------

    const char *IntervalBlock::ID = "Interval";


    IntervalBlock::IntervalBlock ( std::istream &in )
      : BasicBlock( in, ID ),
        intervals_( 0 ),
        good_( false ),
        dimw_( 0 )
    {
      if( !isactive() )
        return;

      getnextline();
      for( double x; getnextentry( x ); ++dimw_ ) ;
      if( dimw_ == 0 )
      {
        DUNE_THROW( DGFException,
                    "Too few coordinates for point p0 in IntervalBlock" );
      }

      reset();
      while( next() ) ;
    }


    int IntervalBlock::getVtx ( int block, std::vector< std::vector< double > > &vtx ) const
    {
      dverb << "reading vertices for interval " << block << "... ";

      const Interval &interval = get( block );

      size_t old_size = vtx.size();
      vtx.resize( old_size + nofvtx( block ) );
      for( size_t i = old_size; i < vtx.size(); ++i )
        vtx[ i ].resize( dimw() );

      size_t m = old_size;
      std::vector< int > i( dimw() );
      const int end = dimw()-1;
      int k = end;
      for( i[ end ] = 0; i[ end ] <= interval.n[ end ]; )
      {
        // go all the way down
        for( ; k > 0; --k )
          i[ k-1 ] = 0;

        assert( m < vtx.size() );
        for( int j = 0; j < dimw(); ++j ) {
          vtx[ m ][ j ] = interval.p[ 0 ][ j ] + double(i[ j ])*interval.h[ j ];
        }
        ++m;

        // increase i[ k ] and go up for all finished loops
        for( ; (++i[ k ] > interval.n[ k ]) && (k < end); ++k ) ;
      }
      assert( m == vtx.size() );

      dverb << "[done]" << std::endl;
      return m - old_size;;
    }


    int IntervalBlock::getHexa ( int block, std::vector< std::vector< unsigned int > > &cubes, int offset ) const
    {
      dverb << "generating cubes for interval " << block << "... ";

      const Interval &interval = get( block );

      const int verticesPerCube = 1 << dimw();

      size_t old_size = cubes.size();
      cubes.resize( old_size + nofhexa( block ) );
      for( size_t i = old_size; i < cubes.size(); ++i )
        cubes[ i ].resize( verticesPerCube );

      size_t m = old_size;
      std::vector< int > i( dimw() );
      const int end = dimw()-1;
      int k = end;
      for( i[ end ] = 0; i[ end ] < interval.n[ end ]; )
      {
        // go all the way down
        for( ; k > 0; --k )
          i[ k-1 ] = 0;

        assert( m < cubes.size() );
        for( int j = 0; j < verticesPerCube; ++j )
        {
          cubes[ m ][ j ] = offset;
          int factor = 1;
          for( int d = 0; d < dimw(); ++d )
          {
            cubes[ m ][ j ] += factor*(i[ d ] + ((j >> d) & 1));
            factor *= interval.n[ d ]+1;
          }
        }
        ++m;

        // increase i[ k ] and go up for all finished loops
        for( ; (++i[ k ] >= interval.n[ k ]) && (k < end); ++k ) ;
      }
      assert( m == cubes.size() );

      dverb << "[done]" << std::endl;
      return m - old_size;
    }


    template< class T >
    void IntervalBlock::parseLine ( std::vector< T > &v )
    {
      getnextline();
      v.resize( dimw_ );
      for( int i = 0; i < dimw_; ++i )
      {
        if( !getnextentry( v[ i ] ) )
          DUNE_THROW( DGFException, "ERROR in " << *this << ": Not enough values." );
      }
    }


    bool IntervalBlock::next ()
    {
      if (linenumber()==noflines()-1)
      {
        good_=false;
        return good_;
      }

      Interval interval;
      parseLine( interval.p[ 0 ] );
      parseLine( interval.p[ 1 ] );
      parseLine( interval.n );

      //find real upper and lower edge and calculate cell width
      interval.h.resize( dimw_ );
      for( int i = 0; i < dimw_; ++i )
      {
        double &left = interval.p[ 0 ][ i ];
        double &right = interval.p[ 1 ][ i ];
        const int &n = interval.n[ i ];

        if( left > right )
        {
          double dummy = left;
          left = right;
          right = dummy;
        }

        interval.h[ i ] = (right - left) / double( n );
        assert( interval.h[ i ] > 0);
      }
      intervals_.push_back( interval );

      dverb << interval << std::endl;

      good_ = true;
      return good_;
    }



    // PeriodicFaceTransformationBlock
    // -------------------------------

    const char *PeriodicFaceTransformationBlock::ID = "PeriodicFaceTransformation";

    PeriodicFaceTransformationBlock
    ::PeriodicFaceTransformationBlock ( std::istream &in, int dimworld )
      : BasicBlock( in, ID )
    {
      while( getnextline() )
      {
        AffineTransformation trafo( dimworld );
        for( int i = 0; i < dimworld; ++i )
        {
          if( i > 0 )
            match( ',' );

          for( int j = 0; j < dimworld; ++j )
          {
            if( !getnextentry( trafo.matrix( i, j ) ) )
            {
              DUNE_THROW( DGFException,
                          "Error in " << *this << ": "
                                      << "Not enough entries in matrix row " << i << "." );
            }
          }
        }

        match( '+' );
        for( int i = 0; i < dimworld; ++i )
        {
          if( !getnextentry( trafo.shift[ i ] ) )
          {
            DUNE_THROW( DGFException,
                        "Error in " << *this << ": "
                                    << "Not enough entries in shift." );
          }
        }

        transformations_.push_back( trafo );
      }
    }


    void PeriodicFaceTransformationBlock::match ( char what )
    {
      char c;
      if( !getnextentry( c ) || (c != what) )
      {
        DUNE_THROW( DGFException,
                    "Error in " << *this << ": " << what << "expected." );
      }
    }

  } // end namespace dgf

} // end namespace Dune
