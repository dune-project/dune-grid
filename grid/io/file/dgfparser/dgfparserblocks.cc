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
      empty = (linecount > 0);
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
          std :: size_t comment = line.find( "%" );
          if( comment != std :: string :: npos )
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

      if( dimworld < 0 )
        dimworld = getDimWorld();
      pdimworld=dimworld;

      if( dimworld <= 0 )
      {
        DUNE_THROW( DGFException,
                    "Error in " << *this << ": "
                                << "Unable to determine dimension of vertices." );
      }
    }


    int VertexBlock :: get ( std :: vector< std :: vector< double > > &points,
                             std :: vector< std :: vector< double > > &params,
                             int &nofp )
    {
      nofp = nofParam;
      reset();

      std :: vector< double > point( dimworld );
      std :: vector< double > param( nofParam );
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
      return 0;
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
        if( n < dimworld )
          point[ n ] = x;
        else if( n-dimworld < nofParam )
          param[ n-dimworld ] = x;
      }
      if( n == dimworld + nofParam )
        return (goodline = true);

      if( n > 0 )
      {
        DUNE_THROW ( DGFException, "Error in " << *this << ": "
                                               << "Wrong number of coordinates and parameters "
                                               << "(got " << n
                                               << ", expected " << (dimworld + nofParam) << ")" );
      }
      return next( point, param );
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


    GridParameterBlock::GridParameterBlock ( std :: istream &in, const bool readOverlapAndBnd )
      : BasicBlock( in, ID ),
        foundFlags_( 0 ),
        _periodic(),
        _overlap(0), // default value
        _noClosure(false), // default value
        _noCopy(true)    // default value
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

      if( readOverlapAndBnd )
      {
        // check overlap
        if( findtoken( "overlap" ) )
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
          foundFlags_ |= foundOverlap;
        }

        // check periodic grid
        if (findtoken("periodic"))
        {
          int x;
          while (getnextentry(x))
          {
            _periodic.insert(x);
          }
          foundFlags_ |= foundPeriodic;
        }
      }

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
        foundFlags_ |= foundClosure;
      }

      // check copies
      if( findtoken( "copies" ) )
      {
        std::string clop;
        if(getnextentry(clop))
        {
          makeupcase(clop);
          if(clop == "NONE")
          {
            _noCopy = true;
          }
          else {
            _noCopy = false;
          }
        }
        foundFlags_ |= foundCopies;
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
