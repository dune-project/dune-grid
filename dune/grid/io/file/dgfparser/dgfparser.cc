// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cstdio>
#if HAVE_MKSTEMP
#include <unistd.h>
#endif

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/blocks/boundarydom.hh>

namespace Dune
{

  static const std::string dgfid( "DGF" );


  class DGFPrintInfo
  {
    std::ofstream out;

  public:
    DGFPrintInfo( const std :: string &name )
      : out( (name+".log").c_str() )
    {
      // out << "Grid generated from file " << name << std::endl;
      out << "DGF parser started" << std::endl;
    }

    void finish ()
    {
      out << "Successful" << std::endl;
    }

    template< class Block >
    void block( const Block& printBlock )
    {
      out << "Using " << printBlock.id() << " block" << std::endl;
    }

    void print( const std :: string &message )
    {
      out << message << std::endl;
    }

    void step1(int dimw,int nofvtx,int nofel)
    {
      out << "Step 1 finished: " << std::endl;
      out << "                 " << dimw << "d" << std::endl;
      out << "                 " << nofvtx << " vertices" << std::endl;
      out << "                 " << nofel << " elements" << std::endl;
    }

    void step2(int bndseg,int totalbndseg,
               int bnddomain,int defsegs,int remaining)
    {
      out << "Step 2 finished: " << std::endl;
      out << "                 " << bndseg
          << " bnd-segs read in BoundarySegment block" << std::endl;
      out << "                 " << totalbndseg
          << " bnd-segs in grid" << std::endl;
      out << "                 " << bnddomain
          << " bnd-segs in a boundary domain" << std::endl;
      out << "                 " << defsegs
          << " bnd-seg assigned default value" << std::endl;
      out << "                 " << remaining
          << " bnd-segs have not been assigned an id!" << std::endl;
    }

    void cube2simplex(DuneGridFormatParser::element_t el)
    {
      if (el == DuneGridFormatParser::General)
        out << "Simplex block found, thus converting "
            << "cube grid to simplex grid" << std::endl;
      else
        out << "Element type should be simplex, thus converting "
            << "cube grid to simplex grid" << std::endl;
    }

    void automatic() {
      out << "Automatic grid generation" << std::endl;
    }
  };



  // DuneGridFormatParser
  // --------------------

  DuneGridFormatParser :: DuneGridFormatParser ( int rank, int size )
    : dimw( -1 ),
      dimgrid( -1 ),
      vtx(0), nofvtx(0), vtxoffset(0), minVertexDistance(1e-12),
      elements(0) , nofelements(0),
      bound(0) , nofbound(0),
      facemap(),
      haveBndParameters( false ),
      element(General),
      simplexgrid(false),
      cube2simplex(false),
      nofvtxparams(0),
      nofelparams(0),
      vtxParams(),
      elParams(),
      info(0),
      rank_(rank),
      size_(size)
  {
    if( (rank_ < 0) || (rank_ >= size_) )
      DUNE_THROW( DGFException, "Invalid rank: " << rank << " (not in [ 0, " << size_ << " [)." );
  }


  // Output to Tetgen/Triangle poly-file
  void DuneGridFormatParser
  ::writeTetgenPoly ( const std::string &prefixname, std::string &extension, std::string &params )
  {
    std::string name = prefixname;
    params = "";
    if( dimw == 2 )
    {
      if (facemap.size()+elements.size()>0)
      {
        extension = ".poly";
        params += " -Ap ";
      }
      else
        extension = ".node";
      name += extension;
      info->print("writing poly file "+name);
      std::ofstream polys(name.c_str());
      writeTetgenPoly(polys);
    }
    else
    {
      if( facemap.size() > 0 && elements.size() == 0 )
      {
        extension = ".poly";
        name += extension;
        info->print("writing poly file "+name);
        std::ofstream polys(name.c_str());
        writeTetgenPoly(polys);
        params += " -p ";
      }
      else
      {
        extension = ".node";
        std::ofstream nodes( (name + extension).c_str() );
        writeTetgenPoly( nodes, false );

        {
          std::ofstream out( (name + ".ele").c_str() );
          dverb << "Writing elements...";
          out << elements.size() << " 4 " << nofelparams << std::endl;
          for (size_t n=0; n<elements.size(); n++) {
            out << n << "   ";
            for (int j=0; j<4; j++) {
              out << elements[n][j] << " ";
            }
            for (int j=0; j<nofelparams; ++j) {
              out << elParams[n][j] << " ";
            }
            out << std::endl;
          }
        }
        {
          std::ofstream out( (name + ".face").c_str() );
          facemap_t :: iterator pos;
          int nr = 0;
          dverb << "Writing boundary faces...";
          out << facemap.size() << " 1 " << std::endl;
          for(pos= facemap.begin(); pos!=facemap.end(); ++pos) {
            out << nr++ << " ";
            for (int i=0; i<pos->first.size(); i++)
              out << pos->first.origKey(i) << " ";
            out << pos->second.first;
            out << std::endl;
          }
        }

        if( elements.size() > 0 )
          params += " -r ";
      }
    }
  }


  // Output to Tetgen/Triangle poly-file
  void DuneGridFormatParser::writeTetgenPoly ( std::ostream &out, bool writeSegments )
  {
    dverb << "Writing vertices...";
    out << nofvtx << " " << dimw << " " << nofvtxparams << " 0" << std::endl;
    for( int n = 0; n < nofvtx; ++n )
    {
      out << n << " ";
      for( int j = 0; j < dimw; ++j )
        out << " " << vtx[ n ][ j ];
      for( int j = 0; j < nofvtxparams; ++j )
        out << " " << vtxParams[ n ][ j ];
      out << std::endl;
    }
    dverb << "done" << std::endl;

    if( writeSegments )
    {
      dverb << "Writing Segments...";
      out << facemap.size() + elements.size()*3 << " 1 " << std::endl;
      int nr = 0;

      // write out triangle faces (this code is purely 2d, but also used for tetgen)
      for( size_t i = 0; i < elements.size(); ++i )
      {
        for( int k = 0; k < 3; ++k )
          out << nr++ << " " << elements[ i ][ (k+1)%3 ] << " " << elements[ i ][ (k+2)%3 ] << " 0" << std::endl;
      }

      // write out boundary segments
      for( facemap_t::iterator pos = facemap.begin(); pos != facemap.end(); ++pos, ++nr )
      {
        if( dimw == 3 )
        {
          out << "1 0 " << pos->second.first << std::endl;
          out << pos->first.size();
        }
        else
          out << nr;
        for( int i = 0; i < pos->first.size(); ++i )
          out << " " << pos->first.origKey( i );
        if( dimw == 2 )
          out << " " << pos->second.first;
        out << std::endl;
      }
      out << "0" << std::endl;

      // write out element parameters (2d only)
      if( nofelparams > 0 )
      {
        if( dimw != 2 )
          DUNE_THROW( InvalidStateException, "Element parameters are not supported by tetgen." );
        out << elements.size() * nofelparams << std::endl;
        int outputLineNum = 0;
        for( size_t i = 0; i < elements.size(); ++i )
        {
          double coord[ 2 ] = { 0, 0 };
          for( int j = 0; j < 3; ++j )
          {
            coord[0] += vtx[ elements[ i ][ j ] ][ 0 ];
            coord[1] += vtx[ elements[ i ][ j ] ][ 1 ];
          }
          coord[0] /= 3.;
          coord[1] /= 3.;
          for( int j = 0; j < nofelparams; ++j )
            out << outputLineNum++ << " " << coord[ 0 ] << " " << coord[ 1 ] << " " << elParams[ i ][ j ] << std::endl;
        }
      }
      else
        out << 0 << std::endl;

      dverb << "done" << std::endl;
    }
  }


  bool DuneGridFormatParser::isDuneGridFormat ( std::istream &input )
  {
    input.clear();
    input.seekg( 0 );
    if( !input )
      return false;

    std::string idline;
    std::getline( input, idline );
    dgf::makeupcase( idline );

    std::string id;
    std::istringstream idstream( idline );
    idstream >> id;

    // compare id to DGF keyword
    return (id == dgfid);
  }


  bool DuneGridFormatParser::isDuneGridFormat ( const std::string &filename )
  {
    std::ifstream input( filename.c_str() );
    return isDuneGridFormat( input );
  }


  bool DuneGridFormatParser::readDuneGrid ( std::istream &gridin, int dimG, int dimW )
  {
    if( !isDuneGridFormat( gridin ) )
    {
      std::cerr << "Couldn't find '" << dgfid << "' keyword."
                << "File is not in DuneGridFormat. Exiting parser..."
                << std::endl;
      return false;
    } // not a DGF file, prehaps native file format

    // initialize variables
    cube2simplex = false;
    simplexgrid = false;
    dimgrid = dimG;
    dimw = dimW;
    vtxoffset = 0;
    nofvtx=0;
    nofelements=0;

    info = new DGFPrintInfo( "dgfparser" );

    dgf :: IntervalBlock interval( gridin );
    dgf :: VertexBlock bvtx( gridin, dimw );

    // generate cartesian grid?
    if ( interval.isactive() )
    {
      info->automatic();
      if( bvtx.isactive() )
      {
        nofvtx = bvtx.get( vtx, vtxParams, nofvtxparams );
        info->block( bvtx );
      }

      info->block( interval );
      if( dimw >= 0 )
      {
        if( dimw != interval.dimw() )
        {
          DUNE_THROW( DGFException,
                      "Error in " << interval << ": Wrong coordinate dimension "
                                  << "(got " << interval.dimw() << ", expected " << dimw << ")" );
        }
      }
      else
        dimw = interval.dimw();
      if( (dimgrid >= 0) && (dimgrid != dimw) )
      {
        DUNE_THROW( DGFException,
                    "Error in " << interval << ": Coordinate dimension differs "
                                << "from grid dimension." );
      }
      else
        dimgrid = dimw;
      cube2simplex = true;

      simplexgrid = (element == Simplex);
      if( element == General )
      {
        dgf :: SimplexBlock bsimplex( gridin, -1, -1, dimgrid );
        simplexgrid = bsimplex.isactive();
        if (simplexgrid)
          info->cube2simplex( element );
      }

      interval.get( vtx, nofvtx, elements, nofelements );
      // nofelements=interval.getHexa(elements);
      if( simplexgrid )
        nofelements = dgf :: SimplexBlock :: cube2simplex( vtx, elements, elParams );

      // remove copies of vertices
      if ( bvtx.isactive() || interval.numIntervals()>1)
        removeCopies();
    }
    else
    {
      if( bvtx.isactive() )
      {
        nofvtx = bvtx.get( vtx, vtxParams, nofvtxparams );
        info->block( bvtx );
        vtxoffset = bvtx.offset();
      }

      dgf :: SimplexBlock bsimplex( gridin, nofvtx, vtxoffset, dimgrid );
      dgf :: CubeBlock bcube( gridin, nofvtx, vtxoffset, dimgrid );
      dgf :: GeneralBlock bgeneral( gridin, nofvtx, vtxoffset, dimgrid );

      if( bcube.isactive() && (element != Simplex) )
      {
        info->block( bcube );
        nofelements = bcube.get( elements, elParams, nofelparams );
        if( bsimplex.isactive() && (element == General) )
        {
          size_t nofsimplex = bsimplex.get( elements, elParams, nofelparams );

          // if no simplices found but simplex block active, convert cubes to simplices
          if( nofsimplex == 0 )
          {
            // make simplex grid
            info->cube2simplex( element );
            nofelements = dgf :: SimplexBlock :: cube2simplex( vtx, elements, elParams );
            simplexgrid = true;
            nofelements += bsimplex.get( elements, elParams, nofelparams );
          }
        }
      }
      else if( bgeneral.isactive() )
      {
        info->block(bgeneral);
        nofelements += bgeneral.get( elements, elParams, nofelparams );
      }
      else
      {
        simplexgrid = true;
        if( bcube.isactive() )
        {
          info->block( bcube );
          info->cube2simplex( element );
          nofelements = bcube.get( elements, elParams, nofelparams );
          // make simplex grid
          nofelements = dgf :: SimplexBlock :: cube2simplex( vtx, elements, elParams );
        }
        if( bsimplex.isactive() )
        {
          info->block(bsimplex);
          nofelements += bsimplex.get( elements, elParams, nofelparams );
          if( dimgrid == 2 && dimw == 2)
          {
            for( size_t i = 0; i < elements.size(); ++i )
              testTriang( i );
          }
        }
      }
    }

    info->step1( dimw, vtx.size(), elements.size() );
    // test for tetgen/triangle block (only if simplex-grid allowed)
    if( (element != Cube) && dgf :: SimplexGenerationBlock( gridin ).isactive() )
    {
      if( !interval.isactive() )
        generateBoundaries( gridin, true );
      info->automatic();
      simplexgrid = true;
      nofelements = 0;
      generateSimplexGrid( gridin );
    }
    generateBoundaries( gridin, !interval.isactive() );

    if( nofelements<=0 )
    {
      std::cerr << "No suitable element block found in dgf file "
        << "for the chosen grid type - "
        << "perhaps you are building a grid supporting cube elements "
        << "from a dgf file containing only a simplex block?" << std::endl;
      DUNE_THROW( DGFException, "Error: No usable block to construct elements found." );
    }

    info->finish();
    delete info;
    info = 0;
    // we made it -
    // although prehaps a few boundary segments are still without id :-<
    return true;
  }


  void DuneGridFormatParser :: removeCopies ()
  {
    std::vector<int> map(vtx.size());
    std::vector<int> shift(vtx.size());
    for (size_t i=0; i<vtx.size(); i++) {
      map[i]=i;
      shift[i]=0;
    }
    nofvtx = vtx.size();
    for (size_t i=0; i<vtx.size(); i++) {
      if ((size_t)map[i]!=i) continue;
      for (size_t j=i+1; j<vtx.size(); j++) {
        double len=std::abs(vtx[i][0]-vtx[j][0]);
        for (int p=1; p<dimw; p++)
          len+=std::abs(vtx[i][p]-vtx[j][p]);
        if (len<minVertexDistance) {
          map[j]=i;
          for (size_t k=j+1; k<vtx.size(); k++)
            shift[k]++;
          nofvtx--;
        }
      }
    }
    for (size_t i=0; i<elements.size(); i++) {
      for (size_t j=0; j<elements[i].size(); j++) {
        elements[i][j]=map[elements[i][j]];
        elements[i][j]-=shift[elements[i][j]];
      }
    }
    for (size_t j=0; j<vtx.size(); j++) {
      vtx[j-shift[j]]=vtx[j];
    }
    vtx.resize(nofvtx);
    assert(vtx.size()==size_t(nofvtx));
  }


  void DuneGridFormatParser
  :: generateBoundaries ( std :: istream &gridin, bool bndseg )
  {
    // **************************************************
    // up to here:
    // filled vertex and element block, now look at boundaries...
    // **************************************************
    // first read macrogrid segments...
    if (bndseg)
    {
      dgf :: BoundarySegBlock segbound(gridin, nofvtx,dimw,simplexgrid);
      if (segbound.isactive())
      {
        info->block(segbound);
        nofbound=segbound.get(facemap,(nofelements>0),vtxoffset);

        // check whether we have boundary parameters
        facemap_t :: iterator pos = facemap.begin();
        for( ; pos != facemap.end(); ++pos )
        {
          if( !pos->second.second.empty() )
            break;
        }
        haveBndParameters = ( pos != facemap.end() );
      }
    }

    // if no boundary elements, return
    if (nofelements==0) return;

    dgf :: BoundaryDomBlock dombound(gridin, dimw);

    if ( !dombound.isactive() && facemap.empty() ) return;

    // now add all boundary faces
    {
      facemap_t :: iterator pos;
      for(int simpl=0; simpl < nofelements ; simpl++)
      {
        const int nFaces = ElementFaceUtil :: nofFaces( dimw, elements[ simpl ]);
        for( int i =0 ; i < nFaces; ++i )
        {
          //DGFEntityKey< unsigned int > key2
          facemap_t :: key_type key2
            = ElementFaceUtil :: generateFace( dimw, elements[ simpl ], i );

          pos=facemap.find(key2);
          if(pos == facemap.end())
          {
            auto& entry = facemap[key2];
            entry.first  = 0;
            entry.second = DGFBoundaryParameter::defaultValue();
          }
          else if (pos->second.first==0 ||
                   pos->first.origKeySet())
          { // face found twice
            facemap.erase(pos);
          }
          else { // use original key as given in key2
            const auto value = pos->second;
            facemap.erase(pos);
            facemap[key2] = value;
          }
        }
      }
    }
    // remove unused boundary faces added through boundaryseg block or cube2simplex conversion
    {
      facemap_t :: iterator pos = facemap.begin();
      while( pos != facemap.end() )
      {
        if( !pos->first.origKeySet() )
        {
          facemap_t :: iterator pos1 = pos;
          ++pos;
          facemap.erase( pos1 );
        }
        else
          ++pos;
      }
    }

    // now try to assign boundary ids...
    int remainingBndSegs = 0;
    int defaultBndSegs = 0;
    int inbnddomain = 0;
    if ( dombound.isactive() )
    {
      info->block( dombound );
      std::vector < std::vector < double > > v;
      for (auto& pos : facemap)
      {
        if( pos.second.first != 0 )
          continue;

        v.resize( pos.first.size() );
        for( int i = 0; i < pos.first.size(); i++ )
          v[ i ] = vtx[ pos.first[ i ] ];
        const dgf::DomainData * data = dombound.contains( v );
        if ( data )
        {
          pos.second.first = data->id();
          pos.second.second = data->parameter();
          data->isDefault() ? defaultBndSegs ++ : inbnddomain++;
        }
      }
    }

    for (const auto& pos : facemap)
    {
      if( pos.second.first == 0 )
        remainingBndSegs++;
    }

    haveBndParameters = ( haveBndParameters || dombound.hasParameter() );

    info->step2(nofbound,facemap.size(),inbnddomain,defaultBndSegs,remainingBndSegs);
  }


  /*************************************************************
     caller to tetgen/triangle
   ****************************************************/
  void DuneGridFormatParser :: generateSimplexGrid ( std :: istream &gridin )
  {
    dgf :: SimplexGenerationBlock para(gridin);
    info->block(para);

    // check whether a dump file name was provided
    std::string name = para.dumpFileName();
    const bool tempFile = name.empty();

    if( para.hasfile() ) // a triangle/tetgen file is provided
      name = para.filename();
    else if( tempFile ) // we have no dump file name and no triangle/tetgen file -> use temporary file
      name = temporaryFileName();

    std::string suffix;
    std::string params;

    if( para.hasfile() )
    {
      if (para.filetype().size()==0) {
        readTetgenTriangle(name);
        return;
      }
      dimw = para.dimension();
      if (dimw!=2 && dimw!=3 && dimw!=-1) {
        DUNE_THROW(DGFException,
                   "SimplexGen can only generate 2d or 3d meshes but not in "
                   << dimw << " dimensions!");
      } else if (dimw==-1) {
        DUNE_THROW(DGFException,
                   "SimplexGen: cannot determine dimension of grid," <<
                   " include parameter DIMENSION in the "
                                                                     << "Simplexgeneration-Block");
      }
    }
    else
      writeTetgenPoly( name, suffix, params );

    int call_nr = 1;
    if( dimw == 2 )
    {
      std::stringstream command;
      command << std::fixed;

      if (para.haspath())
        command << para.path() << "/";
      command << "triangle -ej " << params;
      if(para.hasfile())
      {
        // name = para.filename();
        suffix = "."+para.filetype();
        command << " " << para.parameter() << " ";
      }
      else
      {
        suffix = "";
      }

      if (para.minAngle()>0)
        command << "-q" << para.minAngle() << " ";
      if (para.maxArea()>0)
        command << "-a" << para.maxArea() << " ";
      command << name << suffix;
      dverb << "Calling : " << command.str() << std::endl;
      info->print("Calling : "+command.str());
      const int status = system( command.str().c_str() );
      if( status != 0 )
      {
        if( status < 0 )
          DUNE_THROW( SystemError, "Unable to call '" << command.str() << "'." );
        else
          DUNE_THROW( DGFException, "Triangle finished unsuccessfully (Command: '" << command.str() << "')." );
      }
      if (para.display())
      {
        std::stringstream paraCommand;
        if (para.haspath())
          paraCommand << para.path() << "/";
        paraCommand << "showme " << name; // << ".1.ele";
        dverb << "Calling : " << paraCommand.str() << std::endl;
        if( system( paraCommand.str().c_str() ) < 0 )
          DUNE_THROW( SystemError, "Unable to call " << paraCommand.str() << "." );
      }
    }
    else if( dimw == 3 )
    {
      { // first call
        std::stringstream command;
        command << std::fixed;
        // std::string suffix;

        if (para.haspath())
          command << para.path() << "/";
        command << "tetgen " << params;

        if( para.hasfile() )
        {
          suffix = "." + para.filetype();
          command << " " << para.parameter() << " ";
        }

        command << name << suffix;
        dverb << "Calling : " << command.str() << std::endl;
        info->print("Calling : "+command.str());

        const int status = system( command.str().c_str() );
        if( status != 0 )
        {
          if( status < 0 )
            DUNE_THROW( SystemError, "Unable to call '" << command.str() << "'." );
          else
            DUNE_THROW( DGFException, "TetGen finished unsuccessfully (Command: '" << command.str() << "')." );
        }
      }
      if (para.minAngle()>0 || para.maxArea()>0)
      { // second call
        info->print("Quality enhancement:");
        call_nr = 2;
        std::stringstream command;
        command << std::fixed;
        if (para.haspath())
          command << para.path() << "/";

        command << "tetgen -r";
        if (para.minAngle()>0)
          command << "q" << para.minAngle();

        if (para.maxArea()>0)
          command << "a" << para.maxArea();

        command << " " << name << ".1";
        dverb << "Calling : " << command.str() << std::endl;
        info->print("Calling : "+command.str());
        if( system(command.str().c_str()) < 0 )
          DUNE_THROW( SystemError, "Unable to call " << command.str() << "." );
      }
      if (para.display())
      {
        std::stringstream command;
        if (para.haspath())
          command << para.path() << "/";

        command << "tetview-linux " << name << "." << call_nr << ".ele";
        dverb << "Calling : " << command.str() << std::endl;
        if( system(command.str().c_str()) < 0 )
          DUNE_THROW( SystemError, "Unable to call " << command.str() << "." );
      }
    }
    std::stringstream polyname;
    polyname << name << "." << call_nr;
    readTetgenTriangle(polyname.str());

    if ( tempFile &&
         name.compare(0,12,"TMPDGFParser")==0 )
    {
      /*
         std::stringstream command;
         command << "rm -f "
              << prefixname << "*.ele "
              << prefixname << "*.edge "
              << prefixname << "*.node "
              << prefixname << "*.poly ";
         if( system(command.str().c_str()) < 0 )
         DUNE_THROW( SystemError, "Unable to call " << command.str() << "." );
       */
    }

    info->print("Automatic grid generation finished");
  }


  void DuneGridFormatParser :: readTetgenTriangle ( const std :: string &name )
  {
    int offset,bnd;
    std::string nodename = name + ".node";
    std::string elename = name + ".ele";
    std::string polyname = name + ((dimw==2) ? ".edge" : ".face");
    dverb << "opening " << nodename << "\n";
    std::ifstream node(nodename.c_str());
    if (!node) {
      DUNE_THROW(DGFException,
                 "could not find file " << nodename
                                        << " prehaps something went wrong with Tetgen/Triangle?");
    }
    dverb << "opening " << elename << "\n";
    std::ifstream ele(elename.c_str());
    if (!ele) {
      DUNE_THROW(DGFException,
                 "could not find file " << elename
                                        << " prehaps something went wrong with Tetgen/Triangle?");
    }
    {
      dverb << "calculating offset from " << name << " .... offset = ";
      node >> nofvtx >> dimw >> nofvtxparams >> bnd;
      // offset is 0 by default
      // the offset it the difference of the first vertex number to zero
      node >> offset;
      dverb << offset << " \n";
      node.seekg(0);
    }
    {
      int tmp;
      // first token is number of vertex which should equal i
      node >> nofvtx >> dimw >> nofvtxparams >> bnd;
      vtx.resize(nofvtx);
      if (nofvtxparams>0)
        vtxParams.resize(nofvtx);
      for (int i=0; i<nofvtx; i++) {
        vtx[i].resize(dimw);
        int nr;
        node >> nr;
        // first token is number of vertex which should equal i
        assert(nr-offset==i);
        for (int v=0; v<dimw; v++)
          node >> vtx[i][v];
        if (nofvtxparams>0) {
          vtxParams[i].resize(nofvtxparams);
          for (int p=0; p<nofvtxparams; p++)
            node >> vtxParams[i][p];
        }
        for (int p=0; p<bnd; p++)
          node >> tmp;
      }
    }
    {
      int tmp;
      ele >> nofelements >> tmp >> nofelparams;
      elements.resize(nofelements);
      if (nofelparams>0)
        elParams.resize(nofelements);
      for (int i=0; i<nofelements; i++) {
        elements[i].resize(dimw+1);
        int nr;
        ele >> nr;
        assert(nr-offset==i);
        for (int v=0; v<dimw+1; v++) {
          int elno;
          ele >> elno;
          elements[i][v] = elno - offset;
        }
        if (nofelparams>0) {
          elParams[i].resize(nofelparams);
          for (int p=0; p<nofelparams; p++)
            ele >> elParams[i][p];
        }
      }
    }
    dverb << "opening " << polyname << "\n";
    std::ifstream poly(polyname.c_str());
    if (poly) {
      /*
         if (dimw==2) {
         poly >> tmp >> tmp >> tmp >> tmp;
         }
       */
      {
        int noffaces,params;
        poly >> noffaces >> params;
        if (params>0) {
          assert(params==1);
          facemap.clear();
          for (int i=0; i<noffaces; i++) {
            std :: vector< unsigned int > p( dimw );
            int nr;
            poly >> nr;
            for (size_t k=0; k<p.size(); k++)
              poly >> p[k];
            poly >> params;
            if (params!=0)
            {
              facemap_t :: key_type key( p, false );
              // read bnd id as positive value (tetgen uses negative values)
              facemap[key].first = std::abs(params);
              facemap[key].second.clear();
            }
          }
        }
      }
    }
  }


  /***************************
     Helper methods mostly only for simplex grids
  ***************************/
  void DuneGridFormatParser
  :: setOrientation ( int use1, int use2, orientation_t orientation )
  {
    if (element == Cube) {
      std::cerr << "Reorientation is only implemented for simplex grid!"
                << std::endl;
      return;
    }
    if( dimgrid == 2 )
    {
      for (int i=0; i<nofelements; i++) {
        if (elements[i].size()!=size_t(dimw+1))
          continue;
        double o=testTriang(i);
        if (o*int(orientation)<0) { // wrong orientation
          // dverb << "Reorientation of simplex " << i << std::endl;
          int tmp=elements[i][use1];
          elements[i][use1] = elements[i][use2];
          elements[i][use2] = tmp;
        }
      }
    }
    else if (dimw==3)
    {
      auto refElem = ReferenceElements< double, 3 >::simplex();
      for (int i=0; i<nofelements; i++)
      {
        if (elements[i].size()!=size_t(dimw+1))
          continue;

        std::vector<double>& p0 = vtx[elements[i][1]];
        std::vector<double>& p1 = vtx[elements[i][2]];
        std::vector<double>& p2 = vtx[elements[i][3]];
        std::vector<double>& q  = vtx[elements[i][0]];

        double n[3];
        n[0] = -((p1[1]-p0[1]) *(p2[2]-p0[2]) - (p2[1]-p0[1]) *(p1[2]-p0[2])) ;
        n[1] = -((p1[2]-p0[2]) *(p2[0]-p0[0]) - (p2[2]-p0[2]) *(p1[0]-p0[0])) ;
        n[2] = -((p1[0]-p0[0]) *(p2[1]-p0[1]) - (p2[0]-p0[0]) *(p1[1]-p0[1])) ;
        double test = n[0]*(q[0]-p0[0])+n[1]*(q[1]-p0[1])+n[2]*(q[2]-p0[2]);
        bool reorient = (test*int(orientation)<0);
        if (reorient)
        {
          // reorient element first
          int tmp=elements[i][use1];
          elements[i][use1] = elements[i][use2];
          elements[i][use2] = tmp;

          // for all faces, check reorient
          for(int k=0; k<refElem.size(1); ++k)
          {
            int numVerts = refElem.size(k,1,dimw);
            std :: vector< unsigned int > face( numVerts );
            // get face vertices
            for(int j=0; j<numVerts; ++j)
            {
              int vx = refElem.subEntity(k,1,j,dimw);
              face[j] = elements[i][vx];
            }

            {
              // key is now right oriented
              facemap_t :: key_type key( face );
              facemap_t :: iterator bndFace = facemap.find( key );
              if( bndFace != facemap.end() )
              {
                // delete old key, and store new key
                BndParam bndId = bndFace->second;
                facemap.erase(bndFace);
                facemap[key] = bndId;
              }
            }
          }

        }
      }
    }
  }


  void DuneGridFormatParser::setRefinement ( int use1, int use2, int is1, int is2 )
  {
    if( use1 > use2 )
      std::swap( use1, use2 );
    if( is1 > is2 )
      std::swap( is1, is2 );

    if( element != Cube )
    {
      for( int i = 0; i < nofelements; ++i )
      {
        if( elements[ i ].size() != size_t( dimgrid + 1 ) )
          continue;

        double maxlen = 0.0;
        int vtx1 = is1;
        int vtx2 = is2;
        if( (vtx1 == -1) || (vtx2 == -1) )
        {
          for( int l = 0; l <= dimgrid; ++l )
          {
            const int idxl = elements[ i ][ l ];
            for( int k = l+1; k <= dimgrid; ++k )
            {
              const int idxk = elements[ i ][ k ];
              double len = 0.0;
              for( int p = 0; p < dimw; ++p )
              {
                const double dist = vtx[ idxk ][ p ] - vtx[ idxl ][ p ];
                len += dist * dist;
              }
              if( len > maxlen )
              {
                vtx1 = l;
                vtx2 = k;
                maxlen = len;
              }
            }
          }
        }

        int swapped = 0;
        if( vtx1 != use1 )
        {
          std::swap( elements[ i ][ vtx1 ], elements[ i ][ use1 ] );
          ++swapped;
        }
        if( vtx2 != use2 )
        {
          std::swap( elements[ i ][ vtx2 ], elements[ i ][ use2 ] );
          ++swapped;
        }
        if( swapped == 1 )
          std::swap( elements[ i ][ use1 ], elements[ i ][ use2 ] );
      }
    }
    else
      std::cerr << "Computing refinement edge is only implemented for 2d simplex grids." << std::endl;
  }


  double DuneGridFormatParser :: testTriang ( int snr )
  {
    double o =
      (vtx[elements[snr][1]][0]-vtx[elements[snr][0]][0])*
      (vtx[elements[snr][2]][1]-vtx[elements[snr][1]][1])-
      (vtx[elements[snr][1]][1]-vtx[elements[snr][0]][1])*
      (vtx[elements[snr][2]][0]-vtx[elements[snr][1]][0]);
    if (fabs(o)<1e-10) {
      DUNE_THROW(DGFException,
                 "Simplex number " << snr << " with vertex numbers "
                                   << "(" << elements[snr][0]
                                   << "," << elements[snr][1]
                                   << "," << elements[snr][2] << ")"
                                   << " has zero volume!");
    }
    return o;
  }


  std :: vector< double > &
  DuneGridFormatParser :: getElParam ( int i, std::vector< double > &coord )
  {
    coord.resize(dimw);
    for (int j=0; j<dimw; j++)
      coord[j]=0.;
    for (int j=0; j<dimw; j++) {
      for (size_t k=0; k<elements[i].size(); ++k) {
        coord[j]+=vtx[elements[i][k]][j];
      }
      coord[j]/=double(elements[i].size());
    }
    return elParams[i];
  }


  std :: vector< double > &
  DuneGridFormatParser :: getVtxParam ( int i, std :: vector< double > &coord )
  {
    coord.resize(dimw);
    for (int j=0; j<dimw; j++)
      coord[j]=0.;
    coord = vtx[i];
    return vtxParams[i];
  }

  inline std::string
  DuneGridFormatParser::temporaryFileName ()
  {
#if HAVE_MKSTEMP
    char filetemp[ FILENAME_MAX ];
    std :: strcpy( filetemp, "dgfparser.XXXXXX" );
    const int fd = mkstemp( filetemp );
    if( fd < 0 )
      DUNE_THROW( IOError, "Unable to create temporary file." );
    close( fd );
    return std :: string( filetemp );
#else
    char buffer[ L_tmpnam ]; // supply buffer to make it thread safe
    return std::string( std::tmpnam( buffer ) );
#endif
  }

} // end namespace Dune
