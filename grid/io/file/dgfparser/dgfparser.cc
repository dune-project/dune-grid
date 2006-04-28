// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune {
  // Output to Alberta macrogridfile (2d/3d)
  inline void DuneGridFormatParser::writeAlberta(std::ostream& out) {
    // writes an output file for in gird type Alberta
    out << "DIM: " << dimw
        << "\n" << "DIM_OF_WORLD: "
        << dimw << "\n"
        << "\nnumber of vertices: "
        << nofvtx <<  "\nnumber of elements: "
        << nofelements << std::endl;
    out <<  "\nvertex coordinates: " <<  std::endl;
    for (int n=0; n<nofvtx; n++) {
      for (int j=0; j<dimw; j++) {
        out << vtx[n][j] << " ";
      }
      out << std::endl;
    }
    out << "\nelement vertices: "  << std::endl;
    for (int n=0; n<nofelements; n++) {
      for (int j=0; j<dimw+1; j++) {
        // riesiger Hack! die make6 methode erzeugt nur jeden zweiten Tetraeder
        // richtig, bei den geraden Tetras muessen die letzten beide Knoten
        // vertauscht werden...
        if (isInterval &&
            dimw==3 && n%2==0 && j>=2)
          if (j==2) out << " " << elements[n][3] << " ";
          else out << " " << elements[n][2] << " ";
        else
          out << " " << elements[n][j] << " ";
      }
      out << std::endl;
    }
    out << "\nelement boundaries: "  << std::endl;
    std::map<EntityKey<int>,int>::iterator pos;
    for( int simpl=0; simpl < nofelements ; simpl++) {
      for (int i =0 ; i< dimw+1 ; i++) {
        EntityKey<int> key2(elements[simpl],dimw,i+1);
        pos=facemap.find(key2);
        if (pos==facemap.end())
          out << "0 ";
        else {
          if (pos->second == 0)
            out << "E ";
          else
            out << pos->second << " ";
        }
      }
      out << " " << std::endl;
    }
  }
  // Output to ALU macrogridfile (3d tetra/hexa)
  inline void DuneGridFormatParser::writeAlu(std::ostream& out) {
    // wirtes an output file in grid type ALU
    if (dimw==3) {
      if (simplexgrid)
        out << "!Tetraeder" << std::endl;
      else
        out << "!Hexaeder" << std::endl;
    }
    if (dimw==2) {
      if (!simplexgrid) {
        std::cerr << "ERROR: ALU can only handle simplex grids in 2d!" << std::endl;
        return;
      }
    }
    std::cout << "Writing vertices...";
    out << nofvtx << std::endl;
    for (int n=0; n<nofvtx; n++) {
      for (int j=0; j<dimw; j++) {
        out << vtx[n][j] << " ";
      }
      out << std::endl;
    }
    std::cout << "done" << std::endl;
    std::cout.flush();
    std::cout << "Writing Simplices...";
    out << nofelements << std::endl;
    for (int n=0; n<nofelements; n++) {
      if (simplexgrid) {
        for (size_t j=0; j<elements[n].size(); j++) {
          out << " " << elements[n][j] << " ";
        }
      } else {
        std::vector<unsigned int> el;
        for (int j=0; j<elements[n].size(); j++)
          el.push_back((elements[n][j]));
        if (el.size()==8) {
          unsigned int tmp = el[2];
          el[2] = el[3];
          el[3] = tmp;
          tmp = el[6];
          el[6] = el[7];
          el[7] = tmp;
        }
        for (size_t j=0; j<el.size(); j++) {
          out << " " << el[j] << " ";
        }
      }
      out << std::endl;
    }
    std::cout << "done" << std::endl;
    std::cout.flush();
    std::cout << "Writing Boundary...";
    out << facemap.size() << std::endl;
    std::map<EntityKey<int>,int>::iterator pos;
    for(pos= facemap.begin(); pos!=facemap.end(); ++pos) {
      if (pos->second == 0)
        out << "E ";
      else
        out << -pos->second << " ";
      if (dimw == 3)
        out << pos->first.size() << " ";
      for (int i=0; i<pos->first.size(); i++)
        out << pos->first.origKey(i) << " ";
      out << std::endl;
    }
    if (dimw == 3)
      for (int n=0; n<nofvtx; n++) {
        out << n << " " << -1 << std::endl;
      }
    std::cout << "done" << std::endl;
    std::cout.flush();
  }
  // read the DGF file and store vertex/element/bound structure
  inline int DuneGridFormatParser::readDuneGrid(std::istream& gridin) {
    std::string id;
    getline(gridin,id);
    makeupcase(id);
    if (id != "DGF")
      return -1;     // not a DGF file, prehaps native file format
    dimw=-1;
    IntervalBlock interval(gridin);

    // first test for tetgen/triangle block (only if simplex-grid allowed)
    if (element!=Cube && SimplexGenerationBlock(gridin).isactive()) {
      simplexgrid=true;
      nofelements=0;
      int ok=generateSimplexGrid(gridin);
      if (!ok) {
        return 0;
      }
    }
    else if(interval.ok()) { // generate cartesian grid?
      isInterval = true;
      dimw = interval.dimw();
      simplexgrid = (element == Simplex);
      if (element == General) {
        SimplexBlock bsimplex(gridin,-1,dimw);
        simplexgrid = bsimplex.isactive();
      } else
        simplexgrid = (element == Simplex);
      nofvtx=interval.getVtx(vtx);
      if(simplexgrid)
        nofelements=interval.getSimplex(elements);
      else
        nofelements=interval.getHexa(elements);
    }
    else { // ok: grid generation by hand...
      VertexBlock bvtx(gridin,dimw);
      if (bvtx.isactive()) {
        nofvtx=bvtx.get(vtx);
      }
      else {
        return 0;
      }
      nofelements=0;
      SimplexBlock bsimplex(gridin,bvtx.nofvertex(),dimw);
      CubeBlock bcube(gridin,bvtx.nofvertex(),dimw);
      if (!bcube.isactive() && !bsimplex.isactive() ||
          (element==Cube && !bcube.isactive()) ) {
        std::cerr << "no element info found..." << std::endl;
        return 0;
      }
      if (bcube.isactive() && element!=Simplex) {
        nofelements=bcube.get(elements);
        if (bsimplex.isactive() || element==General) {
          // make simplex grid
          nofelements=bsimplex.cube2simplex(vtx,elements);
          simplexgrid=true;
        }
        else
          simplexgrid=false;
      }
      else {
        if (bsimplex.isactive()) {
          simplexgrid=true;
          nofelements+=bsimplex.get(elements);
          for (int i=0; i<elements.size(); i++) {
            if (testTriang(i)==0) {
              std::cerr << "Found an error in description of Simplex no. "
                        << i << std::endl;
              return 0;
            }
          }
        } else {
          nofelements=bcube.get(elements);
          // make simplex grid
          nofelements=bsimplex.cube2simplex(vtx,elements);
          simplexgrid=true;
          isInterval = true; // needed by AlbertaGrid to write correct simplex info
        }
      }
      if (nofelements<=0) {
        std::cerr << "An Error occured while reading element information ";
        std::cerr << "from the DGF file!" << std::endl;
        return 0;
      }
      // now read macrogrid segments... (works at the moment only for simplex)
      BoundarySegBlock segbound(gridin, nofvtx,dimw);
      if (segbound.isactive()) {
        std::vector<int> bound(dimw+1);
        for (nofbound=0; segbound.ok(); segbound.next(), nofbound++) {
          for (int j=0; j<dimw+1; j++) {
            bound[j] = segbound[j];
          }
          EntityKey<int> key(bound,dimw,1);
          facemap[key] = segbound[0];
        }
      }
    }
    // **************************************************
    // up to here:
    // filled vertex and element block, now look at boundaries...
    // **************************************************
    // first generate a map with all macroboundary segments...
    std::map<EntityKey<int>,int>::iterator pos;
    if(isInterval && !simplexgrid)
      interval.getCubeBoundary(facemap);
    else {
      // at the moment: automatic generation only for simplex grids
      for(int simpl=0; simpl < nofelements ; simpl++) {
        for (int i =0 ; i< dimw+1 ; i++) {
          EntityKey<int> key2(elements[simpl],dimw,i+1);
          key2.orientation(elements[simpl][i],vtx);
          pos=facemap.find(key2);
          if(pos== facemap.end()) {
            facemap[key2]=0;
          }
          else if(pos->second==0 ) {
            facemap.erase(pos);
          }
          else {
            int value = pos->second;
            facemap.erase(pos);
            facemap[key2]=value;
          }
        }
      }
    }
    // now try to assign boundary ids...
    BoundaryDomBlock dombound(gridin, dimw);
    if (dombound.isactive()) {
      for (; dombound.ok(); dombound.next()) {
        for(pos=facemap.begin(); pos!=facemap.end(); ++pos) {
          if(pos->second == 0) {
            // if an edge of a simplex is inside the domain it has the value zero
            bool isinside=true;
            for (int i=0; i<pos->first.size(); i++) {
              if (!dombound.inside(vtx[pos->first[i]])) {
                isinside=false;
                break;
              }
            }
            if (isinside)
              pos->second = dombound.id();
          }
        }
      }
      // now assign default value to remaining segments - if one was given:
      if (dombound.defaultValueGiven()) {
        for(pos=facemap.begin(); pos!=facemap.end(); ++pos) {
          if(pos->second == 0) {
            pos->second = dombound.defaultValue();
          }
        }
      }
    }
    // we made it -
    // although prehaps a few boundary segments are still without id :-<
    return 1;
  }
  /*************************************************************
     caller to tetgen/triangle
   ****************************************************/
  inline int DuneGridFormatParser::generateSimplexGrid(std::istream& gridin) {
    VertexBlock bvtx(gridin,dimw);
    IntervalBlock interval(gridin);
    if (!interval.isactive() && !bvtx.isactive()) {
      std::cerr << "No vertex information found!" << std::endl;
      return 0;
    }
    nofvtx = 0;
    if (interval.ok()) {
      std::cout << "Reading verticies from IntervalBlock" << std::flush;
      dimw = interval.dimw();
      nofvtx += interval.getVtx(vtx);
      std::cout << "Done." << std::endl;
    }
    if (bvtx.ok()) {
      std::cout << "Reading verticies from VertexBlock" << std::flush;
      nofvtx += bvtx.get(vtx);
      std::cout << "Done." << std::endl;
    }
    if (dimw!=2 && dimw!=3) {
      std::cerr << "SimplexGen can only generate 2d or 3d meshes but not in "
                << dimw << " dimensions!" << std::endl;
      return 0;
    }
    std::string name = "/tmp/gridparsertmpfile.nodelists";
    {
      std::string tmpname = name;
      tmpname += ".node";
      std::ofstream nodes(tmpname.c_str());
      nodes << nofvtx << " " << dimw << " 0 1" << std::endl;
      for (int n=0; n<nofvtx; n++) {
        nodes << n << " ";
        for (int j=0; j<dimw; j++) {
          nodes << vtx[n][j] << " ";
        }
        nodes << "1";
        nodes << std::endl;
      }
    }
    SimplexGenerationBlock para(gridin);
    int call_nr = 1;
    if (dimw==2) {
      std::stringstream command;
      if (para.haspath())
        command << para.path() << "/Triangle/";
      command << "/triangle ";
      if (para.minAngle()>0)
        command << "-q" << para.minAngle() << " ";
      if (para.maxArea()>0)
        command << "-a" << para.maxArea() << " ";
      command << name << ".node";
      std::cout << "Calling : " << command.str() << std::endl;
      system(command.str().c_str());
      if (para.display()) {
        std::stringstream command;
        if (para.haspath())
          command << para.path() << "/Triangle/";
        command << "showme " << name << ".1.ele";
        std::cout << "Calling : " << command.str() << std::endl;
        system(command.str().c_str());
      }
    }
    else if (dimw==3) {
      { // first call
        std::stringstream command;
        if (para.haspath())
          command << para.path() << "/TetGen/";

        command << "tetgen ";
        command << name << ".node";
        std::cout << "Calling : " << command.str() << std::endl;
        system(command.str().c_str());
      }
      if (para.minAngle()>0 || para.maxArea()>0)
      { // second call
        call_nr = 2;
        std::stringstream command;
        if (para.haspath())
          command << para.path() << "/TetGen/";
        command << "tetgen -r";
        if (para.minAngle()>0)
          command << "q" << para.minAngle();
        if (para.maxArea()>0)
          command << "a" << para.maxArea();
        command << " " << name << ".1.node";
        std::cout << "Calling : " << command.str() << std::endl;
        system(command.str().c_str());
      }
      if (para.display()) {
        std::stringstream command;
        if (para.haspath())
          command << para.path() << "/TetGen/";
        command << "tetview-linux " << name << "." << call_nr << ".ele";
        std::cout << "Calling : " << command.str() << std::endl;
        system(command.str().c_str());
      }
    }
    {
      std::stringstream nodename;
      nodename << name << "." << call_nr << ".node";
      int tmp,params;
      std::ifstream node(nodename.str().c_str());
      node >> nofvtx >> tmp >> tmp >> params;
      vtx.resize(nofvtx);
      for (int i=0; i<nofvtx; i++) {
        vtx[i].resize(dimw);
        int nr;
        node >> nr;
        for (int v=0; v<dimw; v++)
          node >> vtx[i][v];
        for (int p=0; p<params; p++)
          node >> tmp;
        assert(nr==i);
      }
    }
    {
      std::stringstream elname;
      elname << name << "." << call_nr << ".ele";
      int tmp;
      std::ifstream ele(elname.str().c_str());
      ele >> nofelements >> tmp >> tmp;
      elements.resize(nofelements);
      for (int i=0; i<nofelements; i++) {
        elements[i].resize(dimw+1);
        int nr;
        ele >> nr;
        for (int v=0; v<dimw+1; v++)
          ele >> elements[i][v];
        assert(nr==i);
      }
    }
    return 1;
  }
  /***************************
     Helper methods mostly only for simplex grids
  ***************************/
  inline void
  DuneGridFormatParser::setOrientation(int fixvtx,
                                       orientation_t orientation) {
    if (element == Cube) {
      std::cerr << "Reorientation is only implemented for 2d simplex grid!"
                << std::endl;
      return;
    }
    if (dimw==2) {
      for (int i=0; i<nofelements; i++) {
        double o=testTriang(i);
        if (o*orientation<0) { // wrong orientation
          // std::cout << "Reorientation of simplex " << i << std::endl;
          int tmp=elements[i][(fixvtx+1)%3];
          elements[i][(fixvtx+1)%3] = elements[i][(fixvtx+2)%3];
          elements[i][(fixvtx+2)%3] = tmp;
        }
      }
    }
    else if (dimw==3) {
      for (int i=0; i<nofelements; i++) {
        std::vector<double>& p0 = vtx[elements[i][(fixvtx+1)%4]];
        std::vector<double>& p1 = vtx[elements[i][(fixvtx+2)%4]];
        std::vector<double>& p2 = vtx[elements[i][(fixvtx+3)%4]];
        std::vector<double>& q  = vtx[elements[i][fixvtx]];
        double n[3];
        n[0] = -((p1[1]-p0[1]) *(p2[2]-p1[2]) - (p2[1]-p1[1]) *(p1[2]-p0[2])) ;
        n[1] = -((p1[2]-p0[2]) *(p2[0]-p1[0]) - (p2[2]-p1[2]) *(p1[0]-p0[0])) ;
        n[2] = -((p1[0]-p0[0]) *(p2[1]-p1[1]) - (p2[0]-p1[0]) *(p1[1]-p0[1])) ;
        double test = n[0]*(q[0]-p0[0])+n[1]*(q[1]-p0[1])+n[2]*(q[2]-p0[2]);
        bool reorient = (test*orientation<0);
        if (reorient) {
          int key1=elements[i][(fixvtx)%4];
          elements[i][(fixvtx)%4]=elements[i][(fixvtx+3)%4];
          elements[i][(fixvtx+3)%4]=key1;
        }
      }
    }
  }
  inline void
  DuneGridFormatParser::setRefinement(int fce) {
    if (dimw!=2 || element == Cube) {
      std::cerr << "Computing refinement vertex is only implemented for 2d simplex grid!"
                << std::endl;
      return;
    }
    for (int i=0; i<nofelements; i++) {
      double maxlen=0;
      int maxface;
      int idx[3] = {elements[i][0],elements[i][1],elements[i][2]};
      for (int k=0; k<dimw+1; k++) {
        double len=sqrt
                      (pow(vtx[idx[(k+2)%3]][0]-vtx[idx[(k+1)%3]][0],2.)+
                      pow(vtx[idx[(k+2)%3]][1]-vtx[idx[(k+1)%3]][1],2.)
                      );
        if (len>maxlen) {
          maxface=k;
          maxlen=len;
        }
      }
      if (maxface!=fce) {
        // std::cout << "Rearranging verticies of simplex " << i
        //	<< " since face " << maxface << " is largest" << std::endl;
        elements[i][fce]=idx[maxface];
        elements[i][(fce+1)%3]=idx[(maxface+1)%3];
        elements[i][(fce+2)%3]=idx[(maxface+2)%3];
      }
    }

  }
  inline double
  DuneGridFormatParser:: testTriang(int snr) {
    double o =
      (vtx[elements[snr][1]][0]-vtx[elements[snr][0]][0])*
      (vtx[elements[snr][2]][1]-vtx[elements[snr][1]][1])-
      (vtx[elements[snr][1]][1]-vtx[elements[snr][0]][1])*
      (vtx[elements[snr][2]][0]-vtx[elements[snr][1]][0]);
    if (fabs(o)<1e-10) {
      std::cerr << "Simplex number " << snr << " with vertex numbers "
                << "(" << elements[snr][0]
                << "," << elements[snr][1]
                << "," << elements[snr][2] << ")"
                << " has zero volume!" << std::endl;
      return 0;
    }
    return o;
  }
}
