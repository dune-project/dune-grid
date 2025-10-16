// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_GMSH_GMSH2PARSER_HH
#define DUNE_GRID_IO_FILE_GMSH_GMSH2PARSER_HH

#include <cstdarg>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune::Impl::Gmsh
{
  // arbitrary dimension, implementation is in specialization
  template< int dimension, int dimWorld = dimension >
  class GmshReaderQuadraticBoundarySegment
  {
  public:
    // empty function since this class does not implement anything
    static void registerFactory() {}
  };

  // quadratic boundary segments in 1d
  /*
     Note the points

     (0)   (alpha)   (1)

     are mapped to the points in global coordinates

     p0 p2 p1

     alpha is determined automatically from the given points.
   */
  template< int dimWorld >
  struct GmshReaderQuadraticBoundarySegment< 2, dimWorld >
    : public Dune::BoundarySegment< 2, dimWorld >
  {
    typedef GmshReaderQuadraticBoundarySegment< 2, dimWorld >  ThisType;
    typedef typename Dune::BoundarySegment< 2, dimWorld > :: ObjectStreamType ObjectStreamType;
    typedef Dune::FieldVector< double, dimWorld > GlobalVector;

    GmshReaderQuadraticBoundarySegment ( const GlobalVector &p0_, const GlobalVector &p1_, const GlobalVector &p2_)
      : p0(p0_), p1(p1_), p2(p2_)
    {
      init();
    }

    GmshReaderQuadraticBoundarySegment( ObjectStreamType& in )
    {
      // key is read before by the factory
      const int bytes = sizeof(double)*dimWorld;
      in.read( (char *) &p0[ 0 ], bytes );
      in.read( (char *) &p1[ 0 ], bytes );
      in.read( (char *) &p2[ 0 ], bytes );
      init();
    }

    static void registerFactory()
    {
      if( key() < 0 )
      {
        key() = Dune::BoundarySegment< 2, dimWorld >::template registerFactory< ThisType >();
      }
    }

    virtual GlobalVector operator() ( const Dune::FieldVector<double,1> &local ) const
    {
      GlobalVector y;
      y = 0.0;
      y.axpy((local[0]-alpha)*(local[0]-1.0)/alpha,p0);
      y.axpy(local[0]*(local[0]-1.0)/(alpha*(alpha-1.0)),p1);
      y.axpy(local[0]*(local[0]-alpha)/(1.0-alpha),p2);
      return y;
    }

    void backup( ObjectStreamType& out ) const
    {
      // backup key to identify object
      out.write( (const char *) &key(), sizeof( int ) );
      // backup data
      const int bytes = sizeof(double)*dimWorld;
      out.write( (const char*) &p0[ 0 ], bytes );
      out.write( (const char*) &p1[ 0 ], bytes );
      out.write( (const char*) &p2[ 0 ], bytes );
    }

  protected:
    void init()
    {
      GlobalVector d1 = p1;
      d1 -= p0;
      GlobalVector d2 = p2;
      d2 -= p1;

      alpha=d1.two_norm()/(d1.two_norm()+d2.two_norm());
      if (alpha<1E-6 || alpha>1-1E-6)
        DUNE_THROW(Dune::IOError, "ration in quadratic boundary segment bad");
    }

    static int& key() {
      static int k = -1;
      return k;
    }

  private:
    GlobalVector p0,p1,p2;
    double alpha;
  };


  // quadratic boundary segments in 2d
  /* numbering of points corresponding to gmsh:

     2

     5  4

     0  3  1

     Note: The vertices 3, 4, 5 are not necessarily at the edge midpoints but can
     be placed with parameters alpha, beta , gamma at the following positions
     in local coordinates:


     2 = (0,1)

     5 = (0,beta)   4 = (1-gamma/sqrt(2),gamma/sqrt(2))

     0 = (0,0)      3 = (alpha,0)                        1 = (1,0)

     The parameters alpha, beta, gamma are determined from the given vertices in
     global coordinates.
   */
  template<>
  class GmshReaderQuadraticBoundarySegment< 3, 3 >
    : public Dune::BoundarySegment< 3 >
  {
    typedef GmshReaderQuadraticBoundarySegment< 3, 3 > ThisType;
    typedef typename Dune::BoundarySegment< 3 > :: ObjectStreamType ObjectStreamType;
  public:
    GmshReaderQuadraticBoundarySegment (Dune::FieldVector<double,3> p0_, Dune::FieldVector<double,3> p1_,
                                        Dune::FieldVector<double,3> p2_, Dune::FieldVector<double,3> p3_,
                                        Dune::FieldVector<double,3> p4_, Dune::FieldVector<double,3> p5_)
      : p0(p0_), p1(p1_), p2(p2_), p3(p3_), p4(p4_), p5(p5_)
    {
      init();
    }

    GmshReaderQuadraticBoundarySegment( ObjectStreamType& in )
    {
      const int bytes = sizeof(double)*3;
      in.read( (char *) &p0[ 0 ], bytes );
      in.read( (char *) &p1[ 0 ], bytes );
      in.read( (char *) &p2[ 0 ], bytes );
      in.read( (char *) &p3[ 0 ], bytes );
      in.read( (char *) &p4[ 0 ], bytes );
      in.read( (char *) &p5[ 0 ], bytes );
      init();
    }

    static void registerFactory()
    {
      if( key() < 0 )
      {
        key() = Dune::BoundarySegment< 3 >::template registerFactory< ThisType >();
      }
    }

    virtual Dune::FieldVector<double,3> operator() (const Dune::FieldVector<double,2>& local) const
    {
      Dune::FieldVector<double,3> y;
      y = 0.0;
      y.axpy(phi0(local),p0);
      y.axpy(phi1(local),p1);
      y.axpy(phi2(local),p2);
      y.axpy(phi3(local),p3);
      y.axpy(phi4(local),p4);
      y.axpy(phi5(local),p5);
      return y;
    }

    void backup( ObjectStreamType& out ) const
    {
      // backup key to identify object in factory
      out.write( (const char*) &key(), sizeof( int ) );
      // backup data
      const int bytes = sizeof(double)*3;
      out.write( (const char*) &p0[ 0 ], bytes );
      out.write( (const char*) &p1[ 0 ], bytes );
      out.write( (const char*) &p2[ 0 ], bytes );
      out.write( (const char*) &p3[ 0 ], bytes );
      out.write( (const char*) &p4[ 0 ], bytes );
      out.write( (const char*) &p5[ 0 ], bytes );
    }

  protected:
    void init()
    {
      using std::sqrt;
      sqrt2 = sqrt(2.0);
      Dune::FieldVector<double,3> d1,d2;

      d1 = p3; d1 -= p0;
      d2 = p1; d2 -= p3;
      alpha=d1.two_norm()/(d1.two_norm()+d2.two_norm());
      if (alpha<1E-6 || alpha>1-1E-6)
        DUNE_THROW(Dune::IOError, "alpha in quadratic boundary segment bad");

      d1 = p5; d1 -= p0;
      d2 = p2; d2 -= p5;
      beta=d1.two_norm()/(d1.two_norm()+d2.two_norm());
      if (beta<1E-6 || beta>1-1E-6)
        DUNE_THROW(Dune::IOError, "beta in quadratic boundary segment bad");

      d1 = p4; d1 -= p1;
      d2 = p2; d2 -= p4;
      gamma=sqrt2*(d1.two_norm()/(d1.two_norm()+d2.two_norm()));
      if (gamma<1E-6 || gamma>1-1E-6)
        DUNE_THROW(Dune::IOError, "gamma in quadratic boundary segment bad");
    }

    static int& key() {
      static int k = -1;
      return k;
    }

  private:
    // The six Lagrange basis function on the reference element
    // for the points given above

    double phi0 (const Dune::FieldVector<double,2>& local) const
    {
      return (alpha*beta-beta*local[0]-alpha*local[1])*(1-local[0]-local[1])/(alpha*beta);
    }
    double phi3 (const Dune::FieldVector<double,2>& local) const
    {
      return local[0]*(1-local[0]-local[1])/(alpha*(1-alpha));
    }
    double phi1 (const Dune::FieldVector<double,2>& local) const
    {
      return local[0]*(gamma*local[0]-(sqrt2-gamma-sqrt2*alpha)*local[1]-alpha*gamma)/(gamma*(1-alpha));
    }
    double phi5 (const Dune::FieldVector<double,2>& local) const
    {
      return local[1]*(1-local[0]-local[1])/(beta*(1-beta));
    }
    double phi4 (const Dune::FieldVector<double,2>& local) const
    {
      return local[0]*local[1]/((1-gamma/sqrt2)*gamma/sqrt2);
    }
    double phi2 (const Dune::FieldVector<double,2>& local) const
    {
      return local[1]*(beta*(1-gamma/sqrt2)-local[0]*(beta-gamma/sqrt2)-local[1]*(1-gamma/sqrt2))/((1-gamma/sqrt2)*(beta-1));
    }

    Dune::FieldVector<double,3> p0,p1,p2,p3,p4,p5;
    double alpha,beta,gamma,sqrt2;
  };


  /** \brief A parser for the Version 2 ASCII Gmsh file format
   */
  template<typename GridType>
  class Gmsh2Parser
  {
  protected:
    // private data
    Dune::GridFactory<GridType>& factory;
    bool verbose;
    bool insert_boundary_segments;
    unsigned int number_of_real_vertices;
    int boundary_element_count;
    int element_count;
    // read buffer
    char buf[512];
    std::string fileName;
    // exported data
    std::vector<int> boundary_id_to_physical_entity;
    std::vector<int> element_index_to_physical_entity;
    std::vector<std::string> physical_entity_names;

    // static data
    static const int dim = GridType::dimension;
    static const int dimWorld = GridType::dimensionworld;
    static_assert( (dimWorld <= 3), "GmshReader requires dimWorld <= 3." );

    // typedefs
    typedef FieldVector< double, dimWorld > GlobalVector;

    // don't use something like
    //   readfile(file, 1, "%s\n", buf);
    // to skip the rest of of the line -- that will only skip the next
    // whitespace-separated word!  Use skipline() instead.
    void readfile(FILE * file, int cnt, const char * format, ...)
#ifdef __GNUC__
    __attribute__((format(scanf, 4, 5)))
#endif
    {
      std::va_list ap;
      va_start(ap, format);
      auto pos = std::ftell(file);
      int c = std::vfscanf(file, format, ap);
      va_end(ap);
      if (c != cnt)
        DUNE_THROW(Dune::IOError, "Error parsing " << fileName << " "
                   "file pos " << pos
                                                   << ": Expected '" << format << "', only read " << c << " entries instead of " << cnt << ".");
    }

    // skip over the rest of the line, including the terminating newline
    void skipline(FILE * file)
    {
      int c;
      do {
        c = std::fgetc(file);
      } while(c != '\n' && c != EOF);
    }

  public:

    Gmsh2Parser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      factory(_factory), verbose(v), insert_boundary_segments(i) {}

    std::vector<int> & boundaryIdMap()
    {
      return boundary_id_to_physical_entity;
    }

    //! Returns a map for the gmsh physical entity id (1-based index) for each entity of codim 0
    std::vector<int> & elementIndexMap()
    {
      return element_index_to_physical_entity;
    }

    //! Returns the names of the gmsh physical entities (0-based index)
    std::vector<std::string> & physicalEntityNames()
    {
      return physical_entity_names;
    }

    void read (const std::string& f)
    {
      if (verbose) std::cout << "Reading " << dim << "d Gmsh grid..." << std::endl;

      // open file name, we use C I/O
      fileName = f;
      FILE* file = fopen(fileName.c_str(),"rb");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

      //=========================================
      // Header: Read vertices into vector
      //         Check vertices that are needed
      //=========================================

      number_of_real_vertices = 0;
      boundary_element_count = 0;
      element_count = 0;

      // process header
      double version_number;
      int file_type, data_size;

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      readfile(file,3,"%lg %d %d\n",&version_number,&file_type,&data_size);
      // 2.2 is not representable as float and leads to problems on i386
      // Hence we use >= 2.00001.
      if( (version_number < 2.0) || (version_number >= 2.20001) ) // 2.2 is not representable as float and leads to problems on i386
        DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files");
      if (verbose) std::cout << "version " << version_number << " Gmsh file detected" << std::endl;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // physical name section
      physical_entity_names.clear();
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$PhysicalNames")==0) {
        int number_of_names;
        readfile(file,1,"%d\n",&number_of_names);
        if (verbose) std::cout << "file contains " << number_of_names << " physical entities" << std::endl;
        physical_entity_names.resize(number_of_names);
        std::string buf_name;
        for( int i = 0; i < number_of_names; ++i ) {
          int edim, id;
          readfile(file,3, "%d %d %s\n", &edim, &id, buf);
          buf_name.assign(buf);
          auto begin = buf_name.find_first_of('\"') + 1;
          auto end = buf_name.find_last_of('\"') - begin;
          physical_entity_names[id-1].assign(buf_name.substr(begin, end));
        }
        readfile(file,1,"%s\n",buf);
        if (strcmp(buf,"$EndPhysicalNames")!=0)
          DUNE_THROW(Dune::IOError, "expected $EndPhysicalNames");
        readfile(file,1,"%s\n",buf);
      }

      // node section
      int number_of_nodes;

      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      readfile(file,1,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;

      // read nodes
      // The '+1' is due to the fact that gmsh numbers node starting from 1 rather than from 0
      std::vector< GlobalVector > nodes( number_of_nodes+1 );
      {
        int id;
        double x[ 3 ];
        for( int i = 1; i <= number_of_nodes; ++i )
        {
          readfile(file,4, "%d %lg %lg %lg\n", &id, &x[ 0 ], &x[ 1 ], &x[ 2 ] );

          if (id > number_of_nodes) {
            DUNE_THROW(Dune::IOError,
                       "Only dense sequences of node indices are currently supported (node index "
                       << id << " is invalid).");
          }

          // just store node position
          for( int j = 0; j < dimWorld; ++j )
            nodes[ id ][ j ] = x[ j ];
        }
        readfile(file,1,"%s\n",buf);
        if (strcmp(buf,"$EndNodes")!=0)
          DUNE_THROW(Dune::IOError, "expected $EndNodes");
      }

      // element section
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      readfile(file,1,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl;

      //=========================================
      // Pass 1: Select and insert those vertices in the file that
      //    actually occur as corners of an element.
      //=========================================

      auto section_element_offset = std::ftell(file);
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          // k == 1: physical entity
          // k == 2: elementary entity (not used here)
          // if version_number < 2.2:
          //   k == 3: mesh partition 0
          // else
          //   k == 3: number of mesh partitions
          //   k => 4: mesh partition k-4
        }
        pass1HandleElement(file, elm_type, renumber, nodes);
      }
      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      if (verbose) std::cout << "number of boundary elements = " << boundary_element_count << std::endl;
      if (verbose) std::cout << "number of elements = " << element_count << std::endl;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");
      boundary_id_to_physical_entity.resize(boundary_element_count);
      element_index_to_physical_entity.resize(element_count);

      //==============================================
      // Pass 2: Insert boundary segments and elements
      //==============================================

      std::fseek(file, section_element_offset, SEEK_SET);
      boundary_element_count = 0;
      element_count = 0;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        int physical_entity = -1;

        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          if (k==1) physical_entity = blub;
        }
        pass2HandleElement(file, elm_type, renumber, nodes, physical_entity);
      }
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);
    }

    /** \brief Process one element during the first pass through the list of all elements
     *
     * Mainly, the method inserts all vertices needed by the current element,
     * unless they have been inserted already for a previous element.
     */
    void pass1HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes)
    {
      // some data about gmsh elements
      const int nDofs[16]      = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, -1, 10, -1, -1, -1, 1};
      const int nVertices[16]  = {-1, 2, 3, 4, 4, 8, 6, 5, 2, 3, -1, 4, -1, -1, -1, 1};
      const int elementDim[16] = {-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, -1, 3, -1, -1, -1, 0};

      // test whether we support the element type
      if ( not (elm_type > 0 && elm_type <= 15         // index in suitable range?
                && (elementDim[elm_type] == dim || elementDim[elm_type] == (dim-1) ) ) )         // real element or boundary element?
      {
        skipline(file);         // skip rest of line if element is unknown
        return;
      }

      // The format string for parsing is n times '%d' in a row
      std::string formatString = "%d";
      for (int i=1; i<nDofs[elm_type]; i++)
        formatString += " %d";
      formatString += "\n";

      // '10' is the largest number of dofs we may encounter in a .msh file
      std::vector<int> elementDofs(10);

      readfile(file,nDofs[elm_type], formatString.c_str(),
               &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]),
               &(elementDofs[3]),&(elementDofs[4]),&(elementDofs[5]),
               &(elementDofs[6]),&(elementDofs[7]),&(elementDofs[8]),
               &(elementDofs[9]));

      // insert each vertex if it hasn't been inserted already
      for (int i=0; i<nVertices[elm_type]; i++)
        if (renumber.find(elementDofs[i])==renumber.end())
        {
          renumber[elementDofs[i]] = number_of_real_vertices++;
          factory.insertVertex(nodes[elementDofs[i]]);
        }

      // count elements and boundary elements
      if (elementDim[elm_type] == dim)
        element_count++;
      else
        boundary_element_count++;

    }



    // generic-case: This is not supposed to be used at runtime.
    template <class E, class V, class V2>
    void boundarysegment_insert(
      const V&,
      const E&,
      const V2&
      )
    {
      DUNE_THROW(Dune::IOError, "tried to create a 3D boundary segment in a non-3D Grid");
    }

    // 3d-case:
    template <class E, class V>
    void boundarysegment_insert(
      const std::vector<FieldVector<double, 3> >& nodes,
      const E& elementDofs,
      const V& vertices
      )
    {
      std::array<FieldVector<double,dimWorld>, 6> v;
      for (int i=0; i<6; i++)
        for (int j=0; j<dimWorld; j++)
          v[i][j] = nodes[elementDofs[i]][j];

      BoundarySegment<dim,dimWorld>* newBoundarySegment
        = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 3, 3 >( v[0], v[1], v[2],
                                                                                           v[3], v[4], v[5] );

      factory.insertBoundarySegment( vertices,
                                     std::shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment) );
    }



    /** \brief Process one element during the second pass through the list of all elements
     *
     * This method actually inserts the element into the grid factory.
     */
    virtual void pass2HandleElement(FILE* file, const int elm_type,
                                    std::map<int,unsigned int> & renumber,
                                    const std::vector< GlobalVector > & nodes,
                                    const int physical_entity)
    {
      // some data about gmsh elements
      const int nDofs[16]      = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, -1, 10, -1, -1, -1, 1};
      const int nVertices[16]  = {-1, 2, 3, 4, 4, 8, 6, 5, 2, 3, -1, 4, -1, -1, -1, 1};
      const int elementDim[16] = {-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, -1, 3, -1, -1, -1, 0};

      // test whether we support the element type
      if ( not (elm_type > 0 && elm_type <= 15         // index in suitable range?
                && (elementDim[elm_type] == dim || elementDim[elm_type] == (dim-1) ) ) )         // real element or boundary element?
      {
        skipline(file);         // skip rest of line if element is unknown
        return;
      }

      // The format string for parsing is n times '%d' in a row
      std::string formatString = "%d";
      for (int i=1; i<nDofs[elm_type]; i++)
        formatString += " %d";
      formatString += "\n";

      // '10' is the largest number of dofs we may encounter in a .msh file
      std::vector<int> elementDofs(10);

      readfile(file,nDofs[elm_type], formatString.c_str(),
               &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]),
               &(elementDofs[3]),&(elementDofs[4]),&(elementDofs[5]),
               &(elementDofs[6]),&(elementDofs[7]),&(elementDofs[8]),
               &(elementDofs[9]));

      // correct differences between gmsh and Dune in the local vertex numbering
      switch (elm_type)
      {
      case 3 :          // 4-node quadrilateral
        std::swap(elementDofs[2],elementDofs[3]);
        break;
      case 5 :          // 8-node hexahedron
        std::swap(elementDofs[2],elementDofs[3]);
        std::swap(elementDofs[6],elementDofs[7]);
        break;
      case 7 :          // 5-node pyramid
        std::swap(elementDofs[2],elementDofs[3]);
        break;
      }

      // renumber corners to account for the explicitly given vertex
      // numbering in the file
      std::vector<unsigned int> vertices(nVertices[elm_type]);

      for (int i=0; i<nVertices[elm_type]; i++)
        vertices[i] = renumber[elementDofs[i]];

      // If it is an element, insert it as such
      if (elementDim[elm_type] == dim) {

        switch (elm_type)
        {
        case 1 :            // 2-node line
          factory.insertElement(Dune::GeometryTypes::line,vertices);
          break;
        case 2 :            // 3-node triangle
          factory.insertElement(Dune::GeometryTypes::triangle,vertices);
          break;
        case 3 :            // 4-node quadrilateral
          factory.insertElement(Dune::GeometryTypes::quadrilateral,vertices);
          break;
        case 4 :            // 4-node tetrahedron
          factory.insertElement(Dune::GeometryTypes::tetrahedron,vertices);
          break;
        case 5 :            // 8-node hexahedron
          factory.insertElement(Dune::GeometryTypes::hexahedron,vertices);
          break;
        case 6 :            // 6-node prism
          factory.insertElement(Dune::GeometryTypes::prism,vertices);
          break;
        case 7 :            // 5-node pyramid
          factory.insertElement(Dune::GeometryTypes::pyramid,vertices);
          break;
        case 9 :            // 6-node triangle
          factory.insertElement(Dune::GeometryTypes::triangle,vertices);
          break;
        case 11 :            // 10-node tetrahedron
          factory.insertElement(Dune::GeometryTypes::tetrahedron,vertices);
          break;
        }

      } else {
        // it must be a boundary segment then
        if (insert_boundary_segments) {

          switch (elm_type)
          {
          case 1 :              // 2-node line
            factory.insertBoundarySegment(vertices);
            break;

          case 2 :              // 3-node triangle
            factory.insertBoundarySegment(vertices);
            break;

          case 3 :            // 4-node quadrilateral
            factory.insertBoundarySegment(vertices);
            break;

          case 15 :             // 1-node point
            factory.insertBoundarySegment(vertices);
            break;

          case 8 : {              // 3-node line
            std::array<FieldVector<double,dimWorld>, 3> v;
            for (int i=0; i<dimWorld; i++) {
              v[0][i] = nodes[elementDofs[0]][i];
              v[1][i] = nodes[elementDofs[2]][i];                    // yes, the renumbering is intended!
              v[2][i] = nodes[elementDofs[1]][i];
            }
            BoundarySegment<dim,dimWorld>* newBoundarySegment
              = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 2, dimWorld >(v[0], v[1], v[2]);
            factory.insertBoundarySegment(vertices,
                                          std::shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment));
            break;
          }
          case 9 : {              // 6-node triangle
            boundarysegment_insert(nodes, elementDofs, vertices);
            break;
          }
          default: {
            DUNE_THROW(Dune::IOError, "GmshReader does not support using element-type " << elm_type << " for boundary segments");
            break;
          }

          }

        }
      }

      // count elements and boundary elements
      if (elementDim[elm_type] == dim) {
        element_index_to_physical_entity[element_count] = physical_entity;
        element_count++;
      } else {
        boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
        boundary_element_count++;
      }

    }

  };

} // namespace Dune::Impl::Gmsh

#endif
