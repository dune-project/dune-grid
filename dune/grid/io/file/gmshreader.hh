// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GMSHREADER_HH
#define DUNE_GMSHREADER_HH

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <stdio.h>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune
{

  /**
     \ingroup Gmsh
     \{
   */

  //! Options for read operation
  struct GmshReaderOptions
  {
    enum GeometryOrder {
      /** @brief edges are straight lines. */
      firstOrder,
      /** @brief quadratic boundary approximation. */
      secondOrder
    };
  };

  namespace {

    // arbitrary dimension, implementation is in specialization
    template< int dimension, int dimWorld = dimension >
    class GmshReaderQuadraticBoundarySegment
    {};

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
      typedef Dune::FieldVector< double, dimWorld > GlobalVector;

      GmshReaderQuadraticBoundarySegment ( const GlobalVector &p0_, const GlobalVector &p1_, const GlobalVector &p2_)
        : p0(p0_), p1(p1_), p2(p2_)
      {
        GlobalVector d1 = p1;
        d1 -= p0;
        GlobalVector d2 = p2;
        d2 -= p1;

        alpha=d1.two_norm()/(d1.two_norm()+d2.two_norm());
        if (alpha<1E-6 || alpha>1-1E-6)
          DUNE_THROW(Dune::IOError, "ration in quadratic boundary segment bad");
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
    public:
      GmshReaderQuadraticBoundarySegment (Dune::FieldVector<double,3> p0_, Dune::FieldVector<double,3> p1_,
                                          Dune::FieldVector<double,3> p2_, Dune::FieldVector<double,3> p3_,
                                          Dune::FieldVector<double,3> p4_, Dune::FieldVector<double,3> p5_)
        : p0(p0_), p1(p1_), p2(p2_), p3(p3_), p4(p4_), p5(p5_)
      {
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

  }   // end empty namespace

  //! dimension independent parts for GmshReaderParser
  template<typename GridType>
  class GmshReaderParserBase
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

    // static data
    static const int dim = GridType::dimension;
    static const int dimWorld = GridType::dimensionworld;
    dune_static_assert( (dimWorld <= 3), "GmshReader requires dimWorld <= 3." );

    // typedefs
    typedef FieldVector< double, dimWorld > GlobalVector;

    // template<class T1, class T2, class T3, class T4, class T5, class T6,
    //          class T7, class T8, class T9>
    // void readfile(FILE * file, int cnt, const char * format,
    //     T1& t1, T2& t2, T3& t3, T4& t4, T5& t5, T6& t6, T7& t7, T8& t8, T9& t9)
    // {
    //     readfile(file, cnt, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9);
    // };

    // don't use something like
    //   readfile(file, 1, "%s\n", buf);
    // to skip the rest of of the line -- that will only skip the next
    // whitespace-seperated word!  Use skipline() instead.
    void readfile(FILE * file, int cnt, const char * format,
                  void* t1, void* t2 = 0, void* t3 = 0, void* t4 = 0,
                  void* t5 = 0, void* t6 = 0, void* t7 = 0, void* t8 = 0,
                  void* t9 = 0, void* t10 = 0)
    {
      off_t pos = ftello(file);
      int c = fscanf(file, format, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
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

    GmshReaderParserBase(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      factory(_factory), verbose(v), insert_boundary_segments(i) {}

    std::vector<int> & boundaryIdMap()
    {
      return boundary_id_to_physical_entity;
    }

    std::vector<int> & elementIndexMap()
    {
      return element_index_to_physical_entity;
    }

    void read (const std::string& f)
    {
      if (verbose) std::cout << "Reading " << dim << "d Gmsh grid..." << std::endl;

      // open file name, we use C I/O
      fileName = f;
      FILE* file = fopen(fileName.c_str(),"r");
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
      if( (version_number < 2.0) || (version_number > 2.2) )
        DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files");
      if (verbose) std::cout << "version " << version_number << " Gmsh file detected" << std::endl;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // node section
      int number_of_nodes;

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      readfile(file,1,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;

      // read nodes
      std::vector< GlobalVector > nodes( number_of_nodes+1 );       // store positions
      {
        int id;
        double x[ 3 ];
        for( int i = 1; i <= number_of_nodes; ++i )
        {
          readfile(file,4, "%d %lg %lg %lg\n", &id, &x[ 0 ], &x[ 1 ], &x[ 2 ] );
          if( id != i )
            DUNE_THROW( Dune::IOError, "Expected id " << i << "(got id " << id << "." );

          // just store node position
          for( int j = 0; j < dimWorld; ++j )
            nodes[ i ][ j ] = x[ j ];
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
      // Pass 1: Renumber needed vertices
      //=========================================

      long section_element_offset = ftell(file);
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity;
        std::vector<int> mesh_partitions;
        if ( version_number < 2.2 )
        {
          mesh_partitions.resize(1);
        }
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if ( version_number < 2.2 )
          {
            if (k==3) mesh_partitions[0] = blub;
          }
          else
          {
            if (k > 3)
              mesh_partitions[k-4] = blub;
            else
              mesh_partitions.resize(blub);
          }
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

      fseek(file, section_element_offset, SEEK_SET);
      boundary_element_count = 0;
      element_count = 0;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        int physical_entity = -1, elementary_entity = -1;
        std::vector<int> mesh_partitions;
        if ( version_number < 2.2 )
        {
          mesh_partitions.resize(1);
        }
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if ( version_number < 2.2 )
          {
            if (k==3) mesh_partitions[0] = blub;
          }
          else
          {
            if (k > 3)
              mesh_partitions[k-4] = blub;
            else
              mesh_partitions.resize(blub);
          }
        }
        pass2HandleElement(file, elm_type, renumber, nodes, physical_entity);
      }
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);
    }

    // dimension dependent routines
    void pass1HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes)
    {
      // some data about gmsh elements
      const int nDofs[12]      = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, -1, 10};
      const int nVertices[12]  = {-1, 2, 3, 4, 4, 8, 6, 5, 2, 3, -1, 4};
      const int elementDim[12] = {-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, -1, 3};

      // test whether we support the element type
      if ( not (elm_type >= 0 && elm_type < 12         // index in suitable range?
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

    virtual void pass2HandleElement(FILE* file, const int elm_type,
                                    std::map<int,unsigned int> & renumber,
                                    const std::vector< GlobalVector > & nodes,
                                    const int physical_entity)
    {
      // some data about gmsh elements
      const int nDofs[12]      = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, -1, 10};
      const int nVertices[12]  = {-1, 2, 3, 4, 4, 8, 6, 5, 2, 3, -1, 4};
      const int elementDim[12] = {-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, -1, 3};

      // test whether we support the element type
      if ( not (elm_type >= 0 && elm_type < 12         // index in suitable range?
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
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          break;
        case 2 :            // 3-node triangle
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          break;
        case 3 :            // 4-node quadrilateral
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim),vertices);
          break;
        case 4 :            // 4-node tetrahedron
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          break;
        case 5 :            // 8-node hexahedron
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim),vertices);
          break;
        case 6 :            // 6-node prism
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism,dim),vertices);
          break;
        case 7 :            // 5-node pyramid
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid,dim),vertices);
          break;
        case 9 :            // 6-node triangle
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          break;
        case 11 :            // 10-node tetrahedron
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
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

          case 8 : {              // 3-node line
            array<FieldVector<double,dim>, 3> v;
            for (int i=0; i<dim; i++) {
              v[0][i] = nodes[elementDofs[0]][i];
              v[1][i] = nodes[elementDofs[2]][i];                    // yes, the renumbering is intended!
              v[2][i] = nodes[elementDofs[1]][i];
            }
            BoundarySegment<dim,dimWorld>* newBoundarySegment
              = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 2, dimWorld >(v[0], v[1], v[2]);
            factory.insertBoundarySegment(vertices,
                                          shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment));
            break;
          }
          case 9 : {              // 6-node triangle
            array<FieldVector<double,dim>, 6> v;
            for (int i=0; i<6; i++)
              for (int j=0; j<dim; j++)
                v[i][j] = nodes[elementDofs[i]][j];

            BoundarySegment<dim,dimWorld>* newBoundarySegment
              = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 3, 3 >(v[0], v[1], v[2],
                                                                                                v[3], v[4], v[5]);

            factory.insertBoundarySegment(vertices,
                                          shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment));
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

  //! Parser for Gmsh data. The code is in various specializations of this class
  template<typename GridType, int dimension>
  class GmshReaderParser
  {};

  //! Parser for Gmsh data. Specialization for 1D
  template<typename GridType>
  class GmshReaderParser<GridType, 1> : public GmshReaderParserBase< GridType >
  {
    typedef GmshReaderParserBase< GridType > Base;
  public:
    GmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      Base(_factory,v,i) {}
  };

  //! Parser for Gmsh data. Specialization for 2D
  template<typename GridType>
  class GmshReaderParser<GridType, 2> : public GmshReaderParserBase< GridType >
  {
    typedef GmshReaderParserBase< GridType > Base;
  public:
    GmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      Base(_factory,v,i) {}
  };

  //! Parser for Gmsh data. Specialization for 3D
  template<typename GridType>
  class GmshReaderParser<GridType, 3> : public GmshReaderParserBase< GridType >
  {
    typedef GmshReaderParserBase< GridType > Base;
  public:
    GmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      Base(_factory,v,i) {}
  };

  /**
     \ingroup Gmsh

     \brief Read Gmsh mesh file

     Read a .msh file generated using Gmsh and construct a grid using the grid factory interface.

     The file format used by gmsh can hold grids that are more general than the simplex grids that
     the gmsh grid generator is able to construct.  We try to read as many grids as possible, as
     long as they are valid files.  You can test this by checking whether gmsh will load the file
     and display its content.

     All grids in a gmsh file live in three-dimensional Euclidean space.  If the world dimension
     of the grid type that you are reading the file into is less than three, the remaining coordinates
     are simply ignored.
   */
  template<typename GridType>
  class GmshReader
  {
  public:
    /** \todo doc me */
    static GridType* read (const std::string& fileName, bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<GridType> factory;

      // create parse object
      GmshReaderParser<GridType,GridType::dimension> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      return factory.createGrid();
    }

    /** \todo doc me */
    static GridType* read (const std::string& fileName,
                           std::vector<int>& boundary_id_to_physical_entity,
                           std::vector<int>& element_index_to_physical_entity,
                           bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<GridType> factory;

      // create parse object
      GmshReaderParser<GridType,GridType::dimension> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      boundary_id_to_physical_entity.swap(parser.boundaryIdMap());
      element_index_to_physical_entity.swap(parser.elementIndexMap());

      return factory.createGrid();
    }

    /** \todo doc me */
    static GridType* read (GridType& grid, const std::string& fileName,
                           bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<GridType> factory(&grid);

      // create parse object
      GmshReaderParser<GridType,GridType::dimension> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      return factory.createGrid();
    }

    /** \todo doc me */
    static GridType* read (GridType& grid, const std::string& fileName,
                           std::vector<int>& boundary_id_to_physical_entity,
                           std::vector<int>& element_index_to_physical_entity,
                           bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<GridType> factory(&grid);

      // create parse object
      GmshReaderParser<GridType,GridType::dimension> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      boundary_id_to_physical_entity.swap(parser.boundaryIdMap());
      element_index_to_physical_entity.swap(parser.elementIndexMap());

      return factory.createGrid();
    }


    /** \todo doc me */
    static void read (Dune::GridFactory<GridType>& factory, const std::string& fileName,
                      bool verbose = true, bool insert_boundary_segments=true)
    {
      // create parse object
      GmshReaderParser<GridType,GridType::dimension> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);
    }

    /** \todo doc me */
    static void read (Dune::GridFactory<GridType>& factory,
                      const std::string& fileName,
                      std::vector<int>& boundary_id_to_physical_entity,
                      std::vector<int>& element_index_to_physical_entity,
                      bool verbose = true, bool insert_boundary_segments=true)
    {
      // create parse object
      GmshReaderParser<GridType,GridType::dimension> parser(factory,verbose,insert_boundary_segments);
      parser.read(fileName);

      boundary_id_to_physical_entity.swap(parser.boundaryIdMap());
      element_index_to_physical_entity.swap(parser.elementIndexMap());
    }
  };

  /** \} */

} // namespace Dune

#endif
