// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

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

  }   // end empty namespace

  //! dimension independent parts for GmshReaderParser
  template<typename GridType>
  class GmshReaderParser
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
      off_t pos = ftello(file);
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

    GmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
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

      off_t section_element_offset = ftello(file);
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          // k == 1: physical entity (not used here)
          // k == 2: elementary entity (not used here either)
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

      fseeko(file, section_element_offset, SEEK_SET);
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

  namespace Gmsh {
    /**
      \ingroup Gmsh
      \brief Option for the Gmsh mesh file reader
      \note Two or more options can be composed using the binary "|"-operator
    **/
    enum class ReaderOptions
    {
      verbose = 1,
      insertBoundarySegments = 2,
      readElementData = 4,
      readBoundaryData = 8
    };

    //! composition operator for reader options
    constexpr ReaderOptions operator | (ReaderOptions a, ReaderOptions b)
    {
      return static_cast<ReaderOptions>(
        static_cast<int>(a) | static_cast<int>(b)
      );
    }

    //! query operator for reader options (is b set in a)
    constexpr bool operator & (ReaderOptions a, ReaderOptions b)
    {
      return static_cast<int>(a) & static_cast<int>(b);
    }

  } // end namespace Gmsh

  /**
     \ingroup Gmsh

     \brief Read Gmsh mesh file

     Read a .msh (version 2) file generated using Gmsh and construct a grid using the grid factory interface.

     The file format used by gmsh can hold grids that are more general than the simplex grids that
     the gmsh grid generator is able to construct.  We try to read as many grids as possible, as
     long as they are valid files.  You can test this by checking whether gmsh will load the file
     and display its content.

     All grids in a gmsh file live in three-dimensional Euclidean space.  If the world dimension
     of the grid type that you are reading the file into is less than three, the remaining coordinates
     are simply ignored.

     \note Recent versions of Gmsh introduced a new .msh file format (version 4) with a different syntax.
     This is currently not supported by GmshReader. One can export to an older .msh version as follows:
       - select File&rarr;Export (or CTRL+E)
       - select file format `.msh`
       - a dialog asks for options
       - select 'Version 2 ASCII' and mark 'Save all elements'

   */
  template<typename GridType>
  class GmshReader
  {
    //! internal general reading method
    /**
     * This method does all the highlevel steering of the reader:
     * - it will register the GmshReader boundary segment implementation with
     *   the factory
     * - it will ensure the reader is called on all ranks (debug mode only)
     * - proceed to construct a parser (rank 0 only)
     * - use the parser to read the grid into the factory (rank 0 only)
     * - move entity and boundary data from the parser into the data vector
     *   arguments, or clear the data vector arguments, depending on rank
     *
     * \note That the parser always reads the data vectors from the files.
     *       However, if insertBoundarySegments is false, no boundary segments
     *       are inserted into the factory, and thus there will be no correct
     *       indexing of the boundarySegmentToPhysicalEntity vector possible.
     *       For this reason, this method is not exposed to the user, and the
     *       interface methods are responsible to ensure that
     *       boundarySegmentToPhysicalEntity is discarded if boundary segments
     *       are not inserted.
     */
    static void doRead(Dune::GridFactory<GridType> &factory,
                       const std::string &fileName,
                       std::vector<int>& boundarySegmentToPhysicalEntity,
                       std::vector<int>& elementToPhysicalEntity,
                       bool verbose, bool insertBoundarySegments)
    {
      // register boundary segment to boundary segment factory for possible load balancing
      // this needs to be done on all cores since the type might not be known otherwise
      GmshReaderQuadraticBoundarySegment< Grid::dimension, Grid::dimensionworld >::registerFactory();

#ifndef NDEBUG
      // check that this method is called on all cores
      factory.comm().barrier();
#endif

      // create parse object and read grid on process 0
      if (factory.comm().rank() == 0)
      {
        GmshReaderParser<Grid> parser(factory,verbose,insertBoundarySegments);
        parser.read(fileName);

        boundarySegmentToPhysicalEntity = std::move(parser.boundaryIdMap());
        elementToPhysicalEntity = std::move(parser.elementIndexMap());
      }
      else
      {
        boundarySegmentToPhysicalEntity = {};
        elementToPhysicalEntity = {};
      }
    }

    //! pass a discarded lvalue argument to a function
    /**
     * This method is intended to be used in function calls that require
     * lvalue arguments, when the caller just wants to pass in temporary
     * variable that is immediately discarded after the return of the
     * function.  It expects an rvalue argument, that is turned into an
     * lvalue.  For instance:
     * ```c++
     * doRead(factory, fileName, discarded(std::vector<int>{}),
     *        discarded(std::vector<int>{}));
     * ```
     * Here, the vectors are constructed as rvalues, passed through
     * `discarded()` which turns them into lvalues, so they can be arguments
     * to `doRead()`.  `doRead()` will fill them with some data, and they
     * will be destroyed at the end of the full-expression containing the
     * function call.
     *
     * \note It is very likely an error to use this outside a function call
     *       argument.
     */
    template<class T>
    static T &discarded(T &&value) { return value; }

    struct DataArg {
      std::vector<int> *data_ = nullptr;
      DataArg(std::vector<int> &data) : data_(&data) {}
      DataArg(const decltype(std::ignore)&) {}
      DataArg() = default;
    };

    struct DataFlagArg : DataArg {
      bool flag_ = false;
      using DataArg::DataArg;
      DataFlagArg(bool flag) : flag_(flag) {}
    };

  public:
    typedef GridType Grid;

    /** \todo doc me
     *
     * \return The return type is a special pointer type that casts into Grid*,
     *    std::unique_ptr<Grid>, and std::shared_ptr<Grid>.  It is scheduled
     *    to be replaced by std::unique_ptr<Grid> eventually.
     */
    static std::unique_ptr<Grid> read (const std::string& fileName, bool verbose = true, bool insertBoundarySegments=true)
    {
      // make a grid factory
      Dune::GridFactory<Grid> factory;

      read(factory, fileName, verbose, insertBoundarySegments);

      return factory.createGrid();
    }

    /**
     * \brief Read Gmsh file, possibly with data
     * \param fileName                        Name of the file to read from.
     * \param boundarySegmentToPhysicalEntity Container to fill with boundary segment
     *                                        physical entity data (if insertBoundarySegments=true)
     * \param elementToPhysicalEntity         Container to fill with element physical entity data
     * \param verbose                         Whether to be chatty
     * \param insertBoundarySegments          Whether boundary segments are inserted into the factory
     *
     * \note When insertBoundarySegments=false there is no way to correctly use the values returned
     *       in boundarySegmentToPhysicalEntity. Make sure to set insertBoundarySegments=true if you
     *       intent to do this. An alternative is to use the other overloads which provide compile-time
     *       checking of the provided parameter combinations.
     *
     * \todo This interface is error-prone and should not be exposed to the user. However, the
     *       compile-time overloads may not provide sufficient runtime flexibility in all cases.
     *       Therefore this interface is kept until a better interface can be agreed on.
     *       See https://gitlab.dune-project.org/core/dune-grid/-/issues/107
     */
    static std::unique_ptr<Grid> read (const std::string& fileName,
                       std::vector<int>& boundarySegmentToPhysicalEntity,
                       std::vector<int>& elementToPhysicalEntity,
                       bool verbose = true, bool insertBoundarySegments=true)
    {
      // make a grid factory
      Dune::GridFactory<Grid> factory;

      doRead(
        factory, fileName, boundarySegmentToPhysicalEntity,
        elementToPhysicalEntity, verbose, insertBoundarySegments
      );

      return factory.createGrid();
    }

    /** \todo doc me */
    static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
                      bool verbose = true, bool insertBoundarySegments=true)
    {
      doRead(
        factory, fileName, discarded(std::vector<int>{}),
        discarded(std::vector<int>{}), verbose, insertBoundarySegments
      );
    }

    //! read Gmsh file, possibly with data
    /**
     * \param factory             The GridFactory to fill.
     * \param fileName            Name of the file to read from.
     * \param boundarySegmentData Container to fill with boundary segment
     *                            physical entity data, or `std::ignore`, or a
     *                            `bool` value.  Boundary segments are
     *                            inserted when a container or `true` is
     *                            given, otherwise they are not inserted.
     * \param elementData         Container to fill with element physical
     *                            entity data, or `std::ignore`.
     * \param verbose             Whether to be chatty.
     *
     * Containers to fill with data must be `std::vector<int>` lvalues.
     * Element data is indexed by the insertion index of the element,
     * boundarySegment data is indexed by the insertion index of the boundary
     * intersection.  These can be obtained from the `factory`, and are lost
     * once the grid gets modified (refined or load-balanced).
     *
     * \note At the moment the data containers are still filled internally,
     *       even if they are ignored.  So not having to pass them is more of
     *       a convenience feature and less of an optimization.  This may
     *       however change in the future.
     */
    static void read (Dune::GridFactory<Grid> &factory,
                      const std::string &fileName,
                      DataFlagArg boundarySegmentData,
                      DataArg elementData,
                      bool verbose=true)
    {
      doRead(
        factory, fileName,
        boundarySegmentData.data_
          ? *boundarySegmentData.data_ : discarded(std::vector<int>{}),
        elementData.data_
          ? *elementData.data_ : discarded(std::vector<int>{}),
        verbose,
        boundarySegmentData.flag_ || boundarySegmentData.data_
      );
    }

    /**
     * \brief Read Gmsh file, possibly with data
     * \param factory                         The GridFactory to fill.
     * \param fileName                        Name of the file to read from.
     * \param boundarySegmentToPhysicalEntity Container to fill with boundary segment
     *                                        physical entity data (if insertBoundarySegments=true)
     * \param elementToPhysicalEntity         Container to fill with element physical entity data
     * \param verbose                         Whether to be chatty
     * \param insertBoundarySegments          Whether boundary segments are inserted into the factory
     *
     * \note When insertBoundarySegments=false there is no way to correctly use the values returned
     *       in boundarySegmentToPhysicalEntity. Make sure to set insertBoundarySegments=true if you
     *       intent to do this. An alternative is to use the other overloads which provide compile-time
     *       checking of the provided parameter combinations.
     *
     * \todo This interface is error-prone and should not be exposed to the user. However, the
     *       compile-time overloads may not provide sufficient runtime flexibility in all cases.
     *       Therefore this interface is kept until a better interface can be agreed on.
     *       See https://gitlab.dune-project.org/core/dune-grid/-/issues/107
     */
    static void read (Dune::GridFactory<Grid>& factory,
                      const std::string& fileName,
                      std::vector<int>& boundarySegmentToPhysicalEntity,
                      std::vector<int>& elementToPhysicalEntity,
                      bool verbose, bool insertBoundarySegments)
    {
      doRead(
        factory, fileName, boundarySegmentToPhysicalEntity,
        elementToPhysicalEntity, verbose, insertBoundarySegments
      );
    }

    //! Dynamic Gmsh reader interface
    //\{

    [[deprecated("Will be removed after 2.8. Either use other constructors or use static methods without constructing an object")]]
    GmshReader() = default;

    using Opts = Gmsh::ReaderOptions;

    static constexpr Opts defaultOpts =
      Opts::verbose | Opts::insertBoundarySegments | Opts::readElementData | Opts::readBoundaryData;

    //! Construct a Gmsh reader object (alternatively use one of the static member functions)

    /**
     * \brief Construct a Gmsh reader object from a file name
     * \param fileName Name of the file to read from.
     * \param options Options of the type `Dune::Gmsh::ReaderOptions`
     *
     * To pass several options, combine them with the |-operator like this
     *
      \code
      using Opt = Dune::Gmsh::ReaderOptions;
      auto reader = Dune::GmshReader("grid.msh", Opt::verbose | Opt::readElementData)
      \endcode
     *
     * Per default the reader has enabled the following options
       - Dune::Gmsh::ReaderOptions::verbose
       - Dune::Gmsh::ReaderOptions::insertBoundarySegments
       - Dune::Gmsh::ReaderOptions::readBoundaryData
       - Dune::Gmsh::ReaderOptions::readElementData
     *
     * Passing any option to the interface will overwrite these defaults.
     *
     * A Dune grid object can be obtained via the `createGrid()` member
     */
    GmshReader(const std::string& fileName,
               Gmsh::ReaderOptions options = defaultOpts)
    {
      gridFactory_ = std::make_unique<Dune::GridFactory<Grid>>();
      readGridFile(fileName, *gridFactory_, options);
    }

    /**
     * \brief Construct a Gmsh reader object from a file name and a grid factory
     * \param fileName Name of the file to read from.
     * \param options Options of the type `Dune::Gmsh::ReaderOptions`
     * \note Use this constructor if you need access to the grid factor,
     *       e.g. for obtaining boundary segment insertion indices.
     */
    GmshReader(const std::string& fileName, GridFactory<Grid>& factory,
               Gmsh::ReaderOptions options = defaultOpts)
    {
      readGridFile(fileName, factory, options);
    }

    //! Access element data (maps element index to Gmsh physical entity)
    const std::vector<int>& elementData () const
    {
      checkElementData();
      return elementIndexToGmshPhysicalEntity_;
    }

    //! Access boundary data (maps boundary segment index to Gmsh physical entity)
    const std::vector<int>& boundaryData () const
    {
      checkBoundaryData();
      return boundarySegmentIndexToGmshPhysicalEntity_;
    }

    /**
     * \brief If element data is available
     * \note This is false if no such data was requested
     */
    bool hasElementData () const
    { return hasElementData_ && !extractedElementData_; }

    /**
     * \brief If boundary data is available
     * \note This is false if no such data was requested
     */
    bool hasBoundaryData () const
    { return hasBoundaryData_ && !extractedBoundaryData_; }

    //! Erase element data from reader and return the data
    std::vector<int> extractElementData ()
    {
      checkElementData();
      extractedElementData_ = true;
      return std::move(elementIndexToGmshPhysicalEntity_);
    }

    //! Erase boundary data from reader and return the data
    std::vector<int> extractBoundaryData ()
    {
      checkBoundaryData();
      extractedBoundaryData_ = true;
      return std::move(boundarySegmentIndexToGmshPhysicalEntity_);
    }

    //! Create the grid
    std::unique_ptr<Grid> createGrid ()
    {
      if (!gridFactory_)
        DUNE_THROW(Dune::InvalidStateException,
          "This GmshReader has been constructed with a Dune::GridFactory. "
          << "This grid factory has been filled with all information to create a grid. "
          << "Please use this factory to create the grid by calling factory.createGrid(). "
          << "Alternatively use the constructor without passing the factory in combination with this member function."
        );

      return gridFactory_->createGrid();
    }

    //\}

  private:
    void checkElementData () const
    {
      if (!hasElementData_)
        DUNE_THROW(Dune::InvalidStateException,
          "This GmshReader has been constructed without the option 'readElementData'. "
          << "Please enable reading element data by passing the option 'Gmsh::ReaderOpts::readElementData' "
          << "to the constructor of this class."
        );

      if (extractedElementData_)
        DUNE_THROW(Dune::InvalidStateException,
          "The element data has already been extracted from this GmshReader "
          << "via a function call to reader.extractElementData(). Use the extraced data or "
          << "read the grid data from file again by constructing a new reader."
        );
    }

    void checkBoundaryData () const
    {
      if (!hasBoundaryData_)
        DUNE_THROW(Dune::InvalidStateException,
          "This GmshReader has been constructed without the option 'readBoundaryData'. "
          << "Please enable reading boundary data by passing the option 'Gmsh::ReaderOpts::readBoundaryData' "
          << "to the constructor of this class."
        );

      if (extractedBoundaryData_)
        DUNE_THROW(Dune::InvalidStateException,
          "The boundary data has already been extracted from this GmshReader "
          << "via a function call to reader.extractBoundaryData(). Use the extraced data or "
          << "read the grid data from file again by constructing a new reader."
        );
    }

    void readGridFile (const std::string& fileName, GridFactory<Grid>& factory, Gmsh::ReaderOptions options)
    {
      const bool verbose = options & Opts::verbose;
      const bool insertBoundarySegments = options & Opts::insertBoundarySegments;
      const bool readBoundaryData = options & Opts::readBoundaryData;
      const bool readElementData = options & Opts::readElementData;

      doRead(
        factory, fileName, boundarySegmentIndexToGmshPhysicalEntity_,
        elementIndexToGmshPhysicalEntity_, verbose,
        readBoundaryData || insertBoundarySegments
      );

      // clear unwanted data
      if (!readBoundaryData)
          boundarySegmentIndexToGmshPhysicalEntity_ = std::vector<int>{};
      if (!readElementData)
          elementIndexToGmshPhysicalEntity_ = std::vector<int>{};

      hasElementData_ = readElementData;
      hasBoundaryData_ = readBoundaryData;
    }

    std::unique_ptr<Dune::GridFactory<Grid>> gridFactory_;

    std::vector<int> elementIndexToGmshPhysicalEntity_;
    std::vector<int> boundarySegmentIndexToGmshPhysicalEntity_;

    bool hasElementData_;
    bool hasBoundaryData_;

    // for better error messages, we keep track of these separately
    bool extractedElementData_ = false;
    bool extractedBoundaryData_ = false;
  };

  /** \} */

} // namespace Dune

#endif
