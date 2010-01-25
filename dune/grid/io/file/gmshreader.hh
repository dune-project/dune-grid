// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GMSHREADER_HH
#define DUNE_GMSHREADER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <cstring>
#include <cstdarg>

#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundarysegment.hh>

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
    class GmshReaderLinearBoundarySegment
    {};

    // arbitrary dimension, implementation is in specialization
    template< int dimension, int dimWorld = dimension >
    class GmshReaderQuadraticBoundarySegment
    {};

    // linear boundary segments in 1d
    template< int dimWorld >
    struct GmshReaderLinearBoundarySegment< 2, dimWorld >
      : public Dune::BoundarySegment< 2, dimWorld >
    {
      typedef Dune::FieldVector< double, dimWorld > GlobalVector;

      GmshReaderLinearBoundarySegment ( const GlobalVector &p0, const GlobalVector &p1 )
        : x_( p0 ), p_( p1 - p0 )
      {}

      virtual GlobalVector operator() ( const Dune::FieldVector< double, 1 > &local ) const
      {
        GlobalVector y = x_;
        y.axpy( local[ 0 ], p_ );
        return y;
      }

    private:
      GlobalVector x_, p_;
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

    // linear boundary segments in 2d
    template<>
    struct GmshReaderLinearBoundarySegment< 3, 3 >
      : public Dune::BoundarySegment< 3 >
    {
      GmshReaderLinearBoundarySegment (Dune::FieldVector<double,3> p0_, Dune::FieldVector<double,3> p1_,
                                       Dune::FieldVector<double,3> p2_)
        : p0(p0_), p1(p1_), p2(p2_)
      {
        //        std::cout << "created boundary segment " << p0 << " | " << p1 << " | " << p2 << std::endl;
      }

      virtual Dune::FieldVector<double,3> operator() (const Dune::FieldVector<double,2>& local) const
      {
        Dune::FieldVector<double,3> y;
        y = 0.0;
        y.axpy(1.0-local[0]-local[1],p0);
        y.axpy(local[0],p1);
        y.axpy(local[1],p2);
        //        std::cout << "eval boundary segment local=" << local << " y=" << y << std::endl;
        return y;
      }

    private:
      Dune::FieldVector<double,3> p0,p1,p2;
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

  //! Parser for Gmsh data. Implementations must be (partially) deimension dependent
  template<typename GridType, int dimension>
  class GmshReaderParser
  {};

  //! dimension independent parts for GmshReaderParser
  template<typename GridType, typename DimImp>
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
    int line;
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

    // dimension dependent routines
    void pass1HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes);

    void pass2HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes,
                            const int physical_entity);

    void readfile(FILE * file, int cnt, const char * format, ...)
    {
      va_list vargs;
      va_start(vargs,format);
      int c = vfscanf(file, format, vargs);
      if (c != cnt)
        DUNE_THROW(Dune::IOError, "Error parsing " << fileName << " line " << line
                                                   << ": Expected '" << format << "', only read " << c << " entries instead of " << cnt << ".");
      va_end(vargs);
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
      line = 0;

      // process header
      int version_number, file_type, data_size;

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      readfile(file,3,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
        DUNE_THROW(Dune::IOError, "can only read version_number==2");
      if (verbose) std::cout << "version 2 Gmsh file detected" << std::endl;
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
        readfile(file,3,"%d %d %d",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity, mesh_partition;
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if (k==3) mesh_partition = blub;
        }
        static_cast<DimImp*>(this)->pass1HandleElement(file, elm_type, renumber, nodes);
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
        readfile(file,3,"%d %d %d",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity, mesh_partition;
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if (k==3) mesh_partition = blub;
        }
        static_cast<DimImp*>(this)->pass2HandleElement(file, elm_type, renumber, nodes, physical_entity);
      }
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);
    }
  };

  //! Parser for Gmsh data. Specialization for 2D
  template<typename GridType>
  class GmshReaderParser<GridType, 2> : public GmshReaderParserBase< GridType, GmshReaderParser<GridType, 2> >
  {
    typedef GmshReaderParserBase< GridType, GmshReaderParser<GridType, 2> > Base;
    typedef typename Base::GlobalVector GlobalVector;
    using Base::dim;
    using Base::dimWorld;
    using Base::buf;
    using Base::factory;
    using Base::verbose;
    using Base::element_index_to_physical_entity;
    using Base::element_count;
    using Base::number_of_real_vertices;
    using Base::boundary_element_count;
    using Base::insert_boundary_segments;
    using Base::boundary_id_to_physical_entity;
    using Base::readfile;

    typedef GmshReaderLinearBoundarySegment< dim, dimWorld > LinearBoundarySegment;
    typedef GmshReaderQuadraticBoundarySegment< dim, dimWorld > QuadraticBoundarySegment;

    friend class GmshReaderParserBase< GridType, GmshReaderParser<GridType, 2> >;

    // dimension dependent routines
    void pass1HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes)
    {
      std::vector<int> simplexVertices(6);
      switch (elm_type)
      {
      case 1 :          // 2-node line
        simplexVertices.resize(2);
        readfile(file,2,"%d %d\n",&(simplexVertices[0]),&(simplexVertices[1]));
        for (int i=0; i<2; i++)
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        boundary_element_count++;
        break;
      case 8 :          // 3-node line
        simplexVertices.resize(3);
        readfile(file,3,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
        for (int i=0; i<2; i++)
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        boundary_element_count++;
        break;
      case 2 :          // 3-node triangle
        simplexVertices.resize(3);
        readfile(file,3,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
        for (int i=0; i<3; i++)
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        element_count++;
        break;
      case 9 :          // 6-node triangle
        simplexVertices.resize(6);
        readfile(file,6,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
        for (int i=0; i<3; i++)
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        element_count++;
        break;
      default :
        fgets(buf,512,file);             // skip rest of line if no triangle
      }
    }

    void pass2HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes,
                            const int physical_entity)
    {
      std::vector<int> simplexVertices(6);
      std::vector<unsigned int> vertices(3);
      switch (elm_type)
      {
      case 1 :          // 2-node line
        // read the vertices, but don't do anything with them.  Linear boundary
        // segments are the default.
        simplexVertices.resize(2);
        readfile(file,2,"%d %d\n",&(simplexVertices[0]),&(simplexVertices[1]));
        vertices.resize(2);
        for (int i=0; i<2; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices
        if (insert_boundary_segments)
          factory.insertBoundarySegment(vertices, shared_ptr<BoundarySegment<dim,dimWorld> >(new LinearBoundarySegment(nodes[simplexVertices[0]],nodes[simplexVertices[1]])));
        boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
        boundary_element_count++;
        break;
      case 8 :          // 3-node line
        simplexVertices.resize(3);
        readfile(file,3,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
        vertices.resize(2);
        for (int i=0; i<2; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices
        if (insert_boundary_segments)
          factory.insertBoundarySegment(vertices, shared_ptr<BoundarySegment<dim,dimWorld> >(new QuadraticBoundarySegment(nodes[simplexVertices[0]],nodes[simplexVertices[2]],nodes[simplexVertices[1]])));
        boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
        boundary_element_count++;
        break;
      case 2 :          // 3-node triangle
        simplexVertices.resize(3);
        readfile(file,3,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
        vertices.resize(3);
        for (int i=0; i<3; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
        element_index_to_physical_entity[element_count] = physical_entity;
        element_count++;
        break;
      case 9 :          // 6-node triangle
        simplexVertices.resize(6);
        readfile(file,6,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
        vertices.resize(3);
        for (int i=0; i<3; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
        element_index_to_physical_entity[element_count] = physical_entity;
        element_count++;
        break;
      default :
        fgets(buf,512,file);             // skip rest of line if no tetrahedron
      }
    }
  public:
    GmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      Base(_factory,v,i) {}
  };

  //! Parser for Gmsh data. Specialization for 3D
  template<typename GridType>
  class GmshReaderParser<GridType, 3> : public GmshReaderParserBase< GridType, GmshReaderParser<GridType, 3> >
  {
    typedef GmshReaderParserBase< GridType, GmshReaderParser<GridType, 3> > Base;
    typedef typename Base::GlobalVector GlobalVector;
    using Base::dim;
    using Base::dimWorld;
    using Base::buf;
    using Base::factory;
    using Base::verbose;
    using Base::element_index_to_physical_entity;
    using Base::element_count;
    using Base::number_of_real_vertices;
    using Base::boundary_element_count;
    using Base::insert_boundary_segments;
    using Base::boundary_id_to_physical_entity;
    using Base::readfile;

    typedef GmshReaderLinearBoundarySegment< dim, dimWorld > LinearBoundarySegment;
    typedef GmshReaderQuadraticBoundarySegment< dim, dimWorld > QuadraticBoundarySegment;

    friend class GmshReaderParserBase< GridType, GmshReaderParser<GridType, 3> >;

    void pass1HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes)
    {
      std::vector<int> simplexVertices(10);
      switch (elm_type)
      {
      case 2 :          // 3-node triangle
        simplexVertices.resize(3);
        readfile(file,3,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
        for (size_t i=0; i<simplexVertices.size(); i++)
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        boundary_element_count++;
        break;
      case 4 :          // 4-node tetrahedron
        simplexVertices.resize(4);
        readfile(file,4,"%d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),&(simplexVertices[3]));
        for (size_t i=0; i<simplexVertices.size(); i++)
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        element_count++;
        break;
      case 9 :          // 6-node triangle
        simplexVertices.resize(6);
        readfile(file,6,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
        for (int i=0; i<3; i++)             // insert only the first three !
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        boundary_element_count++;
        break;
      case 11 :          // 10-node tetrahedron
        simplexVertices.resize(10);
        readfile(file,10,"%d %d %d %d %d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]),
                 &(simplexVertices[6]),&(simplexVertices[7]),&(simplexVertices[8]),&(simplexVertices[9]));
        for (int i=0; i<4; i++)             // insert only the first four !
          if (renumber.find(simplexVertices[i])==renumber.end())
          {
            renumber[simplexVertices[i]] = number_of_real_vertices++;
            factory.insertVertex(nodes[simplexVertices[i]]);
          }
        element_count++;
        break;
      default :
        fgets(buf,512,file);             // skip rest of line if no tetrahedron
      }
    }

    void pass2HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes,
                            const int physical_entity)
    {
      std::vector<int> simplexVertices(10);
      std::vector<unsigned int> vertices(10);
      switch (elm_type)
      {
      case 2 :          // 3-node triangle
        // Read the vertices but don't do anything with them.  Linear boundary segments
        // are the default anyways.
        simplexVertices.resize(3);
        readfile(file,3,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
        vertices.resize(3);
        for (int i=0; i<3; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices

        if (insert_boundary_segments)
          factory.insertBoundarySegment(vertices,shared_ptr<BoundarySegment<dim,dimWorld> >(new LinearBoundarySegment(nodes[simplexVertices[0]],nodes[simplexVertices[1]],nodes[simplexVertices[2]])));
        boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
        boundary_element_count++;
        break;
      case 9 :          // 6-node triangle
        simplexVertices.resize(6);
        readfile(file,6,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
        vertices.resize(3);
        for (int i=0; i<3; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices first three vertices
        if (insert_boundary_segments)
          factory.insertBoundarySegment(vertices,shared_ptr<BoundarySegment<dim,dimWorld> >(new QuadraticBoundarySegment(nodes[simplexVertices[0]],nodes[simplexVertices[1]],nodes[simplexVertices[2]],nodes[simplexVertices[3]],nodes[simplexVertices[4]],nodes[simplexVertices[5]])));
        boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
        boundary_element_count++;
        break;
      case 4 :          // 4-node tetrahedron
        simplexVertices.resize(4);
        readfile(file,4,"%d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),&(simplexVertices[3]));
        vertices.resize(4);
        for (int i=0; i<4; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
        element_index_to_physical_entity[element_count] = physical_entity;
        element_count++;
        break;
      case 11 :          // 10-node tetrahedron
        simplexVertices.resize(10);
        readfile(file,10,"%d %d %d %d %d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]),
                 &(simplexVertices[6]),&(simplexVertices[7]),&(simplexVertices[8]),&(simplexVertices[9]));
        vertices.resize(4);
        for (int i=0; i<4; i++)
          vertices[i] = renumber[simplexVertices[i]];               // renumber vertices
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
        element_index_to_physical_entity[element_count] = physical_entity;
        element_count++;
        break;
      default :
        fgets(buf,512,file);             // skip rest of line if no tetrahedron
      }
    }
  public:
    GmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      Base(_factory,v,i) {}
  };

  /**
     \ingroup Gmsh

     \brief Read Gmsh mesh file

     Read a .msh file generated using Gmsh and construct a grid using the grid factory interface.
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
    template<typename T>
    static GridType* read (const std::string& fileName,
                           std::vector<T>& boundary_id_to_physical_entity,
                           std::vector<T>& element_index_to_physical_entity,
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
    template<typename T>
    static GridType* read (GridType& grid, const std::string& fileName,
                           std::vector<T>& boundary_id_to_physical_entity,
                           std::vector<T>& element_index_to_physical_entity,
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
    template<typename T>
    static void read (Dune::GridFactory<GridType>& factory,
                      const std::string& fileName,
                      std::vector<T>& boundary_id_to_physical_entity,
                      std::vector<T>& element_index_to_physical_entity,
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
