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

  // arbitrary dimension, implementation is in specialization
  template< int dimension, int dimWorld = dimension >
  class GmshReaderLinearBoundarySegment
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



  // arbitrary dimension, implementation is in specialization
  template<typename GridType, int dimension>
  class GmshReaderImp
  {};

  template< class GridType >
  struct GmshReaderImp< GridType, 2 >
  {
    static const int dim = GridType::dimension;
    static const int dimWorld = GridType::dimensionworld;

    dune_static_assert( (dimWorld <= 3), "GmshReader requires dimWorld <= 3." );

    typedef FieldVector< double, dimWorld > GlobalVector;

    typedef GmshReaderLinearBoundarySegment< dim, dimWorld > LinearBoundarySegment;
    typedef GmshReaderQuadraticBoundarySegment< dim, dimWorld > QuadraticBoundarySegment;

    template<typename T>
    static void read (Dune::GridFactory<GridType>& factory, const std::string& fileName,
                      bool verbose, bool insert_boundary_segments,
                      std::vector<T>& boundary_id_to_physical_entity,
                      std::vector<T>& element_index_to_physical_entity)
    {
      std::cout << "Reading 2d Gmsh grid..." << std::endl;

      // open file name, we use C I/O
      FILE* file = fopen(fileName.c_str(),"r");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

      // a read buffer
      char buf[512];

      //=========================================
      // Pass 1: Read vertices into vector
      //         Check vertices that are needed
      //         Renumber needed vertices
      //=========================================

      int boundary_element_count = 0;
      int element_count = 0;

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      int version_number, file_type, data_size;
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
        DUNE_THROW(Dune::IOError, "can only read version_number==2");
      if (verbose) std::cout << "version 2 Gmsh file detected" << std::endl;
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      int number_of_nodes;
      fscanf(file,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;
      std::vector< GlobalVector > nodes( number_of_nodes+1 ); // store positions
      for( int i = 1; i <= number_of_nodes; ++i )
      {
        int id;
        double x[ 3 ];
        if( fscanf( file, "%d %lg %lg %lg\n", &id, &x[ 0 ], &x[ 1 ], &x[ 2 ] ) != 4 )
          DUNE_THROW( Dune::IOError, "Unable to read vertex." );
        //          if (verbose) std::cout << id << " " << x << " " << y << " " << z << std::endl;
        if( id != i )
          DUNE_THROW( Dune::IOError, "Expected id " << i << "(got id " << id << "." );

        // just store node position
        for( int j = 0; j < dimWorld; ++j )
          nodes[ i ][ j ] = x[ j ];
      }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");

      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      fscanf(file,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl;
      unsigned int number_of_real_vertices=0;  // count number of vertices that are really needed
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity, mesh_partition;
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          fscanf(file,"%d",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if (k==3) mesh_partition = blub;
        }
        std::vector<int> simplexVertices(6);
        switch (elm_type)
        {
        case 1 :    // 2-node line
          simplexVertices.resize(2);
          fscanf(file,"%d %d\n",&(simplexVertices[0]),&(simplexVertices[1]));
          for (int i=0; i<2; i++)
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              //                   std::cout << "real vertex " << number_of_real_vertices << " at " << nodes[simplexVertices[i]] << " old number was " << simplexVertices[i] << std::endl;
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          boundary_element_count++;
          break;
        case 8 :    // 3-node line
          simplexVertices.resize(3);
          fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
          for (int i=0; i<2; i++)
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          boundary_element_count++;
          break;
        case 2 :    // 3-node triangle
          simplexVertices.resize(3);
          fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
          for (int i=0; i<3; i++)
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          element_count++;
          break;
        case 9 :    // 6-node triangle
          simplexVertices.resize(6);
          fscanf(file,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
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
          fgets(buf,512,file);     // skip rest of line if no triangle
        }
      }
      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      if (verbose) std::cout << "number of boundary elements = " << boundary_element_count << std::endl;
      if (verbose) std::cout << "number of elements = " << element_count << std::endl;
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");
      boundary_id_to_physical_entity.resize(boundary_element_count);
      element_index_to_physical_entity.resize(element_count);

      //==============================================
      // Pass 2: Insert boundary segments and elements
      //==============================================

      // go to beginning of file
      rewind(file);
      boundary_element_count = 0;
      element_count = 0;

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
        DUNE_THROW(Dune::IOError, "can only read version_number==2");
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      fscanf(file,"%d\n",&number_of_nodes);
      for (int i=1; i<=number_of_nodes; i++)
      {
        int id;
        double x,y,z;
        fscanf(file,"%d %lg %lg %lg\n",&id,&x,&y,&z);
        if (id!=i)
          DUNE_THROW(Dune::IOError, "expected id " << i);
      }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");

      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      fscanf(file,"%d\n",&number_of_elements);
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity, mesh_partition;
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          fscanf(file,"%d",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if (k==3) mesh_partition = blub;
        }
        std::vector<int> simplexVertices(6);
        std::vector<unsigned int> vertices(3);
        switch (elm_type)
        {
        case 1 :    // 2-node line

          // read the vertices, but don't do anything with them.  Linear boundary
          // segments are the default.
          simplexVertices.resize(2);
          fscanf(file,"%d %d\n",&(simplexVertices[0]),&(simplexVertices[1]));
          vertices.resize(2);
          for (int i=0; i<2; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices
          if (insert_boundary_segments)
            factory.insertBoundarySegment(vertices, new LinearBoundarySegment(nodes[simplexVertices[0]],nodes[simplexVertices[1]]));
          boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
          boundary_element_count++;
          break;
        case 8 :    // 3-node line
          simplexVertices.resize(3);
          fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
          vertices.resize(2);
          for (int i=0; i<2; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices
          if (insert_boundary_segments)
            factory.insertBoundarySegment(vertices, new QuadraticBoundarySegment(nodes[simplexVertices[0]],nodes[simplexVertices[2]],nodes[simplexVertices[1]]));
          boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
          boundary_element_count++;
          break;
        case 2 :    // 3-node triangle
          simplexVertices.resize(3);
          fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
          vertices.resize(3);
          for (int i=0; i<3; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          element_index_to_physical_entity[element_count] = physical_entity;
          element_count++;
          break;
        case 9 :    // 6-node triangle
          simplexVertices.resize(6);
          fscanf(file,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
          vertices.resize(3);
          for (int i=0; i<3; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          element_index_to_physical_entity[element_count] = physical_entity;
          element_count++;
          break;
        default :
          fgets(buf,512,file);     // skip rest of line if no tetrahedron
        }
      }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);
    }
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

      //         alpha = beta = 0.5;
      //         gamma = 0.5*sqrt2;
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



  template<typename GridType>
  class GmshReaderImp<GridType,3>
  {
  public:

    template<typename T>
    static void read (Dune::GridFactory<GridType>& factory, const std::string& fileName,
                      bool verbose, bool insert_boundary_segments,
                      std::vector<T>& boundary_id_to_physical_entity,
                      std::vector<T>& element_index_to_physical_entity)
    {
      // the grid dimension
      const int dim = GridType::dimension;

      // open file name, we use C I/O
      FILE* file = fopen(fileName.c_str(),"r");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

      // a read buffer
      char buf[512];

      //=========================================
      // Pass 1: Read vertices into vector
      //         Check vertices that are needed
      //         Renumber needed vertices
      //=========================================

      int boundary_element_count = 0;
      int element_count = 0;

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      int version_number, file_type, data_size;
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
        DUNE_THROW(Dune::IOError, "can only read version_number==2");
      if (verbose) std::cout << "version 2 Gmsh file detected" << std::endl;
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      int number_of_nodes;
      fscanf(file,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;
      std::vector<Dune::FieldVector<double,dim> > nodes(number_of_nodes+1); // store positions
      for (int i=1; i<=number_of_nodes; i++) // note: ids are starting with 1!
      {
        int id;
        double x,y,z;
        fscanf(file,"%d %lg %lg %lg\n",&id,&x,&y,&z);
        //          if (verbose) std::cout << id << " " << x << " " << y << " " << z << std::endl;
        if (id!=i)
          DUNE_THROW(Dune::IOError, "expected id " << i);

        Dune::FieldVector<double,dim> position;
        position[0] = x; position[1] = y; position[2] = z;
        nodes[i] = position;   // just store node position
      }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");

      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      fscanf(file,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl;
      unsigned int number_of_real_vertices=0;  // count number of vertices that are really needed
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++) // note: ids are starting with 1!
      {
        int id, elm_type, number_of_tags;
        fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity, mesh_partition;
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          fscanf(file,"%d",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if (k==3) mesh_partition = blub;
        }
        std::vector<int> simplexVertices(10);
        switch (elm_type)
        {
        case 2 :    // 3-node triangle
          simplexVertices.resize(3);
          fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
          for (size_t i=0; i<simplexVertices.size(); i++)
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          boundary_element_count++;
          break;

        case 4 :    // 4-node tetrahedron
          simplexVertices.resize(4);
          fscanf(file,"%d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),&(simplexVertices[3]));
          for (size_t i=0; i<simplexVertices.size(); i++)
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          element_count++;
          break;

        case 9 :    // 6-node triangle
          simplexVertices.resize(6);
          fscanf(file,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
          for (int i=0; i<3; i++)     // insert only the first three !
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          boundary_element_count++;
          break;

        case 11 :    // 10-node tetrahedron
          simplexVertices.resize(10);
          fscanf(file,"%d %d %d %d %d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]),
                 &(simplexVertices[6]),&(simplexVertices[7]),&(simplexVertices[8]),&(simplexVertices[9]));
          for (int i=0; i<4; i++)     // insert only the first four !
            if (renumber.find(simplexVertices[i])==renumber.end())
            {
              renumber[simplexVertices[i]] = number_of_real_vertices++;
              factory.insertVertex(nodes[simplexVertices[i]]);
            }
          element_count++;
          break;

        default :
          fgets(buf,512,file);     // skip rest of line if no tetrahedron
        }
      }
      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      if (verbose) std::cout << "number of boundary elements = " << boundary_element_count << std::endl;
      if (verbose) std::cout << "number of elements = " << element_count << std::endl;
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");
      boundary_id_to_physical_entity.resize(boundary_element_count);
      element_index_to_physical_entity.resize(element_count);

      //==============================================
      // Pass 2: Insert boundary segments and elements
      //==============================================

      // go to beginning of file
      rewind(file);
      boundary_element_count = 0;
      element_count = 0;

      // process header
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      fscanf(file,"%d %d %d\n",&version_number,&file_type,&data_size);
      if (version_number!=2)
        DUNE_THROW(Dune::IOError, "can only read version_number==2");
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // node section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      fscanf(file,"%d\n",&number_of_nodes);
      for (int i=1; i<=number_of_nodes; i++)
      {
        int id;
        double x,y,z;
        fscanf(file,"%d %lg %lg %lg\n",&id,&x,&y,&z);
        if (id!=i)
          DUNE_THROW(Dune::IOError, "expected id " << i);
      }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndNodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndNodes");

      // element section
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      fscanf(file,"%d\n",&number_of_elements);

      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        fscanf(file,"%d %d %d",&id,&elm_type,&number_of_tags);
        int physical_entity, elementary_entity, mesh_partition;
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          fscanf(file,"%d",&blub);
          if (k==1) physical_entity = blub;
          if (k==2) elementary_entity = blub;
          if (k==3) mesh_partition = blub;
        }
        std::vector<int> simplexVertices(10);
        std::vector<unsigned int> vertices(10);
        switch (elm_type)
        {
        case 2 :    // 3-node triangle

          // Read the vertices but don't do anything with them.  Linear boundary segments
          // are the default anyways.
          simplexVertices.resize(3);
          fscanf(file,"%d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]));
          vertices.resize(3);
          for (int i=0; i<3; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices

          if (insert_boundary_segments)
            factory.insertBoundarySegment(vertices,new GmshReaderLinearBoundarySegment<3>(nodes[simplexVertices[0]],nodes[simplexVertices[1]],nodes[simplexVertices[2]]));
          boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
          boundary_element_count++;
          break;

        case 9 :    // 6-node triangle
          simplexVertices.resize(6);
          fscanf(file,"%d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]));
          vertices.resize(3);
          for (int i=0; i<3; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices first three vertices

          if (insert_boundary_segments)
            factory.insertBoundarySegment(vertices,new GmshReaderQuadraticBoundarySegment<3>(nodes[simplexVertices[0]],nodes[simplexVertices[1]],nodes[simplexVertices[2]],nodes[simplexVertices[3]],nodes[simplexVertices[4]],nodes[simplexVertices[5]]));
          boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
          boundary_element_count++;
          break;

        case 4 :    // 4-node tetrahedron
          simplexVertices.resize(4);
          fscanf(file,"%d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),&(simplexVertices[3]));
          vertices.resize(4);
          for (int i=0; i<4; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          element_index_to_physical_entity[element_count] = physical_entity;
          element_count++;
          break;

        case 11 :    // 10-node tetrahedron
          simplexVertices.resize(10);
          fscanf(file,"%d %d %d %d %d %d %d %d %d %d\n",&(simplexVertices[0]),&(simplexVertices[1]),&(simplexVertices[2]),
                 &(simplexVertices[3]),&(simplexVertices[4]),&(simplexVertices[5]),
                 &(simplexVertices[6]),&(simplexVertices[7]),&(simplexVertices[8]),&(simplexVertices[9]));
          vertices.resize(4);
          for (int i=0; i<4; i++)
            vertices[i] = renumber[simplexVertices[i]];     // renumber vertices
          factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),vertices);
          element_index_to_physical_entity[element_count] = physical_entity;
          element_count++;
          break;

        default :
          fgets(buf,512,file);     // skip rest of line if no tetrahedron
        }
      }
      fscanf(file,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);
    }
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

      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;

      GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,verbose,insert_boundary_segments,
                                                        boundary_id_to_physical_entity,
                                                        element_index_to_physical_entity);

      return factory.createGrid();
    }

    /** \todo doc me */
    template<typename T>
    static GridType* read (const std::string& fileName, std::vector<T>& boundary_id_to_physical_entity,
                           std::vector<T>& element_index_to_physical_entity, bool verbose = true,
                           bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<GridType> factory;

      GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,verbose,insert_boundary_segments,
                                                        boundary_id_to_physical_entity,
                                                        element_index_to_physical_entity);
      return factory.createGrid();
    }

    /** \todo doc me */
    static GridType* read (GridType& grid, const std::string& fileName,
                           bool verbose = true, bool insert_boundary_segments=true)
    {
      // make a grid factory
      Dune::GridFactory<GridType> factory(&grid);

      // store dummy
      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;

      GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,verbose,insert_boundary_segments,
                                                        boundary_id_to_physical_entity,
                                                        element_index_to_physical_entity);
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

      GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,verbose,insert_boundary_segments,
                                                        boundary_id_to_physical_entity,
                                                        element_index_to_physical_entity);
      return factory.createGrid();
    }


    /** \todo doc me */
    static void read (Dune::GridFactory<GridType>& factory, const std::string& fileName, bool verbose = true,
                      bool insert_boundary_segments=true)
    {
      // store dummy
      std::vector<int> boundary_id_to_physical_entity;
      std::vector<int> element_index_to_physical_entity;

      GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,verbose,insert_boundary_segments,
                                                        boundary_id_to_physical_entity,
                                                        element_index_to_physical_entity);
    }

    /** \todo doc me */
    template<typename T>
    static void read (Dune::GridFactory<GridType>& factory,
                      const std::string& fileName,
                      std::vector<T>& boundary_id_to_physical_entity,
                      std::vector<T>& element_index_to_physical_entity,
                      bool verbose = true, bool insert_boundary_segments=true)
    {
      GmshReaderImp<GridType,GridType::dimension>::read(factory,fileName,verbose,insert_boundary_segments,
                                                        boundary_id_to_physical_entity,
                                                        element_index_to_physical_entity);
    }
  };

  /** \} */

} // namespace Dune

#endif
