// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MACROGRIDPARSER_HH
#define DUNE_MACROGRIDPARSER_HH

#include <iostream>
#include <fstream>

#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <assert.h>
#include <cmath>

#include "dgfparserblocks.hh"

/**
 * @file
 * @brief  Classes for constructing macrogrids from files in the
 *         Dune Grid Format (DGF)
 *
 * Usage: the class Dune::MacroGrid returns a
 * pointer to a initialized grid. For details on the DGF file format
 * see the @link DuneGridFormatParser macro grid parser@endlink.
 * Note that the caller must explicitly delete the generated grid!
 *
 * The parsing of the DGF file is performed through the Dune::DuneGridFormatParser
 * class.
 *
 * @author Andreas Dedner
 */

namespace Dune {
  //! \brief The %DuneGridFormatParser class: reads a DGF file and stores
  //! build information in vector structures used by the MacroGrid class.
  class DuneGridFormatParser {
  public:
    //! default constructor which does nothing
    DuneGridFormatParser() :  vtx(0), elements(0) , bound(0) {}
    typedef enum {Simplex,Cube,General} element_t;
    typedef enum {counterclockwise=1,clockwise=-1} orientation_t;
    //! \brief method which reads the dgf file
    //!
    //! fills the vtx,element, and bound vectors
    //! \return: \n
    //! -1: no DGF file (prehaps a native file format) \n
    //!  0: some error during parsing of file \n
    //!  1: everything fine
    inline int readDuneGrid(std::istream&);
    //! method to write macrogridfiles in alu format (cam be used without dune)
    inline void writeAlu(std::ostream&);
    //! method to write macrogridfiles in alberta format (cam be used without dune)
    inline void writeAlberta(std::ostream&);
  protected:
    // dimension of world and problem: set through the readDuneGrid() method
    int dimw;
    // vector of vertex coordinates
    std::vector < std::vector <double> > vtx;
    int nofvtx;
    // vector of elements
    std::vector < std::vector <int> > elements;
    int nofelements;
    // vector of boundary segments + identifier
    std::vector < std::vector <int> > bound;
    int nofbound;
    // map to generate and find boundary segments
    std::map<EntityKey<int>,int> facemap;
    // set by generator depending on element type wanted
    element_t element;
    // set by the readDuneGrid method depending
    // on what type the elements were generated
    bool simplexgrid;
    // true if grid is generated using the intervall Block
    bool isInterval;
    // call to tetgen/triangle
    inline int generateSimplexGrid(std::istream&);
    inline int readTetgenTriangle(std::string);
    // helper methods
    inline void setOrientation(int fixvtx,orientation_t orientation=counterclockwise);
    inline void setRefinement(int);
    double testTriang(int snr);
  };
  //! \brief Class for constructing grids from DGF files.
  //!
  //! The constructor of the class is given the filename of the DGF file and
  //! the class exports a template method
  //! operator GridType* () which
  //! reads a file in the DGF file and constructes a pointer to a
  //! GridType instance. A pointer to a grid of type GridType is constructed
  //! as follows:
  //! @code
  //! GridType* grid = Dune::MacroGrid(filename,MPI_COMM_WORLD);
  //! @endcode
  class MacroGrid : protected DuneGridFormatParser {
  public:
    //! constructor given the name of a DGF file
    MacroGrid(const char* filename,int MPICOMM=-1) :
      DuneGridFormatParser(),
      filename_(filename),
      MPICOMM_(MPICOMM) {}
    //! conversion method from a DGF file to a GridType* instance
    template <class GridType>
    inline operator GridType* () {
      return Impl<GridType>::generate(*this,filename_,MPICOMM_);
    }
  private:
    //! class with one static member which has to be defined for each
    //! grid implementation
    template <class GT>
    class Impl {
    public:
      static GT* generate(MacroGrid& mg,const char* filename,int MPICOMM=-1);
    };
    const char* filename_;
    int MPICOMM_;
  };
  template <class GridType>
  class GridPtr : private MacroGrid {
  public:
    //! constructor given the name of a DGF file
    GridPtr(const char* filename,int MPICOMM=-1) :
      MacroGrid(filename,MPICOMM),
      grid_(*this) {}
    ~GridPtr() {
      delete grid_;
    }
    GridType& operator*() {
      return *grid_;
    }
    GridType* operator->() {
      return grid_;
    }
  private:
    GridType* grid_;
  };
}
#include "dgfparser.cc"

/*!
     @defgroup DuneGridFormatParser The Dune Grid Format (DGF)
     @ingroup Grid
     @brief Classes for reading a macrogrid file in the dune
     macrogrid format (dgf)

     @section General
     <!--=========-->
     The DGF format allows the simple description of macrogrids, which
     can be used to construct Dune grids independent of the underlying
     implementation. Due to the generality of the approach only a subset of
     the description language for each grid implementation is available through
     this interface.

     @section Usage
     <!---------------------------------------------->
     There are two ways of constructing Dune grids using DGF files:

     -# By including the file gridtype.hh from the
        dune/grid/utility/macrogridparser
        directory: \n
        by defining one of the symbols
        \c ALBERTAGRID ,
        \c ALUGRID_CUBE ,
        \c ALUGRID_SIMPLEX ,
        \c SGRID , or
        \c YASPGRID
        and the integer
        \c DUNE_PROBLEM_DIM
        one obtains a definition of the type \c GridType and
        the required header files for the desired grid and the macrogrid parser
        are included.
     -# By directly including one of the files macrogrid*.hh
        from the dune/grid/io/file/dgfparser directory: \n
        in this case only the required header files for the
        desired grid and the macrogrid parser are included but no typedef is made.

     After a grid type (denoted with \c GridType in the following)
     is selected in some way, the grid can be constructed either by calling
     @code
     GridType* grid = Dune::MacroGrid(filename,MPI_COMM_WORLD);
     @endcode
     here \c filename is the name of the dgf file. The method returns a
     pointer to a \c GridType instance whose macrogrid is described through the
     dgf file.

     Remarks:
     -# The last argument \c MPI_COMM_WORLD
        defaults to -1, a value which can be used in a serial run.
        If gridtype.hh is included, then macros
        \c MPISTART , \c MPIEND  and the integer variable
        \c MPI_COMM_WORLD
        are defined - in the serial case the macros are empty whereas
        in the parallel case, i.e., when \c HAVE_MPI_CPP is set,
        \c MPI_INIT(&argc, &argv)
        is called through \c MPISTART  and in addition
        \c myrank  and \c mysize
        are set through the corresponding MPI commands. See gridtype.hh
        for more details.
     -# If the file given through the first argument is not a dgf file
        a suitable constructure on the \c GridType class is called - if
        one is available.
     -# The caller must free the allocated memory by calling \c delete \c grid

     @section FORMAT Format Description
     <!--=========-->
     We assume in the following that the type
     \c GridType  is suitably defined, denoting a dune grid type
     in \c dimworld  space dimension. In the following
     a point in dimworld space is called a vector.

     In general dgf files consists of one or more blocks, each starting with a
     keyword and ending with a line starting with a # symbol.
     In a block each full line is parsed from the beginning up to the
     first occupance of a \% symbol, which can be used to include comments.
     Trailing whitespaces are ignored
     during line parsing. Also the parsing is not case sensitive.

     Some example files are given below (\ref EXAMPLES).

     @subsection START First line
     <!---------------------------------------------->
     DGF files must start with the keyword \b DGF. Files passed to the
     grid parser not starting with this keyword are directly passed to a
     suitable constructure in the GridType class.

     @subsection BLOCKS Blocks
     <!---------------------------------------------->
     In the following all blocks are briefly described in alphabetical order.
     The construction of the grid is detailed below. In general lines
     are parser sequentially, with the exception of lines starting
     with a special keyword which are treated separately.

     - \b Boundarydomain  \n
       Each line consists of an integer greater than zero and two
       vectors describing an interval in \c dimworld space. The first entry of
       the line is a boundary identifier.
       - A special keyword is
         \b default
         followed by a positive integer which is to be used
         as a default boundary identifier.
     - \b Boundarysegments \n
       Each line consists of a positive integer denoting a boundary
       id and a set of vertex indices (see \b Vertex block)
       describing one boundary patch of the macrogrid.
       \e Note: so far only simplex grid patches can be prescribed in this way.
     - \b Cube  \n
       Each line consists of \c dimworld<SUP>2</SUP> vertex indices (see \b Vertex block)
       describing one cube
       of the macrogrid. \e Note that the ordering of the local vertices
       has to follow the %Dune
       @link GridReferenceElements reference elements@endlink.
     - \b Interval \n
       The first two lines give the lower and upper vector of an interval,
       the third line consists of \c dimworld integers. used to determine the
       initial partition of the interval into equal sized cubes.
     - \b Simples \n
       Each line consists of \c dimworld+1 vertex indices (see \b Vertex block)
       describing one simplex
       of the macrogrid.
       \e Note that in 2d no ordering of local vertices is required, in
       3d each should conform with the %Dune
       @link GridReferenceElements reference elements@endlink.
     - \b Simplexgenerator \n
       Using this block a simplex grid can be automatically generated using
       one of the freely available grid generation tools
       Tetgen (http://tetgen.berlios.de) for \c dimworld=3 or
       Triangle (http://www.cs.cmu.edu/~quake/triangle.html) for \c dimworld=2.
       For more detail see \ref Simplexgeneration.
      .
     - \b Vertex \n
       Each line consists of a vector representing a vertex of the
       macrogrid. The vertices are consecutively numbered starting from zero -
       these numbers are called vertex indices and can be referenced in the
       \b Simplex, \b Cube, and \b Boundarysegment blocks.

     @section CONSTR The Grid Construction Process
     <!---------------------------------------------->
     The construction of the grid depends on the type of elements the grid
     given by \c GridType can handle.

     @subsection CONSTRCART Cartesian grids
     (Dune::SGrid , or
      Dune::YaspGrid )

     The grid is constructed using only the information from the
     \b Interval  block.

     @subsection CONSTRSIMPL Simplex grids
     (Dune::AlbertaGrid, or
      Dune::ALUSimplexGrid<3,3>, and
      Dune::ALUSimplexGrid<2,2>)

     The verticies and elements of the grid are constructed in one of the
     following four stages:
     -# The parser firsts looks for the \b Simplex generation  block.
        If this block is found, a grid is generated via Triangle/Tetgen
        using as  verticies the union of the verticies defined by the
        \b Interval  and the \b Vertex  block.
     -# If no \b Simplexgeneration  block is present,
        the file is parser for
        an \b Interval  block; if present a Cartesian grid is build and
        each element is partitioned either into two triangles or into six
        tetrahedron.
     -# If neither a \b Simplexgeneration  nor a \b Interval
        block is found, the grid is generated using the information
        from the \b Vertex  and the \b Simplex  blocks.
        Note that no specific ordering of the local numbering of each simplex
        is required.
     -# If so far no element information has been generated, the parser
        tries to first generate a cube grid using information from
        the \b Vertex and the \b Cube  blocks, If this is successful,
        each cube is split into simplex elements.
     .
     If the simplex generation process was successful, boundary ids are
     assigned to all boundary faces of the macrogrid.
     -# If the macrogrid was constructed in the third step of the generation
        process detailed above, then the file is first parsed for
        \b Boundarysegment  block. Here boundary ids can be
        individually assigned to each boundary segment of the macrogrid.
     -# All Boundary segments which have not
        yet been assigned an identifier and lie inside the
        interval defined by the first line of the \b Boundarydomain
        block are assigned the corresponding id. This process is then
        repeated with the remaining boundary segments using the following
        intervals defined in the \b Boundarydomain  block.
        If after this process boundary segments without id remain,
        then the default id is used if one was specified in the
        \b Boundarydomain  block - if no default was given, the
        behavior of the parser is not well defined.
     .
     \b Remark:
     -# Dune::Alberta with \c dimworld=3: \n
        if Tetgen is used to construct a
        tetrahedral grid for Dune::Alberta then the bisection routine does
        not necessarily terminate. This problem does not occur
        if the grid is constructed using the \b Interval block. If the
        simplicies are defined through the \b Simplex block then the second
        vertex of each element is taken as refinement vertex.
     -# \c dimworld=2: \n
        the refinement vertex is always chosen to be the vertex opposite the
        longest edge in each triangle.

     @subsection CONSTRCUBE Cube grids
     (Dune::ALUCubeGrid<3,3>)

     The grid is constructed using the information from the
     \b Interval  block, if present; otherwise the \b Vertex and \b Cube
     blocks are used.

     The boundary ids are assigned in the same manner as for simplex grids
     described above.

     @subsection CONSTRMIXED Mixed grids
     (\c Dune::UGGrid )

     Note that in this version only grids consisting of one element type are
     constructed even if the implemented grid allows for mixed elements.

     The verticies and elements of the grid are constructed in one of the
     following three stages:
     -# The parsers firsts looks for the \b Simplexgeneration  block.
        If this block is found, a simplex grid is generate via Triangle/Tetgen
        using as verticies the union of the verticies defined by the
        \b Interval  and the \b Vertex  block.
     -# If no \b Simplexgeneration  block is present,
        the file is parser for
        an \b Interval  block; if present a Cartesian grid is build;
        a simplex grid is generated if a \c Simplex block is present
        otherwise a cube grid is constructed.
     -# If neither a \b Simplexgeneration  nor a \b Interval
        block is found, the grid is generated using the information
        from the \b Vertex block; cube elements are constructed if the
        \b Cube block is present otherwise the \b Simplex blocks is read.
        If both a \b Cube and a \b Simplex block is found, then only
        the element information from the \b Cube block is used and each
        cube is split into simplex elements so that
        a simplex grid is constructed.
     .
     Boundary ids are assigned in the same manner as described for
     simplex grids.

   @section OPEN Still to be done
   -# There should be a mechanism to fix the desired refinement edge
      for simplex grids. In 2d an automatic selection is performed using the
      longest edge, a similar mechanism for 3d has to be implemented.
      It should be possible to turn off this automatic selection.
   -# For triangular elements it is not necessary to give the local
      numbering of the node conform to the %Dune reference element. This
      approach should be extended to tetrahedrons and to cube elements.
   -# Boundary segments for cube grids are not yet implemented.
   -# There is no possibility to single out individual boundary if
      an automatic grid generation process is used, i.e., the
      boundary segment block is ignored.

   @section EXAMPLES Examples
   <!--=========-->
   In two space dimensions:
   - \ref dgfexample1
   - \ref dgfexample2
   - \ref dgfexample3

   In three space dimensions:
   - \ref dgfexample4
   - \ref dgfexample5

   \section dgfexample1 Manual Grid Construction
   A tessellation of the unit square into six simplex entities.
   Some boundary segments on the lower and right
   boundary are given their own id the remaining are
   given a default id.

   @include examplegrid1s.dgf

   \image html  examplegrid1s.png "The resulting grid"
   \image latex examplegrid1s.eps "The resulting grid" width=\textwidth

   A tessellation into cubes using the same verticies as before,

   @include examplegrid1c.dgf

   \image html  examplegrid1c.png "The resulting grid"
   \image latex examplegrid1c.eps "The resulting grid" width=\textwidth

   Using the last input file with a simplex grid or by adding an empty
   \c Simplex block
   @code
   Simplex
 #
   @endcode
   leads to the following macro triangulation

   \image html  examplegrid1cs.png "The resulting grid"
   \image latex examplegrid1cs.eps "The resulting grid" width=\textwidth

   \section dgfexample2 Automated Grid Construction
   Automatic tessellation using Triangle,
   with verticies defined as in the example \ref dgfexample1:

   @include examplegrid1gen.dgf

   \image html  examplegrid1gen.png "The resulting grid"
   \image latex examplegrid1gen.eps "The resulting grid" width=\textwidth

   The quality of the grid can be enhanced by adding the line
   @code
   min-angle 30
   @endcode
   in the \c Simplexgenerator block

   \image html  examplegrid1genangle.png "The resulting grid"
   \image latex examplegrid1genangle.eps "The resulting grid" width=\textwidth

   Automatic tessellation using Triangle,
   with verticies are defined on a Cartesian grid with two additional
   verticies in the top right corner and one vertex outside the unit square.

   All boundary are given a default id.

   @include examplegrid2.dgf

   \image html  examplegrid2.png "The resulting grid"
   \image latex examplegrid2.eps "The resulting grid" width=\textwidth

   Adding some quality enhancement.
   The boundaries are numbered counterclockwise starting at the left
   boundary from one to four.

   @include examplegrid3.dgf

   \image html  examplegrid3.png "The resulting grid"
   \image latex examplegrid3.eps "The resulting grid" width=\textwidth

   Using both quality enhancement and a maximum area restriction.
   The bottom boundary is given the id 1 all other boundaries have
   id 2; here we do not use a default value.

   @include examplegrid4.dgf

   \image html  examplegrid4.png "The resulting grid"
   \image latex examplegrid4.eps "The resulting grid" width=\textwidth

   \section dgfexample3 Interval Domain
   A automatic tessellation of the unit square using
   a Cartesian Grid.
   All boundaries have id 1.

   @include examplegrid5.dgf

   \image html  examplegrid5c.png "The resulting grid using SGrid<2,2>"
   \image latex examplegrid5c.eps "The resulting grid using SGrid<2,2>" width=\textwidth
   \image html  examplegrid5s.png "The resulting grid using AlbertaGrid<2,2>"
   \image latex examplegrid5s.eps "The resulting grid using AlbertaGrid<2,2>" width=\textwidth

   If UGGrid<2,2> is used the result would be the same as for SGrid<2,2>.
   If an empty \c Simplex Block
   @code
   Simplex
 #
   @endcode
   is added than the same simplex grid as for AlbertaGrid<2,2> would be
   constructed.

   \section dgfexample4 Interval Domain and Automated Grid Generation
   An automatic tessellation of the unit square using
   a Cartesian Grid is shown.
   All boundaries have id 1 except boundary segment on the lower boundary
   which have value 2.

   First we use only the interval block:

   @include examplegrid6.dgf

   \image html  examplegrid6c.png "The resulting grid using ALUCubeGrid<3,3>"
   \image latex examplegrid6c.eps "The resulting grid using ALUCubeGrid<3,3>" width=\textwidth
   \image html  examplegrid6s.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex examplegrid6s.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   Now the verticies are still defined through the interval block; the simplicies
   are constructed using Tetgen (note the comment symbol \% in the
   \b Simplexgeneration block):

   @include examplegrid7.dgf

   \image html  examplegrid7.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex examplegrid7.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   Note that in the grid would be the same as in the above example if
   ALUCubeGrid<3,3> where used.

   Now we can also include some quality enhancement:

   First: \b min-angle = 1.2
          (remove the first \% in the \b Simplexgeneration block)
   \image html  examplegrid7angle.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex examplegrid7angle.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   Second: \b min-angle = 1.2 and \b max-area = 0.1
           (remove both \% in the \b Simplexgeneration block)
   \image html  examplegrid7area.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex examplegrid7area.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   \section dgfexample5 Importing Grids Generated by Tetgen/Triangle

   Here a .mesh file used to generate a 3d simplex grid through Tetgen. The
   mesh file is taken from
   http://www-c.inria.fr/Eric.Saltel/download/download.php

   @include examplegrid9.dgf

   \image html  BBEETH1M.d_cut.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex BBEETH1M.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   \image html  Orb_cut.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex Orb_cut.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   \image html  bunny.p65.param_skin.png "The resulting grid using ALUSimplexGrid<3,3>"
   \image latex bunny.p65.param_skin.eps "The resulting grid using ALUSimplexGrid<3,3>" width=\textwidth

   \section Simplexgeneration Using Tetgen/Triangle
       The freely available simplex grid generators are direcltly
       called via system
       call through the dgfparser.
       Therefore one should either add the path containing the executables of Triangle
       and/or Tetgen to the environment variable PATH or use the path
       option described below. One can either use Tetgen/Triangle directly to generate
       a .node and .ele file or one can prescribe vertices direcly in th DGF file which
       will be used to generate the grid:
       -# For the first approach use the token \b file to give a filemane
          and the type of the file (e.g. node, mesh, ply...). Example:
          \code
          file name mesh
          \endcode
          If the filetype is givn, it will be appended to \b name and this will be
          passed to Tetgen/Triangle. Additional parameters for the grid generators
          can be given by the the \b parameter token.
          If no file type is given it is assumed that in a previous run of Tetgen/Triangle
          files name.node and name.ele were generated and these will be used
          to described the verticies and elements of the Dune grid.
       -# For the second approach the vertex block and the interval blocks (if present)
          are used to generate a .node file which is passed to Tetgen/Triangle.
       .

       Some identifiers can be
       used to influence the quality of the generated mesh:
       -  \b max-area
          followed by a positive real number used as an upper bound for the
          area of all simplicies of the mesh.
       -  \b min-angle
          followed by a positive number. In 2d this limits the angles in
          resulting mesh from below; in 3d this bounds the radius-edge ratio
          from above.
       .
       If these identifiers are present then in 3d the grid is generated in two
       steps: first tetgen is called with only a .node file and then
       tetgen is run a second time including the -r switch and
       the -a and -q switch with the parameters given in the dgf file. In 2d the
       -a and -q switches are added directly in the first run.

       The remaining identifiers are
       -  The identifier \b path
          (followed by a path name) can be used
          to give search path for Triangle/Tetgen.
       -  The identifier \b display
          followed by 1 can be used to get a first impression
          of the resulting mesh using the visualization tools distributed with
          Triangle/Tetgen.
       .
       Download
       - Tetgen http://tetgen.berlios.de
       - Triangle http://www.cs.cmu.edu/~quake/triangle.html
 */


#endif
