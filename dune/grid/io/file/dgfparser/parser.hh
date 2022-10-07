// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_DUNEGRIDFORMATPARSER_HH
#define DUNE_DGF_DUNEGRIDFORMATPARSER_HH

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <dune/grid/io/file/dgfparser/entitykey.hh>

namespace Dune
{

  class DGFPrintInfo;

  /*! @brief Contains types for additional features
   */
  struct DGFBoundaryParameter
  {
    //! type of additional boundary parameters
    typedef std::string type;

    //! default constructor
    static const type &defaultValue ()
    {
      static type value;
      return value;
    }

    //! copy from string
    static type convert ( const std::string & parameter )
    {
      return parameter;
    }

    //! delimiter
    static const char delimiter = ':';
  };

  //! \brief The %DuneGridFormatParser class: reads a DGF file and stores
  //! build information in vector structures used by the MacroGrid class.
  class DuneGridFormatParser
  {
  public:
    typedef enum {Simplex,Cube,General} element_t;

    typedef enum {counterclockwise=1,clockwise=-1} orientation_t;

    //! constructor
    DuneGridFormatParser ( int rank, int size );

    /** \brief check whether a stream is in DUNE grid format
     *
     *  Actually checks whether the stream starts with the keyword 'DGF'.
     *
     *  \param  input  std::istream to check
     *
     *  \note The stream must support seeking.
     *
     *  \return whether the keyword 'DGF' was found
     */
    static bool isDuneGridFormat ( std::istream &input );

    /** \brief check whether a file is in dune grid format
     *
     *  This is just a convenience method. It calls isDuneGridFormat
     *  with a std::ifstream.
     *
     *  \param  filename  file to check
     *
     *  \return whether the keyword 'DGF' was found
     */
    static bool isDuneGridFormat ( const std::string &filename );

    /** \brief parse dune grid format from stream
     *
     *  This method actually fills the vtx, element, and bound vectors.
     *
     *  \param      input  std::istream to read the grid from
     *  \param[in]  dimG   dimension of the grid (i.e., Grid::dimension)
     *  \param[in]  dimW   dimension of the world (i.e., Grid::dimensionworld)
     *
     *  \note The stream must support seeking.
     *
     *  \returns whether reading succeeded
     */
    bool readDuneGrid( std::istream &input, int dimG, int dimW );

    //! method to write in Tetgen/Triangle Poly Format
    void writeTetgenPoly ( const std::string &, std::string &, std::string & );

    void writeTetgenPoly ( std::ostream & out, const bool writeSegments = true );

  protected:
    void generateBoundaries ( std::istream &, bool );

    // call to tetgen/triangle
    void generateSimplexGrid ( std::istream & );
    void readTetgenTriangle ( const std::string & );

    // helper methods
    void removeCopies ();

    void setOrientation ( int use1, int use2,
                          orientation_t orientation=counterclockwise );

    void setRefinement ( int use1, int use2, int is1=-1, int is2=-1 );

    double testTriang ( int snr );

    std::vector< double > & getElParam ( int i, std::vector< double > & coord );

    std::vector< double > & getVtxParam ( int i, std::vector< double > & coord );

    static std::string temporaryFileName ();

    // dimension of world and problem: set through the readDuneGrid() method
    int dimw, dimgrid;

    // vector of vertex coordinates
    std::vector < std::vector < double > > vtx;

    int nofvtx;

    int vtxoffset;

    double minVertexDistance; // min. L^1 distance of distinct points

    // vector of elements
    std :: vector< std :: vector< unsigned int > > elements;

    int nofelements;

    // vector of boundary segments + identifier
    std::vector < std::vector < int > > bound;

    int nofbound;

    // map to generate and find boundary segments
    typedef DGFBoundaryParameter::type BoundaryParameter;
    typedef std::pair < int, BoundaryParameter > BndParam;
    typedef std::map< DGFEntityKey< unsigned int >, BndParam > facemap_t;
    facemap_t facemap;

    // true if parameters on a boundary found
    bool haveBndParameters;

    // set by generator depending on element type wanted
    element_t element;

    // set by the readDuneGrid method depending
    // on what type the elements were generated
    bool simplexgrid;

    // true if grid is generated using the intervall Block
    bool cube2simplex;

    // parameters on elements
    int nofvtxparams,nofelparams;

    std::vector< std::vector< double > > vtxParams,elParams;

    // write information about generation process
    DGFPrintInfo * info;

    std::vector < double > emptyParam_;


  private:
    int rank_;
    int size_;

    template< class GridType >
    friend struct DGFGridFactory;

    template< class GridType >
    friend struct DGFBaseFactory;

  };

} // end namespace Dune

#endif
