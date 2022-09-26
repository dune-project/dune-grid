// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_IO_FILE_DGFPARSER_DGFWRITER_HH
#define DUNE_GRID_IO_FILE_DGFPARSER_DGFWRITER_HH

/** \file
 *  \brief write a GridView to a DGF file
 *  \author Martin Nolte
 */

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <fstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/rangeutilities.hh>
#include <dune/common/typeutilities.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  /** \class DGFWriter
   *  \ingroup DuneGridFormatParser
   *  \brief write a GridView to a DGF file
   *
   *  The DGFWriter allows create a DGF file from a given GridView. It allows
   *  for the easy creation of file format converters.
   *
   *  \tparam  GV  GridView to write in DGF format
   */
  template< class GV >
  class DGFWriter
  {
    typedef DGFWriter< GV > This;

  public:
    /** \brief type of grid view */
    typedef GV GridView;
    /** \brief type of underlying hierarchical grid */
    typedef typename GridView::Grid Grid;

    /** \brief dimension of the grid */
    static const int dimGrid = GridView::dimension;

  private:
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim< 0 >::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename Element::EntitySeed ElementSeed;

    typedef typename IndexSet::IndexType Index;

  public:
    /** \brief constructor
     *
     *  \param[in]  gridView  grid view to operate on
     */
    DGFWriter ( const GridView &gridView )
      : gridView_( gridView )
    {}

    /**
     * \brief write the GridView into a std::ostream
     *
     * \param      gridout       std::ostream to write the grid to
     * \param[in]  newElemOrder  vector providing a new ordering for the elements in the given GridView
     * \param[in]  boundaryData  callable attaching boundary data to each intersection
     * \param[in]  addParams     additional data to write to dgf file, such as projections etc.
     *                           (defaults to an emoty data stream)
     **/
    template< class BoundaryData >
    void write ( std::ostream &gridout, const std::vector< Index > &newElemOrder, BoundaryData &&boundaryData, const std::stringstream &addParams = std::stringstream() ) const;

    /**
     * \brief write the GridView to a file
     *
     * \param      gridout       std::ostream to write the grid to
     * \param[in]  boundaryData  callable attaching boundary data to each intersection
     * \param[in]  addParams     additional data to write to dgf file, such as projections, etc.
     *                           (defaults to an empty data stream)
     **/
    template< class BoundaryData >
    void write ( std::ostream &gridout, BoundaryData &&boundaryData, const std::stringstream &addParams = std::stringstream() ) const;

    /**
     * \brief write the GridView into a std::ostream
     *
     * \param      gridout       std::ostream to write the grid to
     * \param[in]  newElemOrder  vector providing a new ordering for the elements in the given GridView
     * \param[in]  addParams     additional data to write to dgf file, such as projections etc.
     *                           (defaults to an emoty data stream)
     **/
    void write ( std::ostream &gridout, const std::vector< Index > &newElemOrder, const std::stringstream &addParams = std::stringstream() ) const
    {
      write( gridout, newElemOrder, [] ( const Intersection &i ) -> int { return boundaryId( i ); }, addParams );
    }

    /**
     * \brief write the GridView into a std::ostream
     *
     * \param      gridout       std::ostream to write the grid to
     * \param[in]  addParams     additional data to write to dgf file, such as projections, etc.
     *                           (defaults to an empty data stream)
     **/
    void write ( std::ostream &gridout, const std::stringstream &addParams = std::stringstream() ) const
    {
      write( gridout, [] ( const Intersection &i ) -> int { return boundaryId( i ); }, addParams );
    }

    /**
     * \brief write the GridView to a file
     *
     * \param[in]  fileName  name of the write to write the grid to
     * \param[in]  args      arguments for the write method with istream
     **/
    template< class... Args >
    auto write ( const std::string &fileName, Args &&... args ) const
      -> std::void_t< decltype( this->write( std::declval< std::ostream & >(), std::declval< Args >()... ) ) >
    {
      std::ofstream gridout( fileName );
      if( gridout )
        write( gridout, std::forward< Args >( args )... );
      else
        std::cerr << "Couldn't open file `"<< fileName << "'!"<< std::endl;
    }

  protected:
    auto elementsSeeds ( const std::vector< Index > &newElemOrder ) const
      -> std::vector< ElementSeed >;

    void writeHeader ( std::ostream &gridout ) const;
    void writeFooter ( std::ostream &gridout ) const;

    auto writeVertices ( std::ostream &gridout ) const
      -> std::vector< Index >;

    void writeElement ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const Element &element, const GeometryType &elementType ) const;

    void writeSimplices ( std::ostream &gridout, const std::vector< Index > &dgfIndices ) const;
    void writeSimplices ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const std::vector< ElementSeed > &elementSeeds ) const;

    void writeCubes ( std::ostream &gridout, const std::vector< Index > &dgfIndices ) const;
    void writeCubes ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const std::vector< ElementSeed > &elementSeeds ) const;

    template< class... Args >
    void writeElements ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const Args &... args ) const;

  private:
    template< class I >
    static auto boundaryId ( const I &i, PriorityTag< 1 > )
      -> std::enable_if_t< std::is_convertible< std::decay_t< decltype( i.impl().boundaryId() ) >, int >::value, int >
    {
      return i.impl().boundaryId();
    }

    template< class I >
    static int boundaryId ( const I &i, PriorityTag< 0 > )
    {
      return 1;
    }

  protected:
    static int boundaryId ( const Intersection &i ) { return boundaryId( i, PriorityTag< 42 >() ); }

  private:
    static int boundaryId ( const Intersection &, int bndId ) { return bndId; }
    static int boundaryId ( const Intersection &i, const std::string & ) { return boundaryId( i ); }
    static int boundaryId ( const Intersection &i, const std::pair< int, std::string > &data ) { return boundrayId( i, data.first ); }

    static void appendBoundaryData ( std::ostream &gridout, int ) { gridout << std::endl; }
    static void appendBoundaryData ( std::ostream &gridout, std::pair< int, std::string > &data ) { appendBoundaryData( gridout, data.second ); }
    static void appendBoundaryData ( std::ostream &gridout, const std::string &s ) { gridout << " : " << s << std::endl; }

  protected:
    template< class BoundaryData >
    void writeBoundaries ( std::ostream &gridout, const std::vector< Index > &dgfIndices, BoundaryData &&boundaryData ) const;

    void writeBoundaries ( std::ostream &gridout, const std::vector< Index > &dgfIndices ) const
    {
      writeBoundaries( gridout, dgfIndices, [] ( const Intersection &i ) -> int { return boundaryId( i ); } );
    }

  protected:
    GridView gridView_;
  };


  template< class GV >
  inline auto DGFWriter< GV >::elementsSeeds ( const std::vector< Index > &newElemOrder ) const
    -> std::vector< ElementSeed >
  {
    const IndexSet &indexSet = gridView_.indexSet();

    const std::size_t orderSize = newElemOrder.size() ;
    std::vector< ElementSeed > elementSeeds( orderSize );

    for( const Element &element : elements( gridView_ ) )
    {
      assert( newElemOrder[ indexSet.index( element ) ] < orderSize );
      elementSeeds[ newElemOrder[ indexSet.index( element ) ] ] = element.seed();
    }

    return elementSeeds;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeHeader ( std::ostream &gridout ) const
  {
    // set the stream to full double precision
    gridout.setf( std::ios_base::scientific, std::ios_base::floatfield );
    gridout.precision( 16 );

    const IndexSet &indexSet = gridView_.indexSet();

    // write DGF header
    gridout << "DGF" << std::endl;
    gridout << "%" << " Elements = " << indexSet.size( 0 ) << "  |  Vertices = " << indexSet.size( dimGrid ) << std::endl;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeFooter ( std::ostream &gridout ) const
  {
    gridout << std::endl << "#" << std::endl;
  }


  template< class GV >
  inline auto DGFWriter< GV >::writeVertices ( std::ostream &gridout ) const
    -> std::vector< Index >
  {
    const IndexSet &indexSet = gridView_.indexSet();

    const Index vxSize = indexSet.size( dimGrid );
    std::vector< Index > dgfIndices( vxSize, vxSize );

    // write all vertices into the "vertex" block
    gridout << std::endl << "VERTEX" << std::endl;
    Index vertexCount = 0;
    for( const Element &element : elements( gridView_ ) )
    {
      for( auto i : range( element.subEntities( dimGrid ) ) )
      {
        const Index vxIndex = indexSet.subIndex( element, i, dimGrid );
        assert( vxIndex < vxSize );
        if( dgfIndices[ vxIndex ] == vxSize )
        {
          dgfIndices[ vxIndex ] = vertexCount++;
          gridout << element.geometry().corner( i ) << std::endl;
        }
      }
    }
    gridout << "#" << std::endl;

    if( vertexCount != vxSize )
      DUNE_THROW( GridError, "IndexSet reports wrong number of vertices." );
    return dgfIndices;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeElement ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const Element &element, const GeometryType &elementType ) const
  {
    // if element's type is not the same as the type to write the return
    if( element.type() != elementType )
      return;

    // write vertex numbers of the element
    const IndexSet &indexSet = gridView_.indexSet();
    for( auto i : range( element.subEntities( Element::dimension ) ) )
      gridout << (i > 0 ? " " : "") << dgfIndices[ indexSet.subIndex( element, i, dimGrid ) ];
    gridout << std::endl;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeSimplices ( std::ostream &gridout, const std::vector< Index > &dgfIndices ) const
  {
    // write all simplices to the "simplex" block
    gridout << std::endl << "SIMPLEX" << std::endl;

    // write all simplex elements
    for( const Element &element : elements( gridView_ ) )
      writeElement( gridout, dgfIndices, element, GeometryTypes::simplex( dimGrid ) );

    // write end marker for block
    gridout << "#" << std::endl;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeSimplices ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const std::vector< ElementSeed > &elementSeeds ) const
  {
    // write all simplices to the "simplex" block
    gridout << std::endl << "SIMPLEX" << std::endl;

    // write all simplex elements
    for( const ElementSeed &seed : elementSeeds )
      writeElement( gridout, dgfIndices, gridView_.grid().entity( seed ), GeometryTypes::simplex( dimGrid ) );

    // write end marker for block
    gridout << "#" << std::endl;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeCubes ( std::ostream &gridout, const std::vector< Index > &dgfIndices ) const
  {
    // write all cubes to the "cube" block
    gridout << std::endl << "CUBE" << std::endl;

    // write all cube elements
    for( const Element &element : elements( gridView_ ) )
      writeElement( gridout, dgfIndices, element, GeometryTypes::cube( dimGrid ) );

    // write end marker for block
    gridout << "#" << std::endl;
  }


  template< class GV >
  inline void DGFWriter< GV >::writeCubes ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const std::vector< ElementSeed > &elementSeeds ) const
  {
    const IndexSet &indexSet = gridView_.indexSet();

    // write all cubes to the "cube" block
    gridout << std::endl << "CUBE" << std::endl;

    // write all cube elements
    for( const ElementSeed &seed : elementSeeds )
      writeElement( gridout, dgfIndices, gridView_.grid().entity( seed ), GeometryTypes::cube( dimGrid ) );

    // write end marker for block
    gridout << "#" << std::endl;
  }


  template< class GV >
  template< class... Args >
  inline void DGFWriter< GV >::writeElements ( std::ostream &gridout, const std::vector< Index > &dgfIndices, const Args &... args ) const
  {
    const IndexSet &indexSet = gridView_.indexSet();

    if( (dimGrid > 1) && (indexSet.size( GeometryTypes::simplex( dimGrid ) ) > 0) )
      writeSimplices( gridout, dgfIndices, args... );

    if( indexSet.size( GeometryTypes::cube( dimGrid ) ) > 0 )
      writeCubes( gridout, dgfIndices, args... );
  }


  template< class GV >
  template< class BoundaryData >
  inline void DGFWriter< GV >::writeBoundaries ( std::ostream &gridout, const std::vector< Index > &dgfIndices, BoundaryData &&boundaryData ) const
  {
    using std::max;

    const IndexSet &indexSet = gridView_.indexSet();

    // write all boundaries to the "boundarysegments" block
    gridout << std::endl << "BOUNDARYSEGMENTS" << std::endl;

    for( const Element &element : elements( gridView_ ) )
    {
      if( !element.hasBoundaryIntersections() )
        continue;

      const auto &refElement = ReferenceElements< typename Grid::ctype, dimGrid >::general( element.type() );
      for( const Intersection &intersection : intersections( gridView_, element ) )
      {
        if( !intersection.boundary() )
          continue;

        const auto data = boundaryData( intersection );
        const int bndId = max( boundaryId( intersection, data ), 1 );

        const int faceNumber = intersection.indexInInside();
        const unsigned int faceSize = refElement.size( faceNumber, 1, dimGrid );
        gridout << bndId << "  ";
        for( auto i : range( faceSize ) )
        {
          const int j = refElement.subEntity( faceNumber, 1, i, dimGrid );
          gridout << " " << dgfIndices[ indexSet.subIndex( element, j, dimGrid ) ];
        }
        appendBoundaryData( gridout, data );
      }
    }
    gridout << "#" << std::endl;
  }


  template< class GV >
  template< class BoundaryData >
  inline void DGFWriter< GV >::write ( std::ostream &gridout, const std::vector< Index > &newElemOrder, BoundaryData &&boundaryData, const std::stringstream &addParams ) const
  {
    writeHeader( gridout );
    auto dgfIndices = writeVertices( gridout );
    writeElements( gridout, dgfIndices, elementSeeds( newElemOrder ) );
    writeBoundaries( gridout, dgfIndices, std::forward< BoundaryData >( boundaryData ) );
    gridout << addParams.str();
    writeFooter( gridout );
  }


  template< class GV >
  template< class BoundaryData >
  inline void DGFWriter< GV >::write ( std::ostream &gridout, BoundaryData &&boundaryData, const std::stringstream &addParams ) const
  {
    writeHeader( gridout );
    auto dgfIndices = writeVertices( gridout );
    writeElements( gridout, dgfIndices );
    writeBoundaries( gridout, dgfIndices, std::forward< BoundaryData >( boundaryData ) );
    gridout << addParams.str();
    writeFooter( gridout );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_IO_FILE_DGFPARSER_DGFWRITER_HH
