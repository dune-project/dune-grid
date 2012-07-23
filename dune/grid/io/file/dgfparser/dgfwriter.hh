// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFWRITER_HH
#define DUNE_DGFWRITER_HH

/** \file
 *  \brief write a GridView to a DGF file
 *  \author Martin Nolte
 */

#include <fstream>
#include <vector>

#include <dune/grid/common/grid.hh>
#include <dune/geometry/referenceelements.hh>

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
    typedef typename GridView::template Codim< 0 >::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim< dimGrid >::EntityPointer VertexPointer;

    typedef typename IndexSet::IndexType Index;

  public:
    /** \brief constructor
     *
     *  \param[in]  gridView  grid view to operate on
     */
    DGFWriter ( const GridView &gridView )
      : gridView_( gridView )
    {}

    /** \brief write the GridView into a std::ostreamm
     *
     *  \param  gridout  std::ostream to write the grid to
     */
    void write ( std::ostream &gridout ) const;

    /** \brief write the GridView to a file
     *
     *  \param[in] fileName  name of the write to write the grid to
     */
    void write ( const std::string &fileName ) const;

  private:
    GridView gridView_;
  };



  template< class GV >
  inline void DGFWriter< GV >::write ( std::ostream &gridout ) const
  {
    // set the stream to full double precision
    gridout.setf( std::ios_base::scientific, std::ios_base::floatfield );
    gridout.precision( 16 );

    const IndexSet &indexSet = gridView_.indexSet();

    // write DGF header
    gridout << "DGF" << std::endl;

    const Index vxSize = indexSet.size( dimGrid );
    std::vector< Index > vertexIndex( vxSize, vxSize );

    gridout << "%" << " Elements = " << indexSet.size( 0 ) << "  |  Vertices = " << vxSize << std::endl;

    // write all vertices into the "vertex" block
    gridout << std::endl << "VERTEX" << std::endl;
    Index vertexCount = 0;
    typedef typename ElementIterator :: Entity Element ;
    const ElementIterator end = gridView_.template end< 0 >();
    for( ElementIterator it = gridView_.template begin< 0 >(); it != end; ++it )
    {
      const Element& element = *it ;
      const int numCorners = element.template count< dimGrid > ();
      for( int i=0; i<numCorners; ++i )
      {
        const Index vxIndex = indexSet.subIndex( element, i, dimGrid );
        assert( vxIndex < vxSize );
        if( vertexIndex[ vxIndex ] == vxSize )
        {
          vertexIndex[ vxIndex ] = vertexCount++;
          gridout << element.geometry().corner( i ) << std::endl;
        }
      }
    }
    gridout << "#" << std::endl;
    if( vertexCount != vxSize )
      DUNE_THROW( GridError, "Index set reports wrong number of vertices." );

    if( dimGrid > 1 )
    {
      // write all simplices to the "simplex" block
      gridout << std::endl << "SIMPLEX" << std::endl;
      for( ElementIterator it = gridView_.template begin< 0 >(); it != end; ++it )
      {
        const Element& element = *it ;
        if( !element.type().isSimplex() )
          continue;

        std::vector< Index > vertices( dimGrid+1 );
        for( size_t i = 0; i < vertices.size(); ++i )
          vertices[ i ] = vertexIndex[ indexSet.subIndex( element, i, dimGrid ) ];

        gridout << vertices[ 0 ];
        for( size_t i = 1; i < vertices.size(); ++i )
          gridout << " " << vertices[ i ];
        gridout << std::endl;
      }
      gridout << "#" << std::endl;
    }

    // write all cubes to the "cube" block
    gridout << std::endl << "CUBE" << std::endl;
    for( ElementIterator it = gridView_.template begin< 0 >(); it != end; ++it )
    {
      const Element& element = *it ;
      if( !element.type().isCube() )
        continue;

      std::vector< Index > vertices( 1 << dimGrid );
      for( size_t i = 0; i < vertices.size(); ++i )
        vertices[ i ] = vertexIndex[ indexSet.subIndex( element, i, dimGrid ) ];

      gridout << vertices[ 0 ];
      for( size_t i = 1; i < vertices.size(); ++i )
        gridout << " " << vertices[ i ];
      gridout << std::endl;
    }
    gridout << "#" << std::endl;

    // write all boundaries to the "boundarysegments" block
    gridout << std::endl << "BOUNDARYSEGMENTS" << std::endl;
    for( ElementIterator it = gridView_.template begin< 0 >(); it != end; ++it )
    {
      const Element& element = *it ;
      if( !it->hasBoundaryIntersections() )
        continue;

      const GenericReferenceElement< void, dimGrid > &refElement
        = GenericReferenceElements< void, dimGrid >::general( element.type() );

      const IntersectionIterator iend = gridView_.iend( element ) ;
      for( IntersectionIterator iit = gridView_.ibegin( element ); iit != iend; ++iit )
      {
        if( !iit->boundary() )
          continue;

        const int boundaryId = iit->boundaryId();
        if( boundaryId <= 0 )
        {
          std::cerr << "Warning: Ignoring nonpositive boundary id: "
                    << boundaryId << "." << std::endl;
          continue;
        }

        const int faceNumber = iit->indexInInside();
        const unsigned int faceSize = refElement.size( faceNumber, 1, dimGrid );
        std::vector< Index > vertices( faceSize );
        for( unsigned int i = 0; i < faceSize; ++i )
        {
          const int j = refElement.subEntity( faceNumber, 1, i, dimGrid );
          vertices[ i ] = vertexIndex[ indexSet.subIndex( element, j, dimGrid ) ];
        }
        gridout << boundaryId << "   " << vertices[ 0 ];
        for( unsigned int i = 1; i < faceSize; ++i )
          gridout << " " << vertices[ i ];
        gridout << std::endl;
      }
    }
    gridout << "#" << std::endl;

    // write the name into the "gridparameter" block
    gridout << std::endl << "GRIDPARAMETER" << std::endl;
    // gridout << "NAME " << gridView_.grid().name() << std::endl;
    gridout << "#" << std::endl;

    gridout << std::endl << "#" << std::endl;
  }


  template< class GV >
  inline void DGFWriter< GV >::write ( const std::string &fileName ) const
  {
    std::ofstream gridout( fileName.c_str() );
    write( gridout );
  }

}

#endif // #ifndef DUNE_DGFWRITER_HH
