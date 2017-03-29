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

    typedef typename ElementIterator :: Entity Element ;
    typedef typename Element :: EntitySeed ElementSeed ;

    typedef typename IndexSet::IndexType Index;

    typedef ReferenceElements< typename Grid::ctype, dimGrid > RefElements;
    typedef typename RefElements::ReferenceElement RefElement;

  public:
    /** \brief constructor
     *
     *  \param[in]  gridView  grid view to operate on
     */
    DGFWriter ( const GridView &gridView )
      : gridView_( gridView )
    {}

    /** \brief write the GridView into a std::ostream
     *
     *  \param  gridout       std::ostream to write the grid to
     *  \param  newElemOrder  vector providing a new ordering for the elements in the given GridView
     *  \param  addParams     additional data to write to dgf file, such as projections
     *                        etc. (defaults to an emoty data stream)
     */
    void write ( std::ostream &gridout,
                 const std::vector< Index >& newElemOrder,
                 const std::stringstream& addParams = std::stringstream() ) const;

    /** \brief write the GridView into a std::ostream
     *
     *  \param  gridout    std::ostream to write the grid to
     */
    void write ( std::ostream &gridout ) const;

    /** \brief write the GridView into a std::ostream
     *
     *  \param  gridout    std::ostream to write the grid to
     *  \param  addParams  additional data to write to dgf file, such as projections
     *                     etc. (defaults to an emoty data stream)
     */
    void write ( std::ostream &gridout,
                 const std::stringstream& addParams ) const;

    /** \brief write the GridView to a file
     *
     *  \param[in] fileName  name of the write to write the grid to
     */
    void write ( const std::string &fileName ) const;

  protected:
    GridView gridView_;

  protected:
    /////////////////////////////////////////////
    //  helper methods
    /////////////////////////////////////////////

    // write all elements of type elementType
    void writeAllElements( const std::vector<ElementSeed>& elementSeeds,
                           const IndexSet& indexSet,
                           const GeometryType& elementType,
                           const std::vector< Index >& vertexIndex,
                           std::ostream &gridout ) const
    {
      if (!elementSeeds.empty()) {
        // perform grid traversal based on new element ordering
        for (const auto& seed : elementSeeds) {
          const Element element = gridView_.grid().entity(seed);
          writeElement(element, indexSet, elementType, vertexIndex, gridout);
        }
      }
      else {
        // perform default grid traversal
        for (const auto& element : elements(gridView_))
          writeElement(element, indexSet, elementType, vertexIndex, gridout);
      }
    }

    // write one element
    void writeElement( const Element& element,
                       const IndexSet& indexSet,
                       const GeometryType& elementType,
                       const std::vector< Index >& vertexIndex,
                       std::ostream &gridout ) const
    {
      // if element's type is not the same as the type to write the return
      if( element.type() != elementType )
        return ;

      // get vertex numbers of the element
      const size_t vxSize = element.subEntities( Element::dimension );
      std::vector<Index> vertices(vxSize);
      for( size_t i = 0; i < vxSize; ++i )
        vertices[ i ] = vertexIndex[ indexSet.subIndex( element, i, dimGrid ) ];

      gridout << vertices[ 0 ];
      for( size_t i = 1; i < vxSize; ++i )
        gridout << " " << vertices[ i ];
      gridout << std::endl;
    }
  };


  template< class GV >
  inline void DGFWriter< GV >::
  write ( std::ostream &gridout,
          const std::vector< Index >& newElemOrder,
          const std::stringstream& addParams ) const
  {
    // set the stream to full double precision
    gridout.setf( std::ios_base::scientific, std::ios_base::floatfield );
    gridout.precision( 16 );

    const IndexSet &indexSet = gridView_.indexSet();

    // vector containing entity seed (only needed if new ordering is given)
    std::vector< ElementSeed > elementSeeds;

    // if ordering was provided
    const size_t orderSize = newElemOrder.size() ;
    if( orderSize == indexSet.size( 0 ) )
    {
      const ElementIterator end = gridView_.template end< 0 >();
      ElementIterator it = gridView_.template begin< 0 >();

      if( it != end )
      {
        elementSeeds.resize( orderSize, (*it).seed() ) ;
        size_t countElements = 0 ;
        for( ; it != end; ++it, ++countElements )
        {
          const Element& element = *it ;
          assert( newElemOrder[ indexSet.index( element ) ] < orderSize );
          elementSeeds[ newElemOrder[ indexSet.index( element ) ] ] = element.seed();
        }

        // make sure that the size of the index set is equal
        // to the number of counted elements
        if( countElements != orderSize )
          DUNE_THROW(InvalidStateException,"DGFWriter::write: IndexSet not consecutive");
      }
    }

    // write DGF header
    gridout << "DGF" << std::endl;

    const Index vxSize = indexSet.size( dimGrid );
    std::vector< Index > vertexIndex( vxSize, vxSize );

    gridout << "%" << " Elements = " << indexSet.size( 0 ) << "  |  Vertices = " << vxSize << std::endl;

    // write all vertices into the "vertex" block
    gridout << std::endl << "VERTEX" << std::endl;
    Index vertexCount = 0;
    const ElementIterator end = gridView_.template end< 0 >();
    for( ElementIterator it = gridView_.template begin< 0 >(); it != end; ++it )
    {
      const Element& element = *it ;
      const int numCorners = element.subEntities( dimGrid );
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
      // only write simplex block if grid view contains simplices
      if( indexSet.size( GeometryTypes::simplex( dimGrid ) ) > 0 )
      {
        // write all simplices to the "simplex" block
        gridout << std::endl << "SIMPLEX" << std::endl;

        // write all simplex elements
        writeAllElements( elementSeeds, indexSet, GeometryTypes::simplex( dimGrid ), vertexIndex, gridout );

        // write end marker for block
        gridout << "#" << std::endl;
      }
    }

    // only write cube block if grid view contains cubes
    if( indexSet.size( GeometryTypes::cube( dimGrid ) ) > 0 )
    {
      // write all cubes to the "cube" block
      gridout << std::endl << "CUBE" << std::endl;

      // write all simplex elements
      writeAllElements( elementSeeds, indexSet, GeometryTypes::cube( dimGrid ), vertexIndex, gridout );

      // write end marker for block
      gridout << "#" << std::endl;
    }

    // add additional parameters given by the user
    gridout << addParams.str() << std::endl;

    gridout << std::endl << "#" << std::endl;
  }

  template< class GV >
  inline void DGFWriter< GV >::
  write ( std::ostream &gridout) const
  {
    // empty vector means no new ordering
    std::vector< Index > noNewOrdering ;
    std::stringstream addParams;
    write( gridout, noNewOrdering, addParams );
  }

  template< class GV >
  inline void DGFWriter< GV >::
  write ( std::ostream &gridout, const std::stringstream& addParams ) const
  {
    // empty vector means no new ordering
    std::vector< Index > noNewOrdering ;
    write( gridout, noNewOrdering, addParams );
  }

  template< class GV >
  inline void DGFWriter< GV >::write ( const std::string &fileName ) const
  {
    std::ofstream gridout( fileName.c_str() );
    if( gridout )
      write( gridout );
    else
      std::cerr << "Couldn't open file `"<< fileName << "'!"<< std::endl;
  }

}

#endif // #ifndef DUNE_DGFWRITER_HH
