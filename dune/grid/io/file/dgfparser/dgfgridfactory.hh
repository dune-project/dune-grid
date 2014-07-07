// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_GRIDFACTORY_HH
#define DUNE_DGF_GRIDFACTORY_HH

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <assert.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/dgfparser/macrogrid.hh>

#include <dune/grid/io/file/dgfparser/parser.hh>
#include <dune/grid/common/intersection.hh>


namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template < class GridImp, class IntersectionImp >
  class Intersection;



  // DGFGridFactory
  // --------------

  template < class G >
  struct DGFGridFactory
  {
    typedef G Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  private:
    typedef typename Grid::template Codim< 0 >::Entity Element;

    typedef typename Grid::template Codim< dimension >::Entity Vertex;

  public:
    Grid *grid()
    {
      return grid_;
    }

    template <class Intersection>
    bool wasInserted(const Intersection &intersection) const
    {
      return intersection.boundary();
    }

    template <class Intersection>
    int boundaryId(const Intersection &intersection) const
    {
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      if( codim == 0 )
        return macroGrid_.nofelparams;
      else if( codim == dimension )
        return macroGrid_.nofvtxparams;
      else
        return 0;
    }

    template < class Entity >
    int numParameters ( const Entity & ) const
    {
      return numParameters< Entity::codimension >();
    }

    std::vector<double>& parameter(const Element &element)
    {
      const typename Element::Geometry &geo = element.geometry();
      DomainType coord( geo.corner( 0 ) );
      for( int i = 1; i < geo.corners(); ++i )
        coord += geo.corner( i );
      coord /= double( geo.corners() );

      InsertOrderIterator it = elInsertOrder_.find( coord );
      if( it != elInsertOrder_.end() )
        return macroGrid_.elParams[ it->second ];
      assert(0);
      return emptyParam;
    }

    std::vector<double>& parameter(const Vertex &vertex)
    {
      const typename Vertex::Geometry &geo = vertex.geometry();
      DomainType coord( geo.corner( 0 ) );

      InsertOrderIterator it = vtxInsertOrder_.find( coord );
      if( it != vtxInsertOrder_.end() )
        return macroGrid_.vtxParams[ it->second ];
      return emptyParam;
    }

    // return true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      return false;
    }

    template< class GG, class II >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection< GG, II > & intersection ) const
    {
      return DGFBoundaryParameter::defaultValue();
    }

  private:
    typedef FieldVector<typename Grid::ctype,Grid::dimensionworld> DomainType;
    struct Compare
    {
      bool operator() ( const DomainType &a, const DomainType &b ) const
      {
        // returns true, if a < b; c[i] < -eps;
        const DomainType c = a - b;
        const double eps = 1e-8;

        for( int i = 0; i < DomainType::dimension; ++i )
        {
          if( c[ i ] <= -eps )
            return true;
          if( c[ i ] >= eps )
            return false;
        }
        return false;
      }
    };
    typedef std::map< DomainType, size_t, Compare > InsertOrderMap;
    typedef typename InsertOrderMap::const_iterator InsertOrderIterator;

    MacroGrid macroGrid_;
    Grid *grid_;
    InsertOrderMap elInsertOrder_;
    InsertOrderMap vtxInsertOrder_;
    std::vector<double> emptyParam;
  };

} // end namespace Dune

#endif
