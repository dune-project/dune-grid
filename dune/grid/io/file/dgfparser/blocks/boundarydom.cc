// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/boundarydom.hh>

namespace Dune
{

  namespace dgf
  {

    // BoundaryDomBlock
    // ----------------

    BoundaryDomBlock :: BoundaryDomBlock(std::istream& in,int cdimworld )
      : BasicBlock(in, "boundarydomain" ),
        dimworld(cdimworld),
        goodline(true),
        p1(cdimworld),
        p2(cdimworld),
        bndid(0),
        bndparameter( DGFBoundaryParameter::defaultValue() ),
        withdefault(false),
        defaultvalue(0),
        defaultparameter( DGFBoundaryParameter::defaultValue() )
    {
      if ( !isactive() )
        return;

      assert( cdimworld > 0 );

      if (findtoken("default"))
      {
        int x;
        if ( getnextentry( x ) )
        {
          if( x <= 0 )
          {
            DUNE_THROW(DGFException,
                       "ERROR in " << *this
                                   << "      non-positive boundary id (" << x << ") read!");
          }

          defaultvalue = x;

          // find parameter
          std::string currentline = line.str();
          std::size_t delimiter = currentline.find( DGFBoundaryParameter::delimiter );
          if( delimiter != std::string::npos )
          {
            defaultparameter =
              DGFBoundaryParameter::convert( currentline.substr( delimiter+1, std::string::npos ) );
          }

          withdefault = true;
        }
      }

      reset();
      next();
    }


    bool BoundaryDomBlock :: next ()
    {
      assert( ok() );
      getnextline();

      if ( linenumber()==noflines() )
      {
        goodline=false;
        return goodline;
      }

      int id;
      bndparameter = DGFBoundaryParameter::defaultValue();

      if( getnextentry( id ) )
      {
        if( id <= 0 )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      non-positive boundary id (" << id << ") read!");
        }
        bndid = id;

        // find delimiter
        std::string currentline = line.str();
        std::size_t delimiter = currentline.find( DGFBoundaryParameter::delimiter );
        if( delimiter != std::string::npos )
        {
          bndparameter =
            DGFBoundaryParameter::convert( currentline.substr( delimiter+1, std::string::npos ) );
        }

        // read vertices
        double x;
        int n = 0;
        while ( getnextentry( x ) )
        {
          if ( 0 <= n && n < dimworld )
            p1.at( n ) = x;
          else if ( dimworld <= n && n < 2*dimworld )
          {
            p2.at( n-dimworld ) = x;
            if ( p2.at( n-dimworld )< p1.at( n-dimworld ) )
            {
              DUNE_THROW(DGFException,
                         "ERROR in " << *this
                                     << "      second coordinate smaller than first coordinate: "
                                     << p2.at(n-dimworld)
                                     << " read but expected value larger or equal to "
                                     << p1.at(n-dimworld));
            }
          }
          n++;
        }

        goodline = ( n == dimworld*2 );
        if ( !goodline )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this
                                 << "      wrong number of coordinates: "
                                 << n << " read but expected " << dimworld);
        }
      }
      else
        next();
      return goodline;
    }


    bool BoundaryDomBlock :: inside ( const std :: vector< double > &v ) const
    {
      assert( v.size() == (size_t)dimworld );
      for( int i = 0; i < dimworld; ++i )
      {
        if( (v[ i ] < p1[ i ]) || (v[ i ] > p2[ i ]) )
          return false;
      }
      return true;
    }

  } // end namespace dgf

} // end namespace Dune
