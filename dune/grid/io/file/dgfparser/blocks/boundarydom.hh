// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_BOUNDARYDOMBLOCK_HH
#define DUNE_DGF_BOUNDARYDOMBLOCK_HH

#include <iostream>
#include <string>
#include <vector>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>


namespace Dune
{

  namespace dgf
  {

    struct DomainData
    {
      typedef DGFBoundaryParameter::type BoundaryParameter;

      DomainData ()
        : id_( 0 ),
          parameter_( DGFBoundaryParameter::defaultValue() ),
          defaultData_( false )
      { }

      ~DomainData ()  { }

      // constructor
      DomainData ( int id, BoundaryParameter parameter, bool defaultData = false )
        : id_( id ),
          parameter_( parameter ),
          defaultData_( defaultData )
      { }

      // return id
      int id () const
      {
        return id_;
      }

      // return true, if additional parameters given
      bool hasParameter () const
      {
        return (!parameter_.empty());
      }

      // return additional parameters
      const BoundaryParameter & parameter () const
      {
        return parameter_;
      }

      // reset data
      void reset ( int id, BoundaryParameter parameter, bool defaultData = false )
      {
        id_ = id;
        parameter_ = parameter;
        defaultData_ = defaultData;
      }

      // returns true if data origins from default boundary domain
      bool isDefault () const
      {
        return defaultData_;
      }

      friend std::ostream & operator<< ( std :: ostream & os, const DomainData & ddata )
      {
        os << "domain data: id = " << ddata.id();
        if( ddata.hasParameter() )
          os << ", parameter = " << ddata.parameter();
        return os;
      }

    private:
      int id_;
      BoundaryParameter parameter_;
      bool defaultData_;

    }; // end struct DomainData


    struct Domain
    {
      // dimension of world coordinates
      const int dimensionworld;

      typedef DGFBoundaryParameter::type BoundaryParameter;

      // constructor
      Domain( std::vector< double > p1, std::vector< double > p2, int id, BoundaryParameter & parameter )
        : dimensionworld( p1.size() ),
          left_( p1 ),
          right_( p2 ),
          data_( id, parameter )
      {
        if( int( p2.size() ) != dimensionworld )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this << "!");
        }
      }

      // constructor
      Domain( std::vector< double > p1, std::vector< double > p2, DomainData & data )
        : dimensionworld( p1.size() ),
          left_( p1 ),
          right_( p2 ),
          data_( data )
      {
        if( int( p2.size() ) != dimensionworld )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this << "!");
        }
      }

      // copy constructor
      Domain ( const Domain & other )
        : dimensionworld( other.dimensionworld ),
          left_( other.left_ ),
          right_( other.right_ ),
          data_( other.data_ )
      {
        if( dimensionworld != other.dimensionworld )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this << "!");
        }
      }

      // assignment
      Domain & operator = ( const Domain & other )
      {
        if( dimensionworld != other.dimensionworld )
        {
          DUNE_THROW(DGFException,
                     "ERROR in " << *this << "!");
        }

        left_ = other.left_;
        right_= other.right_;
        data_= other.data_;
        return *this;
      }

      // return true if point is contained in boundary domain
      template< class Vector >
      bool contains ( const Vector & x ) const
      {
        bool ret = true;
        for( int i = 0; i < dimensionworld; ++i )
        {
          if( x[ i ] < left_[ i ] || x[ i ] > right_[ i ] )
            ret = false;
        }
        return ret;
      }

      const DomainData & data () const
      {
        return data_;
      }

      // for error messages
      friend std::ostream & operator<< ( std :: ostream &os, const Domain & domain )
      {
        os << "domain: " << std::endl;
        os << "left = ";
        for( int i = 0; i < domain.dimensionworld; ++i )
          os << domain.left_[ i ] << "  ";
        os << std::endl;
        os << "right = ";
        for( int i = 0; i < domain.dimensionworld; ++i )
          os << domain.right_[ i ] << "  ";
        os << std::endl;
        os << domain.data();
        return os;
      }

    private:
      std::vector< double > left_, right_;
      DomainData data_;

    };

    class BoundaryDomBlock
      : public BasicBlock
    {
      typedef DGFBoundaryParameter::type BoundaryParameter;

      // the dimension of the vertices (is given  from user)
      int dimworld_;

      // internal counter
      int counter_;

      // default values if given
      DomainData * default_;

      // storage for all domains;
      int ndomains_;
      std::vector< Domain > domains_;

    public:
      // initialize vertex block and get first vertex
      BoundaryDomBlock ( std::istream & in, int cdimworld );

      // destructor
      ~BoundaryDomBlock ()
      {
        if( default_ )
          delete default_;
      }

      // go to next domain in block
      bool next ()
      {
        counter_++;
        return ( counter_ < ndomains_ );
      }

      // return domain
      const Domain & domain () const
      {
        return domains_.at( counter_ );
      }

      // return true if default is given
      bool hasDefaultData () const
      {
        return bool( default_ );
      }

      // return default data
      const DomainData * defaultData () const
      {
        return default_;
      }

      // return true if any boundary domain block has
      // additional parameters
      bool hasParameter () const;

      void reset ()
      {
        BasicBlock::reset();
        counter_ = -1;
      }

      // return true while block is active
      bool ok ()
      {
        return ( counter_ <= ndomains_ );
      }

      // return data if all vectors in array are contained within
      // a single domain
      template< class Vector >
      const DomainData * contains ( const std::vector< Vector > & v ) const
      {
        std::vector< int > index( ndomains_ );
        for( int i = 0; i <  ndomains_; ++i)
          index[ i ] = i;

        size_t N = v.size();
        for( size_t i = 0; i <  N; ++i )
        {
          if( index.empty() )
            break;

          const int n = index.size();
          assert( n > 0 );
          for( int j = n-1; j >= 0; --j )
          {
            bool inside = domains_[ index[ j ] ].contains( v[ i ] );
            if( !inside )
              index.erase( index.begin() + j );
          }
        }

        // check wheter no boundary domain found
        if( index.empty() )
          return default_;

        // check for ambiguity
        if( index.size() > 1 )
          dwarn << "WARNING: ambiguous boundary domain assignment, use first boundary domain in list" << std::endl;

        return &domains_[ index[ 0 ] ].data();
      }

    private:
      void readBlock ();
    };

  } // end namespace dgf

} // end namespace Dune

#endif
