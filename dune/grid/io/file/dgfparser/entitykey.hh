// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFEnTITYKEY_HH
#define DUNE_DGFEnTITYKEY_HH

#include <iostream>
#include <vector>

#include <dune/grid/io/file/dgfparser/dgfexception.hh>

namespace Dune
{

  // DGFEntityKey
  // ------------

  template< class A >
  struct DGFEntityKey
  {
    DGFEntityKey ( const std :: vector< A > &key, bool setOrigKey = true );
    DGFEntityKey ( const std::vector< A > &key,
                   int N, int offset, bool setOrigKey = true );
    DGFEntityKey ( const DGFEntityKey< A > &k );

    DGFEntityKey< A > &operator= ( const DGFEntityKey< A > &k );

    inline const A &operator[] ( int i ) const;
    inline bool operator < ( const DGFEntityKey< A > &k ) const;

    void orientation ( int base, std :: vector< std :: vector< double > > &vtx );
    void print( std :: ostream &out = std :: cerr ) const;

    inline bool origKeySet () const;
    inline const A &origKey ( int i ) const;
    inline int size () const;

  private:
    std :: vector< A > key_, origKey_;
    bool origKeySet_;
  };


  template< class A >
  inline const A &DGFEntityKey< A > :: operator[] ( int i ) const
  {
    return key_[ i ];
  }


  template< class A >
  inline bool DGFEntityKey< A > :: operator< ( const DGFEntityKey< A > &k ) const
  {
    // assert(k.key_.size()==key_.size());
    return key_ < k.key_;
  }


  template< class A >
  inline bool DGFEntityKey< A > :: origKeySet () const
  {
    return origKeySet_;
  }


  template< class A >
  inline const A &DGFEntityKey< A > :: origKey ( int i ) const
  {
    return origKey_[ i ];
  }


  template< class A >
  inline int DGFEntityKey< A > :: size () const
  {
    return key_.size();
  }



  // ElementFaceUtil
  // ---------------

  struct ElementFaceUtil
  {
    inline static int nofFaces ( int dim, const std::vector< unsigned int > &element );
    inline static int faceSize ( int dim, bool simpl );

    static DGFEntityKey< unsigned int >
    generateFace ( int dim, const std::vector< unsigned int > &element, int f );

  private:
    template< int dim >
    static DGFEntityKey< unsigned int >
    generateCubeFace( const std::vector< unsigned int > &element, int f );

    template< int dim >
    static DGFEntityKey< unsigned int >
    generateSimplexFace ( const std::vector< unsigned int > &element, int f );
  };


  inline int ElementFaceUtil::nofFaces ( int dim, const std::vector< unsigned int > &element )
  {
    switch( dim )
    {
    case 1 :
      return 2;
    case 2 :
      switch( element.size() )
      {
      case 3 :
        return 3;
      case 4 :
        return 4;
      default :
        return -1;
      }
    case 3 :
      switch( element.size() )
      {
      case 4 :
        return 4;
      case 8 :
        return 6;
      default :
        return -1;
      }
    default :
      return -1;
    }
  }


  inline int ElementFaceUtil::faceSize( int dim, bool simpl )
  {
    switch( dim )
    {
    case 1 :
      return 1;
    case 2 :
      return 2;
    case 3 :
      return (simpl ? 3 : 4);
    default :
      return -1;
    }
  }

} //end namespace Dune

// inlcude inline implementation
#include "entitykey_inline.hh"
#endif
