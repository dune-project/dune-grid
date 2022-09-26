// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ENTITYKEY_INLINE_HH
#define DUNE_ENTITYKEY_INLINE_HH

#include <algorithm>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>

namespace Dune
{

  // DGFEntityKey
  // ------------

  template< class A >
  inline DGFEntityKey< A > :: DGFEntityKey ( const std :: vector< A > &key, bool setOrigKey )
    : key_( key.size() ),
      origKey_( key.size() ),
      origKeySet_( setOrigKey )
  {
    for (size_t i=0; i<key_.size(); i++)
    {
      key_[i]=key[i];
      origKey_[i]=key_[i];
    }
    std :: sort( key_.begin(), key_.end() );
  }


  template< class A >
  inline DGFEntityKey< A > :: DGFEntityKey ( const std :: vector< A > &key,
                                             int N, int offset, bool setOrigKey )
    : key_( N ),
      origKey_( N ),
      origKeySet_( setOrigKey )
  {
    for (size_t i=0; i<key_.size(); i++)
    {
      key_[i]=key[(i+offset)%key.size()];
      origKey_[i]=key[(i+offset)%key.size()];
    }
    std :: sort( key_.begin(), key_.end() );
  }


  template< class A >
  inline DGFEntityKey< A > :: DGFEntityKey ( const DGFEntityKey< A > &k )
    : key_( k.key_.size() ),
      origKey_( k.key_.size() ),
      origKeySet_( k. origKeySet_ )
  {
    for (size_t i=0; i<key_.size(); i++)
    {
      key_[i]=k.key_[i];
      origKey_[i]=k.origKey_[i];
    }
  }


  template< class A >
  inline DGFEntityKey< A > &DGFEntityKey< A > :: operator= ( const DGFEntityKey< A > &k )
  {
    assert(key_.size()==k.key_.size());
    for (size_t i=0; i<key_.size(); i++) {
      key_[i]=k.key_[i];
      origKey_[i]=k.origKey_[i];
    }
    origKeySet_ = k.origKeySet_;
    return *this;
  }


  template< class A >
  inline void DGFEntityKey< A >
  :: orientation ( int base, std::vector< std :: vector< double > > &vtx )
  {
    if (key_.size()==3)  {
      assert( (size_t) origKey_[0] < vtx.size() );
      std::vector<double>& p0 = vtx[origKey_[0]];
      assert( (size_t) origKey_[1] < vtx.size() );
      std::vector<double>& p1 = vtx[origKey_[1]];
      assert( (size_t) origKey_[2] < vtx.size() );
      std::vector<double>& p2 = vtx[origKey_[2]];
      assert( (size_t) base < vtx.size() );
      std::vector<double>& q  = vtx[base];
      double n[3];
      n[0] = (p1[1]-p0[1])*(p2[2]-p0[2])-(p2[1]-p0[1])*(p1[2]-p0[2]);
      n[1] = (p1[2]-p0[2])*(p2[0]-p0[0])-(p2[2]-p0[2])*(p1[0]-p0[0]);
      n[2] = (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]);
      double test = n[0]*(q[0]-p0[0])+n[1]*(q[1]-p0[1])+n[2]*(q[2]-p0[2]);
      bool reorient = (test>0);
      if (reorient) {
        A key1=origKey_[1];
        origKey_[1]=origKey_[2];
        origKey_[2]=key1;
      }
    }
  }


  template< class A >
  inline void DGFEntityKey< A > :: print ( std :: ostream &out ) const
  {
    for( size_t i = 0; i < key_.size(); ++i )
      out << key_[ i ] << " ";
    out << std :: endl;
  }


  // ElementFaceUtil
  // ---------------

  template< int dim >
  inline DGFEntityKey< unsigned int >
  ElementFaceUtil::generateCubeFace
    ( const std::vector< unsigned int > &element, int f )
  {
    auto refCube = ReferenceElements< double, dim >::cube();
    const unsigned int size = refCube.size( f, 1, dim );
    std::vector< unsigned int > k( size );
    for( unsigned int i = 0; i < size; ++ i )
      k[ i ] = element[ refCube.subEntity( f, 1, i, dim ) ];
    return DGFEntityKey< unsigned int >( k );
  }


  template< int dim >
  inline DGFEntityKey< unsigned int >
  ElementFaceUtil :: generateSimplexFace
    ( const std :: vector< unsigned int > &element, int f )
  {
    auto refSimplex = ReferenceElements< double, dim >::simplex();
    const unsigned int size = refSimplex.size( f, 1, dim );
    std :: vector< unsigned int > k( size );
    for( unsigned int i = 0; i < size; ++i )
      k[ i ] = element[ refSimplex.subEntity( f, 1, i, dim ) ];
    return DGFEntityKey< unsigned int >( k );
  }


  inline DGFEntityKey< unsigned int >
  ElementFaceUtil::generateFace ( int dim, const std::vector< unsigned int > &element, int f )
  {
    if( element.size() == size_t(dim+1) )
    {
      // Simplex element
      switch( dim )
      {
      case 3 :
        return generateSimplexFace< 3 >( element, f );
      case 2 :
        return generateSimplexFace< 2 >( element, f );
      case 1 :
        return generateSimplexFace< 1 >( element, f );
      default :
        DUNE_THROW( NotImplemented, "ElementUtil::generateFace not implemented for dim = " << dim << "." );
      }
    }
    else
    {
      // Cube element
      switch( dim )
      {
      case 3 :
        return generateCubeFace< 3 >( element, f );
      case 2 :
        return generateCubeFace< 2 >( element, f );
      case 1 :
        return generateCubeFace< 1 >( element, f );
      default :
        DUNE_THROW( NotImplemented, "ElementUtil::generateFace not implemented for dim = " << dim << "." );
      }
    }
  }

} // end namespace Dune

#endif // DUNE_ENTITYKEY_INLINE_HH
