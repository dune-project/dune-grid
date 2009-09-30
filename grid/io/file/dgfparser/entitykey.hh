// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFEnTITYKEY_HH
#define DUNE_DGFEnTITYKEY_HH

#include <iostream>
#include <vector>

#include <dune/grid/alugrid/3d/topology.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

namespace Dune
{

  // DGFEntityKey
  // ------------

  template< class A > class DGFEntityKey
  {
    std :: vector< A > key_, origKey_;
    bool origKeySet_;

  public:
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

  class ElementFaceUtil
  {
  public:
    inline static int nofFaces ( int dimw, std :: vector< unsigned int > &element );
    inline static int faceSize ( int dimw, bool simpl );

    static DGFEntityKey< unsigned int >
    generateFace ( int dimw, const std :: vector< unsigned int > &element, int f );

  private:
    template< int dimworld >
    static DGFEntityKey< unsigned int >
    generateCubeFace( const std :: vector< unsigned int > &element, int f );

    template< int dimworld >
    static DGFEntityKey< unsigned int >
    generateSimplexFace ( const std :: vector< unsigned int > &element, int f );
  };


  inline int ElementFaceUtil
  :: nofFaces( int dimw, std :: vector< unsigned int > &element )
  {
    if (dimw==1)
      return 2;
    else if (dimw==2)
      switch (element.size()) {
      case 3 : return 3; break;
      case 4 : return 4; break;
      }
    else if (dimw==3)
      switch (element.size()) {
      case 4 : return 4; break;
      case 8 : return 6; break;
      }
    return -1;
  }


  inline int ElementFaceUtil :: faceSize( int dimw, bool simpl )
  {
    if (dimw==1)
      return 1;
    else if (dimw==2)
      return 2;
    else if (dimw==3)
      return ((simpl) ? 3 : 4);
    return -1;
  }

} //end namespace Dune

// inlcude inline implementation
#include "entitykey_inline.hh"
#endif
