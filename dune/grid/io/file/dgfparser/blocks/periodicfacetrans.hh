// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_PERIODICFACETRANSBLOCK_HH
#define DUNE_DGF_PERIODICFACETRANSBLOCK_HH

#include <iostream>
#include <vector>

#include <dune/grid/io/file/dgfparser/blocks/basic.hh>


namespace Dune
{

  namespace dgf
  {

    // PeriodicFaceTransformationBlock
    // -------------------------------

    struct PeriodicFaceTransformationBlock
      : public BasicBlock
    {
      template< class T >
      class Matrix;

      struct AffineTransformation;

    private:
      std::vector< AffineTransformation > transformations_;

      // copy not implemented
      PeriodicFaceTransformationBlock ( const PeriodicFaceTransformationBlock & );

    public:
      // initialize block and get dimension of world
      PeriodicFaceTransformationBlock ( std::istream &in, int dimworld );

      const AffineTransformation &transformation ( int i ) const
      {
        assert( i < numTransformations() );
        return transformations_[ i ];
      }

      int numTransformations () const
      {
        return transformations_.size();
      }

    private:
      void match ( char what );
    };


    // PeriodicFaceTransformationBlock::Matrix
    // ---------------------------------------

    template< class T >
    class PeriodicFaceTransformationBlock::Matrix
    {
      int rows_;
      int cols_;
      std::vector< T > fields_;

    public:
      Matrix ( int rows, int cols )
        : rows_( rows ),
          cols_( cols ),
          fields_( rows * cols )
      {}

      const T &operator() ( int i, int j ) const
      {
        return fields_[ i * cols_ + j ];
      }

      T &operator() ( int i, int j )
      {
        return fields_[ i * cols_ + j ];
      }

      int rows () const
      {
        return rows_;
      }

      int cols () const
      {
        return cols_;
      }
    };


    // PeriodicFaceTransformationBlock::AffineTransformation
    // -----------------------------------------------------

    struct PeriodicFaceTransformationBlock::AffineTransformation
    {
      Matrix< double > matrix;
      std::vector< double > shift;

      explicit AffineTransformation ( int dimworld )
        : matrix( dimworld, dimworld ),
          shift( dimworld )
      {}
    };


    inline std::ostream &
    operator<< ( std::ostream &out, const PeriodicFaceTransformationBlock::AffineTransformation &trafo )
    {
      for( int i = 0; i < trafo.matrix.rows(); ++i )
      {
        out << (i > 0 ? ", " : "");
        for( int j = 0; j < trafo.matrix.cols(); ++j )
          out << (j > 0 ? " " : "") << trafo.matrix( i, j );
      }
      out << " +";
      for( unsigned int i = 0; i < trafo.shift.size(); ++i )
        out << " " << trafo.shift[ i ];
      return out;
    }

  } // end namespace dgf

} // end namespace Dune

#endif
