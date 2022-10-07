// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_PROJECTIONBLOCK_HH
#define DUNE_DGF_PROJECTIONBLOCK_HH

#include <map>

#include <dune/grid/common/boundaryprojection.hh>
#include <dune/grid/io/file/dgfparser/blocks/basic.hh>

namespace Dune
{

  namespace dgf
  {

    // ProjectionBlock
    // ---------------

    class ProjectionBlock
      : public BasicBlock
    {
      struct Token
      {
        friend std::ostream &operator<< ( std::ostream &, const Token & );

        enum Type
        {
          string, number,
          defaultKeyword, functionKeyword, segmentKeyword,
          sqrtKeyword, sinKeyword, cosKeyword, piKeyword,
          comma,
          equals,
          openingParen, closingParen, openingBracket, closingBracket, normDelim,
          additiveOperator, multiplicativeOperator, powerOperator,
          endOfLine
        };

        Type type;
        char symbol;
        std::string literal;
        double value;

        void setSymbol ( const Type &t, char c )
        {
          type = t;
          symbol = c;
        }
      };

      friend std::ostream &operator<< ( std::ostream &, const Token & );

    public:
      struct Expression;

      typedef std::shared_ptr< Expression > ExpressionPointer;
      typedef std::pair< ExpressionPointer, std::string > ExpressionPair ;

    private:
      template< int dimworld >
      class BoundaryProjection;

      void registerProjectionFactory( const int dimworld );

      static const char* blockId() { return "Projection"; }
    public:
      ProjectionBlock ( std::istream &in, int dimworld );

      template< int dimworld >
      const DuneBoundaryProjection< dimworld > *defaultProjection () const
      {
        if( defaultFunction_.first )
        {
          return new BoundaryProjection< dimworld >( defaultFunction_ );
        }
        else
          return 0;
      }

      size_t numBoundaryProjections () const
      {
        return boundaryFunctions_.size();
      }

      const std::vector< unsigned int > &boundaryFace ( const size_t i ) const
      {
        assert( i < numBoundaryProjections() );
        return boundaryFunctions_[ i ].first;
      }

      template< int dimworld >
      const DuneBoundaryProjection< dimworld > *boundaryProjection ( const size_t i ) const
      {
        assert( i < numBoundaryProjections() );
        return new BoundaryProjection< dimworld >( boundaryFunctions_[ i ].second );
      }

      ExpressionPointer function ( const std::string &name ) const
      {
        const FunctionMap::const_iterator it = functions_.find( name );
        return (it != functions_.end() ? it->second.first : 0);
      }

      ExpressionPair lastFunctionInserted () const
      {
        assert( ! functions_.empty() );
        return functions_.begin()->second;
      }

      static ProjectionBlock::ExpressionPair createExpression( const std::string& funcexpr, const int dimworld )
      {
        std::stringstream str;
        str << blockId() << std::endl;
        str << funcexpr << std::endl;
        str << "#" << std::endl;
        ProjectionBlock problock( str, dimworld );
        return problock.lastFunctionInserted();
      }



    private:
      void parseFunction ( const std::string& exprname );
      ExpressionPointer parseBasicExpression ( const std::string &variableName );
      ExpressionPointer parsePostfixExpression ( const std::string &variableName );
      ExpressionPointer parseUnaryExpression ( const std::string &variableName );
      ExpressionPointer parsePowerExpression ( const std::string &variableName );
      ExpressionPointer parseMultiplicativeExpression ( const std::string &variableName );
      ExpressionPointer parseExpression ( const std::string &variableName );
      void parseDefault ();
      void parseSegment ();

      void matchToken ( const Token::Type &type, const std::string &message );
      void nextToken ();

      static char lowerCase ( char c )
      {
        return ((c >= 'A') && (c <= 'Z') ? c + ('a' - 'A') : c);
      }

    protected:
      typedef std::map< std::string, ExpressionPair >     FunctionMap;
      typedef std::pair< std::vector< unsigned int >, ExpressionPair > BoundaryFunction;

      using BasicBlock::line;

      Token token;
      FunctionMap functions_;
      ExpressionPair defaultFunction_;
      std::vector< BoundaryFunction > boundaryFunctions_;
    };


    std::ostream &operator<< ( std::ostream &out, const ProjectionBlock::Token &token );


    struct ProjectionBlock::Expression
    {
      typedef std::vector< double > Vector;

      virtual ~Expression ()
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const = 0;
    };


    template< int dimworld >
    class ProjectionBlock::BoundaryProjection
      : public DuneBoundaryProjection< dimworld >
    {
      typedef DuneBoundaryProjection< dimworld > Base;
      typedef BoundaryProjection < dimworld >    This;
      typedef typename Base :: ObjectStreamType  ObjectStreamType;

    public:
      typedef typename Base::CoordinateType CoordinateType;

      BoundaryProjection ( const ExpressionPair& exprpair )
        : expression_( exprpair.first ),
          expressionName_( exprpair.second )
      {}

      BoundaryProjection( ObjectStreamType& buffer )
      {
        int size = 0;
        buffer.read( (char *) &size, sizeof(int) );
        expressionName_.resize( size );
        buffer.read( (char *) expressionName_.c_str(), size );
        expression_ = ProjectionBlock::createExpression( expressionName_, dimworld ).first;
      }

      virtual CoordinateType operator() ( const CoordinateType &global ) const override
      {
        std::vector< double > x( dimworld );
        for( int i = 0; i < dimworld; ++i )
          x[ i ] = global[ i ];
        std::vector< double > y;
        expression_->evaluate( x, y );
        CoordinateType result;
        for( int i = 0; i < dimworld; ++i )
          result[ i ] = y[ i ];
        return result;
      }

      // backup name of expression that should allow to recreate this class
      virtual void backup( std::stringstream& buffer ) const override
      {
        buffer.write( (const char *) &key(), sizeof( int ));
        int size = expressionName_.size();
        buffer.write( (const char *) &size, sizeof(int) );
        buffer.write( expressionName_.c_str(), size );
      }

      static void registerFactory()
      {
        if( key() < 0 )
        {
          key() = Base::template registerFactory< This >();
        }
      }

    protected:
      static int& key ()
      {
        static int k = -1;
        return k;
      }

      ExpressionPointer expression_;
      std::string expressionName_;
    };

    inline void ProjectionBlock::registerProjectionFactory( const int dimworld )
    {
      if( dimworld == 3 ) {
        BoundaryProjection< 3 > :: registerFactory();
      }
      else if ( dimworld == 2 ) {
        BoundaryProjection< 2 > :: registerFactory();
      }
      else if ( dimworld == 1 ) {
        BoundaryProjection< 1 > :: registerFactory();
      }
      else {
        DUNE_THROW(NotImplemented,"ProjectionBlock::registerProjectionFactory not implemented for dimworld = " << dimworld);
      }
    }

  }

}

#endif // #ifndef DUNE_DGF_PROJECTIONBLOCK_HH
