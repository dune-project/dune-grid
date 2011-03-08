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

    private:
      template< int dimworld >
      struct BoundaryProjection;

    public:
      ProjectionBlock ( std::istream &in, int dimworld );

      template< int dimworld >
      const DuneBoundaryProjection< dimworld > *defaultProjection () const
      {
        if( defaultFunction_ != 0 )
          return new BoundaryProjection< dimworld >( defaultFunction_ );
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

      const Expression *function ( const std::string &name ) const
      {
        const FunctionMap::const_iterator it = functions_.find( name );
        return (it != functions_.end() ? it->second : 0);
      }

    private:
      void parseFunction ();
      const Expression *parseBasicExpression ( const std::string &variableName );
      const Expression *parsePostfixExpression ( const std::string &variableName );
      const Expression *parseUnaryExpression ( const std::string &variableName );
      const Expression *parsePowerExpression ( const std::string &variableName );
      const Expression *parseMultiplicativeExpression ( const std::string &variableName );
      const Expression *parseExpression ( const std::string &variableName );
      void parseDefault ();
      void parseSegment ();

      void matchToken ( const Token::Type &type, const std::string &message );
      void nextToken ();

      static char lowerCase ( char c )
      {
        return ((c >= 'A') && (c <= 'Z') ? c + ('a' - 'A') : c);
      }

    protected:
      typedef std::map< std::string, const Expression * > FunctionMap;
      typedef std::pair< std::vector< unsigned int >, const Expression * > BoundaryFunction;

      using BasicBlock::line;

      Token token;
      FunctionMap functions_;
      const Expression *defaultFunction_;
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

    public:
      typedef typename Base::CoordinateType CoordinateType;

      BoundaryProjection ( const Expression *expression )
        : expression_( expression )
      {}

      virtual CoordinateType operator() ( const CoordinateType &global ) const
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

    private:
      const Expression *expression_;
    };

  }

}

#endif // #ifndef DUNE_DGF_PROJECTIONBLOCK_HH
