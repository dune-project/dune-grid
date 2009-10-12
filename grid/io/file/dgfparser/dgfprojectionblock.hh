// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGF_PROJECTIONBLOCK_HH
#define DUNE_DGF_PROJECTIONBLOCK_HH

#include <map>

#include <dune/grid/common/boundaryprojection.hh>
#include <dune/grid/io/file/dgfparser/dgfparserblocks.hh>

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
          defaultKeyword, functionKeyword, sqrtKeyword, sinKeyword, cosKeyword, piKeyword,
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
      struct NumberExpression;
      struct VariableExpression;
      struct FunctionCallExpression;
      struct BracketExpression;
      struct MinusExpression;
      struct NormExpression;
      struct SqrtExpression;
      struct SinExpression;
      struct CosExpression;
      struct PowerExpression;
      struct SumExpression;
      struct DifferenceExpression;
      struct ProductExpression;
      struct QuotientExpression;

      template< int dimworld >
      struct BoundaryProjection;

    public:
      static const char *ID;

      ProjectionBlock ( std::istream &in, int dimworld );

      template< int dimworld >
      const DuneBoundaryProjection< dimworld > *defaultProjection () const
      {
        if( defaultFunction_ != 0 )
          return new BoundaryProjection< dimworld >( defaultFunction_ );
        else
          return 0;
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

      void matchToken ( const Token::Type &type, const std::string &message );
      void nextToken ();

      static char lowerCase ( char c )
      {
        return ((c >= 'A') && (c <= 'Z') ? c + ('a' - 'A') : c);
      }

    protected:
      typedef std::map< std::string, const Expression * > FunctionMap;

      using BasicBlock::line;

      Token token;
      FunctionMap functions_;
      const Expression *defaultFunction_;
    };


    std::ostream &operator<< ( std::ostream &out, const ProjectionBlock::Token &token );


    struct ProjectionBlock::Expression
    {
      typedef std::vector< double > Vector;

      virtual ~Expression ()
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const = 0;
    };


    struct ProjectionBlock::NumberExpression
      : public Expression
    {
      explicit NumberExpression ( const double value )
        : value_( value )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      double value_;
    };


    struct ProjectionBlock::VariableExpression
      : public Expression
    {
      virtual void evaluate ( const Vector &argument, Vector &result ) const;
    };


    struct ProjectionBlock::FunctionCallExpression
      : public Expression
    {
      FunctionCallExpression ( const Expression *function, const Expression *expression )
        : function_( function ),
          expression_( expression )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *function_;
      const Expression *expression_;

      mutable Vector tmp_;
    };


    struct ProjectionBlock::BracketExpression
      : public Expression
    {
      BracketExpression ( const Expression *expression, size_t field )
        : expression_( expression ),
          field_( field )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *expression_;
      size_t field_;
    };


    struct ProjectionBlock::MinusExpression
      : public Expression
    {
      explicit MinusExpression ( const Expression *expression )
        : expression_( expression )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *expression_;
    };


    struct ProjectionBlock::NormExpression
      : public Expression
    {
      explicit NormExpression ( const Expression *expression )
        : expression_( expression )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *expression_;
    };


    struct ProjectionBlock::SqrtExpression
      : public Expression
    {
      explicit SqrtExpression ( const Expression *expression )
        : expression_( expression )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *expression_;
    };


    struct ProjectionBlock::SinExpression
      : public Expression
    {
      explicit SinExpression ( const Expression *expression )
        : expression_( expression )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *expression_;
    };


    struct ProjectionBlock::CosExpression
      : public Expression
    {
      explicit CosExpression ( const Expression *expression )
        : expression_( expression )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *expression_;
    };


    struct ProjectionBlock::PowerExpression
      : public Expression
    {
      PowerExpression ( const Expression *exprA, const Expression *exprB )
        : exprA_( exprA ),
          exprB_( exprB )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *exprA_;
      const Expression *exprB_;

      mutable Vector tmp_;
    };


    struct ProjectionBlock::SumExpression
      : public Expression
    {
      explicit SumExpression ( const Expression *exprA, const Expression *exprB )
        : exprA_( exprA ),
          exprB_( exprB )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *exprA_;
      const Expression *exprB_;

      mutable Vector tmp_;
    };


    struct ProjectionBlock::DifferenceExpression
      : public Expression
    {
      explicit DifferenceExpression ( const Expression *exprA, const Expression *exprB )
        : exprA_( exprA ),
          exprB_( exprB )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *exprA_;
      const Expression *exprB_;

      mutable Vector tmp_;
    };


    struct ProjectionBlock::ProductExpression
      : public Expression
    {
      explicit ProductExpression ( const Expression *exprA, const Expression *exprB )
        : exprA_( exprA ),
          exprB_( exprB )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *exprA_;
      const Expression *exprB_;

      mutable Vector tmp_;
    };


    struct ProjectionBlock::QuotientExpression
      : public Expression
    {
      explicit QuotientExpression ( const Expression *exprA, const Expression *exprB )
        : exprA_( exprA ),
          exprB_( exprB )
      {}

      virtual void evaluate ( const Vector &argument, Vector &result ) const;

    private:
      const Expression *exprA_;
      const Expression *exprB_;
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
