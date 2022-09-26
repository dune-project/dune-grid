// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/math.hh>

#include <dune/grid/io/file/dgfparser/blocks/projection.hh>

namespace Dune
{

  namespace dgf
  {

    namespace Expr
    {

      struct ConstantExpression
        : public ProjectionBlock::Expression
      {
        explicit ConstantExpression ( const Vector &value )
          : value_( value )
        {}

        explicit ConstantExpression ( const double &value )
          : value_( 1, value )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        Vector value_;
      };


      struct VariableExpression
        : public ProjectionBlock::Expression
      {
        virtual void evaluate ( const Vector &argument, Vector &result ) const;
      };


      struct FunctionCallExpression
        : public ProjectionBlock::Expression
      {
        FunctionCallExpression ( const ProjectionBlock::ExpressionPointer& function,
                                 const ProjectionBlock::ExpressionPointer& expression )
          : function_( function ),
            expression_( expression )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer function_;
        ProjectionBlock::ExpressionPointer expression_;

        mutable Vector tmp_;
      };


      struct VectorExpression
        : public ProjectionBlock::Expression
      {
        explicit VectorExpression ( const std::vector< ProjectionBlock::ExpressionPointer >& expressions )
          : expressions_( expressions )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        std::vector< ProjectionBlock::ExpressionPointer > expressions_;
      };


      struct BracketExpression
        : public ProjectionBlock::Expression
      {
        BracketExpression ( const ProjectionBlock::ExpressionPointer& expression, size_t field )
          : expression_( expression ),
            field_( field )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer expression_;
        size_t field_;
      };


      struct MinusExpression
        : public ProjectionBlock::Expression
      {
        explicit MinusExpression ( const ProjectionBlock::ExpressionPointer& expression )
          : expression_( expression )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer expression_;
      };


      struct NormExpression
        : public ProjectionBlock::Expression
      {
        explicit NormExpression ( const ProjectionBlock::ExpressionPointer& expression )
          : expression_( expression )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer expression_;
      };


      struct SqrtExpression
        : public ProjectionBlock::Expression
      {
        explicit SqrtExpression ( const ProjectionBlock::ExpressionPointer& expression )
          : expression_( expression )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer  expression_;
      };


      struct SinExpression
        : public ProjectionBlock::Expression
      {
        explicit SinExpression ( const ProjectionBlock::ExpressionPointer& expression )
          : expression_( expression )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer expression_;
      };


      struct CosExpression
        : public ProjectionBlock::Expression
      {
        explicit CosExpression ( const ProjectionBlock::ExpressionPointer& expression )
          : expression_( expression )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer expression_;
      };


      struct PowerExpression
        : public ProjectionBlock::Expression
      {
        PowerExpression ( const ProjectionBlock::ExpressionPointer& exprA,
                          const ProjectionBlock::ExpressionPointer& exprB )
          : exprA_( exprA ),
            exprB_( exprB )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer exprA_;
        ProjectionBlock::ExpressionPointer exprB_;

        mutable Vector tmp_;
      };


      struct SumExpression
        : public ProjectionBlock::Expression
      {
        SumExpression ( const ProjectionBlock::ExpressionPointer& exprA,
                        const ProjectionBlock::ExpressionPointer& exprB )
          : exprA_( exprA ),
            exprB_( exprB )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer exprA_;
        ProjectionBlock::ExpressionPointer exprB_;

        mutable Vector tmp_;
      };


      struct DifferenceExpression
        : public ProjectionBlock::Expression
      {
        DifferenceExpression ( const ProjectionBlock::ExpressionPointer& exprA,
                               const ProjectionBlock::ExpressionPointer&  exprB )
          : exprA_( exprA ),
            exprB_( exprB )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer exprA_;
        ProjectionBlock::ExpressionPointer exprB_;

        mutable Vector tmp_;
      };


      struct ProductExpression
        : public ProjectionBlock::Expression
      {
        ProductExpression ( const ProjectionBlock::ExpressionPointer& exprA,
                            const ProjectionBlock::ExpressionPointer&  exprB )
          : exprA_( exprA ),
            exprB_( exprB )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer exprA_;
        ProjectionBlock::ExpressionPointer exprB_;

        mutable Vector tmp_;
      };


      struct QuotientExpression
        : public ProjectionBlock::Expression
      {
        QuotientExpression ( const ProjectionBlock::ExpressionPointer& exprA,
                             const ProjectionBlock::ExpressionPointer&  exprB )
          : exprA_( exprA ),
            exprB_( exprB )
        {}

        virtual void evaluate ( const Vector &argument, Vector &result ) const;

      private:
        ProjectionBlock::ExpressionPointer exprA_;
        ProjectionBlock::ExpressionPointer exprB_;
      };



      void ConstantExpression::evaluate ( const Vector & /* argument */, Vector &result ) const
      {
        result = value_;
      }


      void VariableExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        result = argument;
      }


      void FunctionCallExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        expression_->evaluate( argument, tmp_ );
        return function_->evaluate( tmp_, result );
      }


      void VectorExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        result.resize( 0 );
        Vector r;
        const auto end = expressions_.end();
        for( auto it = expressions_.begin(); it != end; ++it )
        {
          (*it)->evaluate( argument, r );
          for( size_t i = 0; i < r.size(); ++i )
            result.push_back( r[ i ] );
        }
      }


      void BracketExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        expression_->evaluate( argument, result );
        if( field_ >= result.size() )
          DUNE_THROW( MathError, "Index out of bounds (" <<  field_ << " not in [ 0, " << result.size() << " [)." );
        result[ 0 ] = result[ field_ ];
        result.resize( 1 );
      }


      void MinusExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        expression_->evaluate( argument, result );
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          result[ i ] *= -1.0;
      }


      void NormExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        using std::sqrt;
        expression_->evaluate( argument, result );
        double normsqr = 0.0;
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          normsqr += result[ i ] * result[ i ];
        result.resize( 1 );
        result[ 0 ] = sqrt( normsqr );
      }


      void SqrtExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        using std::sqrt;
        expression_->evaluate( argument, result );
        if( result.size() != 1 )
          DUNE_THROW( MathError, "Cannot calculate square root of a vector." );
        result[ 0 ] = sqrt( result[ 0 ] );
      }


      void SinExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        expression_->evaluate( argument, result );
        if( result.size() != 1 )
          DUNE_THROW( MathError, "Cannot calculate the sine of a vector." );
        result[ 0 ] = sin( result[ 0 ] );
      }


      void CosExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        expression_->evaluate( argument, result );
        if( result.size() != 1 )
          DUNE_THROW( MathError, "Cannot calculate the cosine of a vector." );
        result[ 0 ] = cos( result[ 0 ] );
      }


      void PowerExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        exprA_->evaluate( argument, result );
        exprB_->evaluate( argument, tmp_ );

        if( (result.size() == 1) && (tmp_.size() == 1) )
          result[ 0 ] = std::pow( result[ 0 ], tmp_[ 0 ] );
        else
          DUNE_THROW( MathError, "Cannot calculate powers of vectors." );
      }


      void SumExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        exprA_->evaluate( argument, result );
        exprB_->evaluate( argument, tmp_ );

        if( result.size() == tmp_.size() )
        {
          const size_t size = result.size();
          for( size_t i = 0; i < size; ++i )
            result[ i ] += tmp_[ i ];
        }
        else
          DUNE_THROW( MathError, "Cannot sum vectors of different size." );
      }


      void DifferenceExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        exprA_->evaluate( argument, result );
        exprB_->evaluate( argument, tmp_ );

        if( result.size() == tmp_.size() )
        {
          const size_t size = result.size();
          for( size_t i = 0; i < size; ++i )
            result[ i ] -= tmp_[ i ];
        }
        else
          DUNE_THROW( MathError, "Cannot sum vectors of different size." );
      }


      void ProductExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        exprA_->evaluate( argument, result );
        exprB_->evaluate( argument, tmp_ );

        if( result.size() == tmp_.size() )
        {
          double product = 0.0;
          const size_t size = result.size();
          for( size_t i = 0; i < size; ++i )
            product += result[ i ] * tmp_[ i ];
          result.resize( 1 );
          result[ 0 ] = product;
        }
        else if( tmp_.size() == 1 )
        {
          const size_t size = result.size();
          for( size_t i = 0; i < size; ++i )
            result[ i ] *= tmp_[ 0 ];
        }
        else if( result.size() == 1 )
        {
          std::swap( result, tmp_ );
          const size_t size = result.size();
          for( size_t i = 0; i < size; ++i )
            result[ i ] *= tmp_[ 0 ];
        }
        else
          DUNE_THROW( MathError, "Cannot multiply non-scalar vectors of different size." );
      }


      void QuotientExpression::evaluate ( const Vector &argument, Vector &result ) const
      {
        exprB_->evaluate( argument, result );
        if( result.size() != 1 )
          DUNE_THROW( MathError, "Cannot divide by a vector." );
        double factor = 1.0 / result[ 0 ];

        exprA_->evaluate( argument, result );
        const size_t size = result.size();
        for( size_t i = 0; i < size; ++i )
          result[ i ] *= factor;
      }

    } // namespace Expr



    // ProjectionBlock
    // ---------------

    ProjectionBlock::ProjectionBlock ( std::istream &in, int dimworld )
      : BasicBlock( in, blockId() ),
        defaultFunction_()
    {
      // for backup and load balancing
      registerProjectionFactory( dimworld );

      while( getnextline() )
      {
        std::string thisLine = line.str();
        nextToken();

        if( token.type == Token::functionKeyword )
        {
          nextToken();
          parseFunction( thisLine );
          // std::cout << "Projection line: '" << thisLine << "'" << std::endl;
        }
        else if( token.type == Token::defaultKeyword )
        {
          nextToken();
          parseDefault();
        }
        else if( token.type == Token::segmentKeyword )
        {
          nextToken();
          parseSegment();
        }
        else if( token.type != Token::endOfLine )
          DUNE_THROW( DGFException, "Error in " << *this << ": Invalid token (" << token << ")." );
        matchToken( Token::endOfLine, "trailing tokens on line." );
      }
    }


    void ProjectionBlock::parseFunction ( const std::string& exprname )
    {
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": function name expected." );
      const std::string functionName = token.literal;
      if( functions_.find( functionName ) != functions_.end() )
        DUNE_THROW( DGFException, "Error in " << *this << ": redeclaration of function " << functionName << "." );
      nextToken();

      matchToken( Token::openingParen, "'(' expected." );
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": variable name expected." );
      const std::string variableName = token.literal;
      nextToken();
      matchToken( Token::closingParen, "')' expected." );

      matchToken( Token::equals, "'=' expected." );
      ExpressionPointer expression = parseExpression( variableName );

      //std::cout << std::endl << "Declaring function: " << functionName << "( " << variableName << " )" << std::endl;
      functions_[ functionName ] = std::make_pair( expression, exprname );
    }


    ProjectionBlock::ExpressionPointer
    ProjectionBlock::parseBasicExpression ( const std::string &variableName )
    {
      ExpressionPointer expression;

      // parenthesized expression
      if( token.type == Token::openingParen )
      {
        nextToken();
        expression = parseExpression( variableName );
        matchToken( Token::closingParen, "')' expected." );
      }
      // vector constant
      else if( token.type == Token::openingBracket )
      {
        std::vector< ExpressionPointer > expressions;
        nextToken();
        while( token.type != Token::closingBracket )
        {
          expressions.push_back( parseExpression( variableName ) );
          if( (token.type != Token::comma) && (token.type != Token::closingBracket) )
          {
            std::cerr << "Warning: Components of vector expressions should be "
                      << "separated by ','." << std::endl;
            std::cerr << "         This separation will be mandatory in future "
                      << "versions." << std::endl;
          }
          if( token.type == Token::comma )
            nextToken();
        }
        nextToken();
        expression.reset( new Expr::VectorExpression( expressions ) );
      }
      // norm expression
      else if( token.type == Token::normDelim )
      {
        nextToken();
        expression.reset( new Expr::NormExpression( parseExpression( variableName ) ) );
        matchToken( Token::normDelim, "'|' expected." );
      }
      // number
      else if( token.type == Token::number )
      {
        expression.reset( new Expr::ConstantExpression( token.value ) );
        nextToken();
      }
      // pi
      else if( token.type == Token::piKeyword )
      {
        expression.reset( new Expr::ConstantExpression( MathematicalConstants< double >::pi() ) );
        nextToken();
      }
      else if( token.type == Token::string )
      {
        if( token.literal != variableName )
        {
          FunctionMap::iterator it = functions_.find( token.literal );
          if( it == functions_.end() )
            DUNE_THROW( DGFException, "Error in " << *this << ": function " << token.literal << " not declared." );
          nextToken();
          matchToken( Token::openingParen, "'(' expected." );
          expression.reset( new Expr::FunctionCallExpression( it->second.first, parseExpression( variableName ) ) );
          matchToken( Token::closingParen, "')' expected." );
        }
        else
        {
          expression.reset( new Expr::VariableExpression );
          nextToken();
        }
      }
      else
        DUNE_THROW( DGFException, "Error in " << *this << ": " << "basic expression expected." );

      return expression;
    }


    ProjectionBlock::ExpressionPointer
    ProjectionBlock::parsePostfixExpression ( const std::string &variableName )
    {
      ProjectionBlock::ExpressionPointer expression = parseBasicExpression( variableName );
      if( token.type == Token::openingBracket )
      {
        nextToken();
        if( (token.type != Token::number) || (double( int( token.value ) ) != token.value) )
          DUNE_THROW( DGFException, "Error in " << *this << ": integral number expected." );
        expression.reset( new Expr::BracketExpression( expression, int( token.value ) ) );
        nextToken();
        matchToken( Token::closingBracket, "']' expected." );
      }
      return expression;
    }


    ProjectionBlock::ExpressionPointer
    ProjectionBlock::parseUnaryExpression ( const std::string &variableName )
    {
      ProjectionBlock::ExpressionPointer expression;

      // unary minus
      if( (token.type == Token::additiveOperator) && (token.symbol == '-') )
      {
        nextToken();
        expression.reset( new Expr::MinusExpression( parsePostfixExpression( variableName ) ) );
      }
      // sqrt
      else if( token.type == Token::sqrtKeyword )
      {
        nextToken();
        expression.reset( new Expr::SqrtExpression( parseUnaryExpression( variableName ) ) );
      }
      // sin
      else if( token.type == Token::sinKeyword )
      {
        nextToken();
        expression.reset( new Expr::SinExpression( parseUnaryExpression( variableName ) ) );
      }
      // cos
      else if( token.type == Token::cosKeyword )
      {
        nextToken();
        expression.reset( new Expr::CosExpression( parseUnaryExpression( variableName ) ) );
      }
      else
        expression = parsePostfixExpression( variableName );

      return expression;
    }


    ProjectionBlock::ExpressionPointer
    ProjectionBlock::parsePowerExpression ( const std::string &variableName )
    {
      ProjectionBlock::ExpressionPointer expression = parseUnaryExpression( variableName );
      while( token.type == Token::powerOperator )
      {
        nextToken();
        expression.reset( new Expr::PowerExpression( expression, parseUnaryExpression( variableName ) ) );
      }
      return expression;
    }


    ProjectionBlock::ExpressionPointer
    ProjectionBlock::parseMultiplicativeExpression ( const std::string &variableName )
    {
      ProjectionBlock::ExpressionPointer expression = parsePowerExpression( variableName );
      while( token.type == Token::multiplicativeOperator )
      {
        const char symbol = token.symbol;
        nextToken();
        if( symbol == '*' )
          expression.reset( new Expr::ProductExpression( expression, parsePowerExpression( variableName ) ) );
        else if( symbol == '/' )
          expression.reset( new Expr::QuotientExpression( expression, parsePowerExpression( variableName ) ) );
        else
          DUNE_THROW( DGFException, "Error in " << *this << ": Internal tokenizer error." );
      }
      return expression;
    }


    ProjectionBlock::ExpressionPointer
    ProjectionBlock::parseExpression ( const std::string &variableName )
    {
      ProjectionBlock::ExpressionPointer expression = parseMultiplicativeExpression( variableName );
      while( token.type == Token::additiveOperator )
      {
        const char symbol = token.symbol;
        nextToken();
        if( symbol == '+' )
          expression.reset( new Expr::SumExpression( expression, parseMultiplicativeExpression( variableName ) ) );
        else if( symbol == '-' )
          expression.reset( new Expr::DifferenceExpression( expression, parseMultiplicativeExpression( variableName ) ) );
        else
          DUNE_THROW( DGFException, "Error in " << *this << ": Internal tokenizer error." );
      }
      return expression;
    }


    void ProjectionBlock::parseDefault ()
    {
      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": function name expected." );
      const std::string functionName = token.literal;
      nextToken();

      //std::cout << std::endl << "Default function: " << functionName << std::endl;
      FunctionMap::iterator it = functions_.find( functionName );
      if( it == functions_.end() )
        DUNE_THROW( DGFException, "Error in " << *this << ": function " << functionName << " not declared." );
      defaultFunction_ = it->second;
    }


    void ProjectionBlock::parseSegment ()
    {
      std::vector< unsigned int > faceId;
      while( token.type == Token::number )
      {
        if( double( (unsigned int)token.value ) != token.value )
          DUNE_THROW( DGFException, "Error in " << *this << ": integral number expected." );
        faceId.push_back( (unsigned int)token.value );
        nextToken();
      }

      if( token.type != Token::string )
        DUNE_THROW( DGFException, "Error in " << *this << ": function name expected." );
      const std::string functionName = token.literal;
      nextToken();

      //std::cout << std::endl << "Boundary projection for face";
      //for( size_t int i = 0; i < faceId.size(); ++i )
      //  std::cout << " " << faceId[ i ];
      //std::cout << ": " << functionName << std::endl;
      FunctionMap::iterator it = functions_.find( functionName );
      if( it == functions_.end() )
        DUNE_THROW( DGFException, "Error in " << *this << ": function " << functionName << " not declared." );
      boundaryFunctions_.push_back( std::make_pair( faceId, it->second ) );
    }


    void ProjectionBlock::matchToken ( const Token::Type &type, const std::string &message )
    {
      if( token.type != type )
        DUNE_THROW( DGFException, "Error in " << *this << ": " << message );
      if( type != Token::endOfLine )
        nextToken();
    }


    void ProjectionBlock::nextToken ()
    {
      int c;

      // eat white space
      while( ((c = line.peek()) == ' ') || (c == '\t') || (c == '\r') )
        line.get();

      // parse string literals
      if( ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z')) )
      {
        token.type = Token::string;
        token.literal = "";
        while( ((c >= 'a') && (c <= 'z')) || ((c >= 'A') && (c <= 'Z')) )
        {
          token.literal += lowerCase( line.get() );
          c = line.peek();
        }

        if( token.literal == "default" )
          token.type = Token::defaultKeyword;
        else if( token.literal == "function" )
          token.type = Token::functionKeyword;
        else if( token.literal == "segment" )
          token.type = Token::segmentKeyword;
        else if( token.literal == "sqrt" )
          token.type = Token::sqrtKeyword;
        else if( token.literal == "sin" )
          token.type = Token::sinKeyword;
        else if( token.literal == "cos" )
          token.type = Token::cosKeyword;
        else if( token.literal == "pi" )
          token.type = Token::piKeyword;
      }
      // parse numeric constant
      else if( (c >= '0') && (c <= '9') )
      {
        token.type = Token::number;
        token.value = 0;
        while( (c >= '0') && (c <= '9') )
        {
          token.value = 10*token.value + double( c - '0' );
          token.literal += char( line.get() );
          c = line.peek();
        }
        if( c == '.' )
        {
          token.literal += line.get();
          c = line.peek();
          double factor = 0.1;
          while( (c >= '0') && (c <= '9') )
          {
            token.value += factor * double( c - '0' );
            token.literal += line.get();
            factor *= 0.1;
            c = line.peek();
          }
        }
      }
      // parse single character tokens
      else if( c == ',' )
        token.setSymbol( Token::comma, line.get() );
      else if( c == '=' )
        token.setSymbol( Token::equals, line.get() );
      else if( c == '(' )
        token.setSymbol( Token::openingParen, line.get() );
      else if( c == ')' )
        token.setSymbol( Token::closingParen, line.get() );
      else if( c == '[' )
        token.setSymbol( Token::openingBracket, line.get() );
      else if( c == ']' )
        token.setSymbol( Token::closingBracket, line.get() );
      else if( c == '|' )
        token.setSymbol( Token::normDelim, line.get() );
      else if( (c == '+') || (c == '-') )
        token.setSymbol( Token::additiveOperator, line.get() );
      else if( c == '*' )
      {
        c = line.get();
        if( (line.peek() == '*') )
        {
          token.type = Token::powerOperator;
          line.get();
        }
        else
          token.setSymbol( Token::multiplicativeOperator, c );
      }
      else if( c == '/' )
        token.setSymbol( Token::multiplicativeOperator, line.get() );
      // parse end of line
      else if( c == std::stringstream::traits_type::eof() )
        token.type = Token::endOfLine;
      else
        DUNE_THROW( DGFException, "Invalid character parsed: code=0x" << std::hex << c << "." );

      //std::cout << " " << token << std::flush;
    }



    std::ostream &operator<< ( std::ostream &out, const ProjectionBlock::Token &token )
    {
      typedef ProjectionBlock::Token Token;
      switch( token.type )
      {
      case Token::string :
        return out << "string [" << token.literal << "]";
      case Token::number :
        return out << "number [" << token.value << "]";
      case Token::defaultKeyword :
        return out << "default";
      case Token::functionKeyword :
        return out << "function";
      case Token::segmentKeyword :
        return out << "segment";
      case Token::sqrtKeyword :
        return out << "sqrt";
      case Token::sinKeyword :
        return out << "sin";
      case Token::cosKeyword :
        return out << "cos";
      case Token::piKeyword :
        return out << "pi";
      case Token::equals :
        return out << "'='";
      case Token::openingParen :
        return out << "'('";
      case Token::closingParen :
        return out << "')'";
      case Token::openingBracket :
        return out << "'['";
      case Token::closingBracket :
        return out << "']'";
      case Token::normDelim :
        return out << "'|'";
      case Token::additiveOperator :
        return out << "addop [" << token.symbol << "]";
      case Token::multiplicativeOperator :
        return out << "mulop [" << token.symbol << "]";
      case Token::powerOperator :
        return out << "powerop" << std::endl;
      case Token::endOfLine :
        return out << "eol";
      default :
        return out << "invalid [" << token.type << "]";
      }
    }

  } // namespace dgf

} // namespace Dune
