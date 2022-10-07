// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYTHON_FUNCTION_SIMPLEGRIDFUNCTION_HH
#define DUNE_PYTHON_FUNCTION_SIMPLEGRIDFUNCTION_HH

#include <cassert>

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    // SimpleLocalFunction
    // -------------------

    template< class GridView, class LocalEvaluator >
    class SimpleLocalFunction
    {
      typedef SimpleLocalFunction< GridView, LocalEvaluator > This;

    public:
      typedef typename GridView::template Codim< 0 >::Entity Element;

      typedef typename Element::Geometry::LocalCoordinate LocalCoordinate;
      typedef std::decay_t< std::result_of_t< LocalEvaluator( Element, LocalCoordinate ) > > Value;

      explicit SimpleLocalFunction ( LocalEvaluator localEvaluator ) : localEvaluator_( std::move( localEvaluator ) ) {}
      ~SimpleLocalFunction() { unbind(); }
      template <class GF>
      explicit SimpleLocalFunction ( const GF &gf ) : SimpleLocalFunction(gf.localEvaluator()) {}

      void bind ( const Element &element ) { element_ = &element; }
      void unbind () { element_ = nullptr; }

      template< class X >
      auto operator() ( const X &x ) const
        -> decltype( std::declval< const LocalEvaluator & >()( std::declval< const Element & >(), x ) )
      {
        return localEvaluator_( element(), x );
      }

      const Element &element () const { assert( element_ ); return *element_; }

    private:
      const Element *element_ = nullptr;
      LocalEvaluator localEvaluator_;
    };



    // SimpleGridFunction
    // ------------------

    template< class GV, class LocalEval >
    class SimpleGridFunction
    {
      typedef SimpleGridFunction< GV, LocalEval > This;

    public:
      typedef LocalEval LocalEvaluator;
      static const unsigned int dimRange = LocalEvaluator::dimRange;
      typedef GV GridView;

      typedef SimpleLocalFunction< GridView, LocalEvaluator > LocalFunction;

      typedef typename LocalFunction::Element Element;
      typedef typename LocalFunction::Value Value;

      SimpleGridFunction ( const GridView &gridView, LocalEvaluator localEvaluator )
        : gridView_( gridView ), localEvaluator_( std::move( localEvaluator ) )
      {}

      const GridView &gridView () const { return gridView_; }

      const LocalEvaluator &localEvaluator () const { return localEvaluator_; }

      friend LocalFunction localFunction ( const This &gf ) { return LocalFunction( gf.localEvaluator_ ); }

    protected:
      const GridView &gridView_;
      LocalEvaluator localEvaluator_;
    };



    // LocalEvaluatorAdapter
    // ---------------------

    template< class Element, class Evaluator >
    struct LocalEvaluatorAdapter
    {
      static const unsigned int dimRange = Evaluator::dimRange;

      typedef typename Element::Geometry::GlobalCoordinate GlobalCoordinate;
      typedef typename Element::Geometry::LocalCoordinate LocalCoordinate;

      typedef decltype( std::declval< Evaluator >()( std::declval< GlobalCoordinate >() ) ) Value;

      LocalEvaluatorAdapter ( Evaluator evaluator ) : evaluator_( std::move( evaluator ) ) {}

      Value operator () ( const GlobalCoordinate &x ) const { return evaluator_( x ); }
      Value operator () ( const Element &element, const LocalCoordinate &x ) const { return evaluator_( element.geometry().global( x ) ); }

    private:
      Evaluator evaluator_;
    };



    // SimpleGlobalGridFunction
    // ------------------------

    template< class GV, class Evaluator >
    class SimpleGlobalGridFunction
      : public SimpleGridFunction< GV, LocalEvaluatorAdapter< typename GV::template Codim< 0 >::Entity, Evaluator > >
    {
      typedef SimpleGlobalGridFunction< GV, Evaluator > This;
      typedef SimpleGridFunction< GV, LocalEvaluatorAdapter< typename GV::template Codim< 0 >::Entity, Evaluator > > Base;

    public:
      typedef typename Base::GridView GridView;
      typedef typename Base::Value Value;

      typedef typename GridView::template Codim< 0 >::Geometry::GlobalCoordinate GlobalCoordinate;

      SimpleGlobalGridFunction ( const GridView &gridView, Evaluator evaluator )
        : Base( gridView, std::move( evaluator ) )
      {}

      Value operator() ( const GlobalCoordinate &x ) const { return localEvaluator_( x ); }

    protected:
      using Base::localEvaluator_;
    };



    // simpleGridFunction
    // ------------------

    template< class GV, class LE,
              class = std::result_of_t< std::decay_t< LE >( typename GV::template Codim< 0 >::Entity, typename GV::template Codim< 0 >::Geometry::LocalCoordinate ) > >
    inline static auto simpleGridFunction ( const GV &gridView, LE &&localEvaluator, PriorityTag< 1 > )
      -> SimpleGridFunction< GV, std::decay_t< LE > >
    {
      return SimpleGridFunction< GV, std::decay_t< LE > >( gridView, std::forward< LE >( localEvaluator ) );
    }

    template< class GV, class E,
              class = std::result_of_t< std::decay_t< E >( typename GV::template Codim< 0 >::Geometry::GlobalCoordinate ) > >
    inline static auto simpleGridFunction ( const GV &gridView, E &&evaluator, PriorityTag< 0 > )
      -> SimpleGlobalGridFunction< GV, std::decay_t< E > >
    {
      return SimpleGlobalGridFunction< GV, std::decay_t< E > >( gridView, std::forward< E >( evaluator ) );
    }

    template< class GridView, class Evaluator >
    inline static auto simpleGridFunction ( const GridView &gridView, Evaluator &&evaluator )
      -> decltype( simpleGridFunction( gridView, std::forward< Evaluator >( evaluator ), PriorityTag< 42 >() ) )
    {
      return simpleGridFunction( gridView, std::forward< Evaluator >( evaluator ), PriorityTag< 42 >() );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_FUNCTION_SIMPLEGRIDFUNCTION_HH
