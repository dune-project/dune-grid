// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MISC_HH
#define DUNE_GENERICGEOMETRY_MISC_HH

#include <dune/common/static_assert.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< bool condition, template< bool > class True, template< bool > class False >
    struct ProtectedIf;

    template< template< bool > class True, template< bool > class False >
    struct ProtectedIf< true, True, False >
      : public True< true >
    {};

    template< template< bool > class True, template< bool > class False >
    struct ProtectedIf< false, True, False >
      : public False< false >
    {};



    template< template< int > class Operation, int first, int last >
    struct ForLoop
    {
      static void apply ()
      {
        Operation< first >::apply();
        ForLoop< Operation, first+1, last >::apply();
      }

      template< class T1 >
      static void apply ( T1 &p1 )
      {
        Operation< first >::apply( p1 );
        ForLoop< Operation, first+1, last >::apply( p1 );
      }

      template< class T1, class T2 >
      static void apply ( T1 &p1, T2 &p2 )
      {
        Operation< first >::apply( p1, p2 );
        ForLoop< Operation, first+1, last >::apply( p1, p2 );
      }

      template< class T1, class T2, class T3 >
      static void apply ( T1 &p1, T2 &p2, T3 &p3 )
      {
        Operation< first >::apply( p1, p2, p3 );
        ForLoop< Operation, first+1, last >::apply( p1, p2, p3 );
      }

      template< class T1, class T2, class T3, class T4 >
      static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4 )
      {
        Operation< first >::apply( p1, p2, p3, p4 );
        ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4 );
      }

      template< class T1, class T2, class T3, class T4, class T5 >
      static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4, T5 &p5 )
      {
        Operation< first >::apply( p1, p2, p3, p4, p5 );
        ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4, p5 );
      }

    private:
      dune_static_assert( (first <= last), "ForLoop: first > last" );
    };

    template< template< int > class Operation, int last >
    struct ForLoop< Operation, last, last >
    {
      static void apply ()
      {
        Operation< last >::apply();
      }

      template< class T1 >
      static void apply ( T1 &p1 )
      {
        Operation< last >::apply( p1 );
      }

      template< class T1, class T2 >
      static void apply ( T1 &p1, T2 &p2 )
      {
        Operation< last >::apply( p1, p2 );
      }

      template< class T1, class T2, class T3 >
      static void apply ( T1 &p1, T2 &p2, T3 &p3 )
      {
        Operation< last >::apply( p1, p2, p3 );
      }

      template< class T1, class T2, class T3, class T4 >
      static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4 )
      {
        Operation< last >::apply( p1, p2, p3, p4 );
      }

      template< class T1, class T2, class T3, class T4, class T5 >
      static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4, T5 &p5 )
      {
        Operation< last >::apply( p1, p2, p3, p4, p5 );
      }
    };

  }

}

#endif
