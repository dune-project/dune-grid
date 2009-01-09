// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MISC_HH
#define DUNE_GENERICGEOMETRY_MISC_HH

#include <dune/common/static_assert.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< unsigned int n >
    struct Faculty
    {
      enum { value = n * Faculty< n-1 > :: value };
    };

    template<>
    struct Faculty< 0 >
    {
      enum { value = 1 };
    };



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
        Operation< first > :: apply();
        ForLoop< Operation, first+1, last > :: apply();
      }

      template< class Type >
      static void apply ( Type &param )
      {
        Operation< first > :: apply( param );
        ForLoop< Operation, first+1, last > :: apply( param );
      }

      template< class Type1, class Type2 >
      static void apply ( Type1 &param1, Type2 &param2 )
      {
        Operation< first > :: apply( param1, param2 );
        ForLoop< Operation, first+1, last > :: apply( param1, param2 );
      }

      template< class Type1, class Type2, class Type3 >
      static void apply ( Type1 &param1, Type2 &param2, Type3 &param3 )
      {
        Operation< first > :: apply( param1, param2, param3 );
        ForLoop< Operation, first+1, last > :: apply( param1, param2, param3 );
      }

    private:
      dune_static_assert( (first <= last), "ForLoop: first > last" );
    };

    template< template< int > class Operation, int last >
    struct ForLoop< Operation, last, last >
    {
      static void apply ()
      {
        Operation< last > :: apply();
      }

      template< class Type >
      static void apply ( Type &param )
      {
        Operation< last > :: apply( param );
      }

      template< class Type1, class Type2 >
      static void apply ( Type1 &param1, Type2 &param2 )
      {
        Operation< last > :: apply( param1, param2 );
      }

      template< class Type1, class Type2, class Type3 >
      static void apply ( Type1 &param1, Type2 &param2, Type3 &param3 )
      {
        Operation< last > :: apply( param1, param2, param3 );
      }
    };

  }

}

#endif
