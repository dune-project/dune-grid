// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_COORDFUNCTION_HH
#define DUNE_GEOGRID_COORDFUNCTION_HH

#include <dune/common/fvector.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ct, unsigned int dimD, unsigned int dimR, class Impl >
  class AnalyticalCoordFunction;

  template< class ct, unsigned int dimR, class Impl >
  class DiscreteCoordFunction;



  // AnalyticalCoordFunctionInterface
  // --------------------------------

  template< class ct, unsigned int dimD, unsigned int dimR, class Impl >
  class AnalyticalCoordFunctionInterface
  {
    typedef AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > This;

    friend class AnalyticalCoordFunction< ct, dimD, dimR, Impl >;

  public:
    typedef This Interface;
    typedef Impl Implementation;

    typedef ct ctype;

    static const unsigned int dimDomain = dimD;
    static const unsigned int dimRange = dimR;

    typedef FieldVector< ctype, dimDomain > DomainVector;
    typedef FieldVector< ctype, dimRange > RangeVector;

  private:
    AnalyticalCoordFunctionInterface ()
    {}

    AnalyticalCoordFunctionInterface ( const This & );
    This &operator= ( const This & );

  public:
    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      return asImp().evaluate( x, y );
    }

  protected:
    const Implementation &asImp () const
    {
      return static_cast< const Implementation & >( *this );
    }

    Implementation &asImp ()
    {
      return static_cast< Implementation & >( *this );
    }
  };



  // AnalyticalCoordFunction
  // -----------------------

  template< class ct, unsigned int dimD, unsigned int dimR, class Impl >
  class AnalyticalCoordFunction
    : public AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl >
  {
    typedef AnalyticalCoordFunction< ct, dimD, dimR, Impl > This;
    typedef AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > Base;

  public:
    typedef typename Base :: DomainVector DomainVector;
    typedef typename Base :: RangeVector RangeVector;

  protected:
    AnalyticalCoordFunction ()
    {}

  private:
    AnalyticalCoordFunction ( const This & );
    This &operator= ( const This & );

    void evaluate ( const DomainVector &x, RangeVector &y ) const;
  };



  // DiscreteCoordFunctionInterface
  // ------------------------------

  template< class ct, unsigned int dimR, class Impl >
  class DiscreteCoordFunctionInterface
  {
    typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > This;

    friend class DiscreteCoordFunction< ct, dimR, Impl >;

  public:
    typedef This Interface;
    typedef Impl Implementation;

    typedef ct ctype;

    static const unsigned int dimRange = dimR;

    typedef FieldVector< ctype, dimRange > RangeVector;

  private:
    DiscreteCoordFunctionInterface ()
    {}

    DiscreteCoordFunctionInterface ( const This & );

    This &operator= ( const This & );

  public:
    template< class Index >
    void evaluate ( const Index &vIdx, RangeVector &y ) const
    {
      return asImp().evaluate( vIdx, y );
    }

  protected:
    const Implementation &asImp () const
    {
      return static_cast< const Implementation & >( *this );
    }

    Implementation &asImp ()
    {
      return static_cast< Implementation & >( *this );
    }
  };



  // DiscreteCoordFunction
  // ---------------------

  template< class ct, unsigned int dimR, class Impl >
  class DiscreteCoordFunction
    : public DiscreteCoordFunctionInterface< ct, dimR, Impl >
  {
    typedef DiscreteCoordFunction< ct, dimR, Impl > This;
    typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > Base;

  public:
    typedef typename Base :: RangeVector RangeVector;

  protected:
    DiscreteCoordFunction ()
    {}

  private:
    DiscreteCoordFunction ( const This & );
    This &operator= ( const This & );

    template< class Index >
    void evaluate ( const Index &vIdx, RangeVector &y ) const;
  };

}

#endif
