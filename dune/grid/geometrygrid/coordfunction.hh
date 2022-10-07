// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_COORDFUNCTION_HH
#define DUNE_GEOGRID_COORDFUNCTION_HH

#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/std/type_traits.hh>

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

  /** \brief Interface class for using an analytical function to define the
   *  geometry of a Dune::GeometryGrid. An implementation should be derived
   *  from Dune::AnalyticalCoordFunction and the evaluate
   *  method mapping \f$ R^d\to R^r \f$ has to be supplied.
   *
   *  \tparam ct coordinate field type (\c ct in Dune::GeometryGrid)
   *  \tparam dimD dimension of the domain of the mapping (\c dimension
   *               in the host grid).
   *  \tparam dimR dimension of the range of the mapping (\c dimensionworld
   *               in Dune::GeometryGrid)
   *  \tparam Impl implementation class (BN trick)
   **/
  template< class ct, unsigned int dimD, unsigned int dimR, class Impl >
  class AnalyticalCoordFunctionInterface
  {
    typedef AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > This;

    friend class AnalyticalCoordFunction< ct, dimD, dimR, Impl >;

  public:
    typedef This Interface;
    typedef Impl Implementation;

    //! field type of the coordinate vector
    typedef ct ctype;

    //! dimension of the range vector (dimensionworld of host grid)
    static const unsigned int dimDomain = dimD;
    //! dimension of the range vector
    static const unsigned int dimRange = dimR;

    //! domain vector for the evaluate method
    typedef FieldVector< ctype, dimDomain > DomainVector;
    //! range vector for the evaluate method
    typedef FieldVector< ctype, dimRange > RangeVector;

  private:
    AnalyticalCoordFunctionInterface () = default;
    AnalyticalCoordFunctionInterface ( const This & ) = default;
    AnalyticalCoordFunctionInterface ( This && ) = default;
    ~AnalyticalCoordFunctionInterface () = default;
    This &operator= ( const This & ) = default;
    This &operator= ( This && ) = default;

    // helper for picking the correct version of evaluate further down
    template<typename F, typename DV>
    using has_operator_parentheses = decltype(std::declval<F>()(std::declval<DV>()));

  public:

#ifdef DOXYGEN

    //! evaluate method for global mapping
    void evaluate ( const DomainVector &x, RangeVector &y ) const;

#else

    template<typename DV>
    std::enable_if_t<
      Std::is_detected<has_operator_parentheses,Impl,DV>::value
      >
    evaluate ( const DV &x, RangeVector &y ) const
    {
      y = asImp()(x);
    }

    template<typename DV>
    std::enable_if_t<
      not Std::is_detected<has_operator_parentheses,Impl,DV>::value
      >
    evaluate ( const DV &x, RangeVector &y ) const
    {
      assert(
        static_cast<void(This::*)(const DomainVector&, RangeVector&) const>(&This::evaluate) !=
        static_cast<void(Impl::*)(const DomainVector&, RangeVector&) const>(&Impl::evaluate) &&
        "You need to implement either operator() or evaluate() in your coordinate function!");
      asImp().evaluate( x, y );
    }

#endif // DOXYGEN

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
  /** @brief Derive an implementation of an analytical coordinate function
   * from this class.
   **/
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
    AnalyticalCoordFunction () = default;
    AnalyticalCoordFunction ( const This & ) = default;
    AnalyticalCoordFunction ( This && ) = default;
    ~AnalyticalCoordFunction () = default;
    This &operator= ( const This & ) = default;
    This &operator= ( This && ) = default;

  private:

  };



  // DiscreteCoordFunctionInterface
  // ------------------------------

  /** \brief Interface class for using a discrete function to define the
   *  geometry of a Dune::GeometryGrid. An implementation should be derived
   *  from Dune::DiscreteCoordinateFunction and the evaluate method taking an entity
   *  of the host grid together with the number of a vertex returns the coordinate in
   *  \f$ R^r \f$
   *  of that corner. The user must ensure continuity of this mapping.
   *  In addition an adapt method is provided which is called whenever
   *  \c adapt() is called on the Dune::GeometryGrid.
   *
   *  \tparam ct coordinate field type (\c ct in Dune::GeometryGrid)
   *  \tparam dimR dimension of the range of the mapping (\c dimensionworld
   *               in Dune::GeometryGrid)
   *  \tparam Impl implementation class (BN trick)
   **/
  template< class ct, unsigned int dimR, class Impl >
  class DiscreteCoordFunctionInterface
  {
    typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > This;

    friend class DiscreteCoordFunction< ct, dimR, Impl >;

  public:
    typedef This Interface;
    typedef Impl Implementation;

    //! field type of the coordinate vector
    typedef ct ctype;

    //! dimension of the range vector
    static const unsigned int dimRange = dimR;

    //! range vector for the evaluate method
    typedef FieldVector< ctype, dimRange > RangeVector;

  private:
    DiscreteCoordFunctionInterface () = default;
    DiscreteCoordFunctionInterface ( const This & ) = default;
    DiscreteCoordFunctionInterface ( This && ) = default;
    ~DiscreteCoordFunctionInterface () = default;
    This &operator= ( const This & ) = default;
    This &operator= ( This && ) = default;

  public:
    /** \brief evaluate method
     *  \param hostEntity an entity of the host grid
     *  \param corner the local number of the corner in the host entity
     *  \param y return value for the coordinate of this corner
     *
     * \note This method needs to work for entities of all codimensions,
     *    not just for elements!
     **/
    template< class HostEntity >
    void evaluate ( const HostEntity &hostEntity, unsigned int corner,
                    RangeVector &y ) const
    {
      asImp().evaluate( hostEntity, corner, y );
    }

    /** \brief method called from grid.adapt() method to allow adaptation
     *  of the discrete coordinate function
     **/
    void adapt ()
    {
      asImp().adapt();
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
  //
  /** @brief Derive an implementation of a discrete coordinate function
   * from this class.
   **/
  template< class ct, unsigned int dimR, class Impl >
  class DiscreteCoordFunction
    : public DiscreteCoordFunctionInterface< ct, dimR, Impl >
  {
    typedef DiscreteCoordFunction< ct, dimR, Impl > This;
    typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > Base;

  public:
    typedef typename Base :: RangeVector RangeVector;

  protected:
    DiscreteCoordFunction () = default;
    DiscreteCoordFunction ( const This & ) = default;
    DiscreteCoordFunction ( This && ) = default;
    ~DiscreteCoordFunction () = default;
    This &operator= ( const This & ) = default;
    This &operator= ( This && ) = default;

    void adapt ()
    {}

  private:
    template< class HostEntity >
    void evaluate ( const HostEntity &hostEntity, unsigned int corner,
                    RangeVector &y ) const;
  };



  namespace GeoGrid
  {

    // isCoordFunctionInterface
    // ------------------------

    template< class CoordFunctionInterface >
    struct isCoordFunctionInterface
    {
      static const bool value = false;
    };

    template< class ct, unsigned int dimD, unsigned int dimR, class Impl >
    struct isCoordFunctionInterface
    < AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > >
    {
      static const bool value = true;
    };

    template< class ct, unsigned int dimR, class Impl >
    struct isCoordFunctionInterface
    < DiscreteCoordFunctionInterface< ct, dimR, Impl > >
    {
      static const bool value = true;
    };



    // isDiscreteCoordFunctionInterface
    // --------------------------------

    template< class CoordFunctionInterface >
    struct isDiscreteCoordFunctionInterface
    {
      static const bool value = false;
    };

    template< class ct, unsigned int dimR, class Impl >
    struct isDiscreteCoordFunctionInterface
    < DiscreteCoordFunctionInterface< ct, dimR, Impl > >
    {
      static const bool value = true;
    };



    // AdaptCoordFunction
    // ------------------

    template< class CoordFunctionInterface >
    struct AdaptCoordFunction
    {
      static void adapt ( CoordFunctionInterface & )
      {}
    };

    template< class ct, unsigned int dimR, class Impl >
    struct AdaptCoordFunction< DiscreteCoordFunctionInterface< ct, dimR, Impl > >
    {
      typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > CoordFunctionInterface;

      static void adapt ( CoordFunctionInterface &coordFunction )
      {
        coordFunction.adapt();
      }
    };

  } // namespace GeoGrid

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_COORDFUNCTION_HH
