// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_GEOMETRY_HH
#define DUNE_GRID_COMMON_GEOMETRY_HH

/** \file
    \brief Wrapper and interface classes for element geometries
 */

#include <cassert>

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/transpose.hh>
#include <dune/common/std/type_traits.hh>

#include <dune/geometry/referenceelements.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim, int dimworld, class ct, class GridFamily >
  class GridDefaultImplementation;



  //*****************************************************************************
  //
  // Geometry
  // forwards the interface to the implementation
  //
  //*****************************************************************************

  /**
     @brief Wrapper class for geometries


     \tparam mydim Dimension of the domain
     \tparam cdim Dimension of the range
     \tparam GridImp Type that is a model of Dune::Grid
     \tparam GeometryImp Class template that is a model of Dune::Geometry

     <H3>Maps</H3>

     A Geometry defines a map \f[ g : D \to W\f] where
     \f$D\subseteq\mathbf{R}^\textrm{mydim}\f$ and
     \f$W\subseteq\mathbf{R}^\textrm{cdim}\f$.
     The domain \f$D\f$ is one of a set of predefined convex polytopes, the
     so-called reference elements (see Dune::ReferenceElement). The dimensionality
     of \f$D\f$ is <tt>mydim</tt>.
     In general \f$\textrm{mydim}\leq\textrm{cdim}\f$, i.e.
     the convex polytope may be mapped to a manifold. Moreover, we require that
     \f$ g\in \left( C^1(D) \right)^\textrm{cdim}\f$ and one-to-one.


     <H3>Engine Concept</H3>

     The Geometry class template wraps an object of type GeometryImp and forwards all member
     function calls to corresponding members of this class. In that sense Geometry
     defines the interface and GeometryImp supplies the implementation.

     \ingroup GIGeometry
   */
  template< int mydim, int cdim, class GridImp, template< int, int, class > class GeometryImp >
  class Geometry
  {
  public:
    /**
     * \brief type of underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    typedef GeometryImp< mydim, cdim, GridImp > Implementation;

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    Implementation &impl () { return realGeometry; }
    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    const Implementation &impl () const { return realGeometry; }

    //! @brief geometry dimension
    constexpr static int mydimension = mydim;

    //! @brief dimension of embedding coordinate system
    constexpr static int coorddimension = cdim;

    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! type of local coordinates
    typedef FieldVector<ctype, mydim> LocalCoordinate;

    //! type of the global coordinates
    typedef FieldVector< ctype, cdim > GlobalCoordinate;

    //! Number type used for the geometry volume
    typedef decltype(std::declval<Implementation>().volume()) Volume;

    /**
     * \brief type of jacobian inverse transposed
     *
     * The exact type is implementation-dependent.
     * However, it is guaranteed to have the following properties:
     * - It satisfies the ConstMatrix interface.
     * - It is copy constructible and copy assignable.
     * .
     */
    typedef typename Implementation::JacobianInverseTransposed JacobianInverseTransposed;

    /**
     * \brief type of jacobian transposed
     *
     * The exact type is implementation-dependent.
     * However, it is guaranteed to have the following properties:
     * - It satisfies the ConstMatrix interface.
     * - It is copy constructible and copy assignable.
     * .
     */
    typedef typename Implementation::JacobianTransposed JacobianTransposed;

  private:

    template<class Implementation_T>
    using JacobianInverseOfImplementation = decltype(typename Implementation_T::JacobianInverse{std::declval<Implementation_T>().jacobianInverse(std::declval<LocalCoordinate>())});

    using JacobianInverseDefault = decltype(transpose(std::declval<JacobianInverseTransposed>()));

    template<class Implementation_T>
    using JacobianOfImplementation = decltype(typename Implementation_T::Jacobian{std::declval<Implementation_T>().jacobian(std::declval<LocalCoordinate>())});

    using JacobianDefault = decltype(transpose(std::declval<JacobianTransposed>()));

    auto jacobianImpl ( const LocalCoordinate &local, std::true_type /*implDetected*/ ) const
    {
      return impl().jacobian(local);
    }

    [[deprecated("Geometry implementatons are required to provide a jacobian(local) method. The default implementation is deprecated and will be removed after release 2.9")]]
    auto jacobianImpl ( const LocalCoordinate &local, std::false_type /*implNotDetected*/ ) const
    {
      return transpose(jacobianTransposed(local));
    }

    auto jacobianInverseImpl ( const LocalCoordinate &local, std::true_type /*implDetected*/ ) const
    {
      return impl().jacobianInverse(local);
    }

    [[deprecated("Geometry implementatons are required to provide a jacobianInverse(local) method. The default implementation is deprecated and will be removed after release 2.9")]]
    auto jacobianInverseImpl ( const LocalCoordinate &local, std::false_type /*implNotDetected*/ ) const
    {
      return transpose(jacobianInverseTransposed(local));
    }

  public:

    /**
     * \brief type of jacobian inverse
     *
     * The exact type is implementation-dependent.
     * However, it is guaranteed to have the following properties:
     * - You can multilply it from the right to a suitable FieldMatrix
     * - It is copy constructible and copy assignable.
     * .
     */
    using JacobianInverse = Std::detected_or_t<JacobianInverseDefault, JacobianInverseOfImplementation, Implementation>;

    /**
     * \brief type of jacobian
     *
     * The exact type is implementation-dependent.
     * However, it is guaranteed to have the following properties:
     * - You can multilply it from the right to a suitable FieldMatrix
     * - It is copy constructible and copy assignable.
     * .
     */
    using Jacobian = Std::detected_or_t<JacobianDefault, JacobianOfImplementation, Implementation>;

    /** \brief Return the type of the reference element. The type can
       be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return impl().type(); }

    /** \brief Return true if the geometry mapping is affine and false otherwise */
    bool affine() const { return impl().affine(); }

    /** \brief Return the number of corners of the reference element.
     *
       Since a geometry is a convex polytope the number of corners is a well-defined concept.
       The method is redundant because this information is also available
       via the reference element. It is here for efficiency and ease of use.
     */
    int corners () const { return impl().corners(); }

    /** \brief Obtain a corner of the geometry
     *
     *  This method is for convenient access to the corners of the geometry. The
     *  same result could be achieved by calling
        \code
        global( referenceElement.position( i, mydimension ) )
        \endcode
     *
     *  \param[in]  i  number of the corner (with respect to the reference element)
     *
     *  \returns position of the i-th corner
     */
    GlobalCoordinate corner ( int i ) const
    {
      return impl().corner( i );
    }

    /** \brief Evaluate the map \f$ g\f$.
       \param[in] local Position in the reference element \f$D\f$
       \return Position in \f$W\f$
     */
    GlobalCoordinate global (const LocalCoordinate& local) const
    {
      return impl().global( local );
    }

    /** \brief Evaluate the inverse map \f$ g^{-1}\f$
       \param[in] global Position in \f$W\f$
       \return Position in \f$D\f$ that maps to global
     */
    LocalCoordinate local (const GlobalCoordinate& global) const
    {
      return impl().local( global );
    }

    /** \brief Return the factor appearing in the integral transformation formula

       Let \f$ g : D \to W\f$ denote the transformation described by the Geometry.
       Then the jacobian of the transformation is defined as the
       \f$\textrm{cdim}\times\textrm{mydim}\f$ matrix
       \f[ J_g(x) = \left( \begin{array}{ccc} \frac{\partial g_0}{\partial x_0} &
       \cdots & \frac{\partial g_0}{\partial x_{n-1}} \\
       \vdots & \ddots & \vdots \\ \frac{\partial g_{m-1}}{\partial x_0} &
       \cdots & \frac{\partial g_{m-1}}{\partial x_{n-1}}
       \end{array} \right).\f]
       Here we abbreviated \f$m=\textrm{cdim}\f$ and \f$n=\textrm{mydim}\f$ for ease of
       readability.

       The integration element \f$\mu(x)\f$ for any \f$x\in D\f$ is then defined as
       \f[ \mu(x) = \sqrt{|\det J_g^T(x)J_g(x)|}.\f]

       \param[in] local Position \f$x\in D\f$
       \return    integration element \f$\mu(x)\f$

       \note Each implementation computes the integration element with optimal
       efficiency. For example in an equidistant structured mesh it may be as
       simple as \f$h^\textrm{mydim}\f$.
     */
    Volume integrationElement (const LocalCoordinate& local) const
    {
      return impl().integrationElement( local );
    }

    /** \brief return volume of geometry */
    Volume volume () const
    {
      return impl().volume();
    }

    /** \brief return center of geometry
     *
     *  Note that this method is still subject to a change of name and semantics.
     *  At the moment, the center is not required to be the centroid of the
     *  geometry, or even the centroid of its corners. This makes the current
     *  default implementation acceptable, which maps the centroid of the
     *  reference element to the geometry.
     *  We may change the name (and semantic) of the method to centroid() if we
     *  find reasonably efficient ways to implement it properly.
     */
    GlobalCoordinate center () const
    {
      return impl().center();
    }

    /** \brief Return the transposed of the Jacobian
     *
     *  The Jacobian is defined in the documentation of
     *  \ref Dune::Geometry::integrationElement "integrationElement".
     *
     *  \param[in]  local  position \f$x\in D\f$
     *
     *  \return \f$J_g^T(x)\f$
     *
     *  \note The exact return type is implementation defined.
     */
    JacobianTransposed jacobianTransposed ( const LocalCoordinate& local ) const
    {
      return impl().jacobianTransposed( local );
    }

    /** \brief Return inverse of transposed of Jacobian
     *
     *  The Jacobian is defined in the documentation of
     *  \ref Dune::Geometry::integrationElement "integrationElement".
     *
     *  \param[in]  local  position \f$x\in D\f$
     *  \return \f$J_g^{-T}(x)\f$
     *
     *  The use of this function is to compute the gradient of some function
     *  \f$f : W \to \textbf{R}\f$ at some position \f$y=g(x)\f$, where
     *  \f$x\in D\f$ and \f$g\f$ the transformation of the Geometry.
     *  When we set \f$\hat{f}(x) = f(g(x))\f$ and apply the chain rule we obtain
     *  \f[\nabla f(g(x)) = J_g^{-T}(x) \nabla \hat{f}(x).\f]
     *
     *  \note In the non-quadratic case \f$\textrm{cdim} \neq \textrm{mydim}\f$, the
     *        pseudoinverse of \f$J_g^T(x)\f$ is returned.
     *        This means that it is inverse for all tangential vectors in
     *        \f$g(x)\f$ while mapping all normal vectors to zero.
     *
     *  \note The exact return type is implementation defined.
     */
    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return impl().jacobianInverseTransposed(local);
    }

    /** \brief Return the Jacobian
     *
     *  The Jacobian is defined in the documentation of
     *  \ref Dune::Geometry::integrationElement "integrationElement".
     *
     *  \param[in]  local  position \f$x\in D\f$
     *
     *  \return \f$J_g(x)\f$
     *
     *  \note The exact return type is implementation defined.
     */
    Jacobian jacobian ( const LocalCoordinate& local ) const
    {
      auto implDetected = Std::is_detected<JacobianOfImplementation, Implementation>{};
      return jacobianImpl(local, implDetected);
    }

    /** \brief Return inverse of Jacobian
     *
     *  The Jacobian is defined in the documentation of
     *  \ref Dune::Geometry::integrationElement "integrationElement".
     *
     *  \param[in]  local  position \f$x\in D\f$
     *  \return \f$J_g^{-1}(x)\f$
     *
     *  The use of this function is to compute the jacobians of some function
     *  \f$f : W \to \textbf{R}\f$ at some position \f$y=g(x)\f$, where
     *  \f$x\in D\f$ and \f$g\f$ the transformation of the Geometry.
     *  When we set \f$\hat{f}(x) = f(g(x))\f$ and apply the chain rule we obtain
     *  \f[\nabla f(g(x)) = J_{\hat{f}}(x) J_g^{-1}(x).\f]
     *
     *  \note In the non-quadratic case \f$\textrm{cdim} \neq \textrm{mydim}\f$, the
     *        pseudoinverse of \f$J_g(x)\f$ is returned.
     *        This means that its transposed is inverse for all tangential vectors in
     *        \f$g(x)\f$ while mapping all normal vectors to zero.
     *
     *  \note The exact return type is implementation defined.
     */
    JacobianInverse jacobianInverse ( const LocalCoordinate &local ) const
    {
      auto implDetected = Std::is_detected<JacobianInverseOfImplementation, Implementation>{};
      return jacobianInverseImpl(local, implDetected);
    }

    //===========================================================
    /** @name Interface for grid implementers
     */
    //@{
    //===========================================================

    //! copy constructor from implementation
    explicit Geometry ( const Implementation &impl )
      : realGeometry( impl )
    {}

    //@}

  protected:

    Implementation realGeometry;
  };



  //************************************************************************
  // GEOMETRY Default Implementations
  //*************************************************************************
  //
  // --GeometryDefault
  //
  //! Default implementation for class Geometry
  template<int mydim, int cdim, class GridImp, template<int,int,class> class GeometryImp>
  class GeometryDefaultImplementation
  {
  public:
    static const int mydimension = mydim;
    static const int coorddimension = cdim;

    // save typing
    typedef typename GridImp::ctype ctype;

    typedef FieldVector< ctype, mydim > LocalCoordinate;
    typedef FieldVector< ctype, cdim > GlobalCoordinate;

    //! Number type used for the geometry volume
    typedef ctype Volume;

    //! type of jacobian inverse transposed
    typedef FieldMatrix< ctype, cdim, mydim > JacobianInverseTransposed;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

    //! type of jacobian inverse
    typedef FieldMatrix< ctype, mydim, cdim > JacobianInverse;

    //! type of jacobian
    typedef FieldMatrix< ctype, cdim, mydim > Jacobian;

    //! return volume of the geometry
    Volume volume () const
    {
      GeometryType type = asImp().type();

      // get corresponding reference element
      auto refElement = referenceElement< ctype , mydim >(type);

      LocalCoordinate localBaryCenter ( 0 );
      // calculate local bary center
      const int corners = refElement.size(0,0,mydim);
      for(int i=0; i<corners; ++i) localBaryCenter += refElement.position(i,mydim);
      localBaryCenter *= (ctype) (1.0/corners);

      // volume is volume of reference element times integrationElement
      return refElement.volume() * asImp().integrationElement(localBaryCenter);
    }

    //! return center of the geometry
    GlobalCoordinate center () const
    {
      GeometryType type = asImp().type();

      // get corresponding reference element
      auto refElement = referenceElement< ctype , mydim >(type);

      // center is (for now) the centroid of the reference element mapped to
      // this geometry.
      return asImp().global(refElement.position(0,0));
    }

    //! Return the Jacobian
    Jacobian jacobian ( const LocalCoordinate& local ) const
    {
      return asImp().jacobianTransposed(local).transposed();
    }

    //! Return inverse of Jacobian
    JacobianInverse jacobianInverse ( const LocalCoordinate &local ) const
    {
      return asImp().jacobianInverseTransposed(local).transposed();
    }

  private:
    //!  Barton-Nackman trick
    GeometryImp<mydim,cdim,GridImp>& asImp () {return static_cast<GeometryImp<mydim,cdim,GridImp>&>(*this);}
    const GeometryImp<mydim,cdim,GridImp>& asImp () const {return static_cast<const GeometryImp<mydim,cdim,GridImp>&>(*this);}
  }; // end GeometryDefault

  template<int cdim, class GridImp, template<int,int,class> class GeometryImp>
  class GeometryDefaultImplementation<0,cdim,GridImp,GeometryImp>
  {
    // my dimension
    constexpr static int mydim = 0;

  public:
    static const int mydimension = mydim;
    static const int coorddimension = cdim;

    // save typing
    typedef typename GridImp::ctype ctype;

    typedef FieldVector< ctype, mydim > LocalCoordinate;
    typedef FieldVector< ctype, cdim > GlobalCoordinate;
    typedef ctype Volume;

    //! type of jacobian inverse transposed
    typedef FieldMatrix< ctype, cdim, mydim > JacobianInverseTransposed;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

    //! type of jacobian inverse
    typedef FieldMatrix< ctype, mydim, cdim > JacobianInverse;

    //! type of jacobian
    typedef FieldMatrix< ctype, cdim, mydim > Jacobian;

    //! return the only coordinate
    FieldVector<ctype, cdim> global (const FieldVector<ctype, mydim>& local) const
    {
      return asImp().corner(0);
    }

    //! return empty vector
    FieldVector<ctype, mydim> local (const FieldVector<ctype, cdim>& ) const
    {
      return FieldVector<ctype, mydim>();
    }

    //! return volume of the geometry
    Volume volume () const
    {
      return Volume(1.0);
    }

    //! return center of the geometry
    FieldVector<ctype, cdim> center () const
    {
      return asImp().corner(0);
    }

    //! Return the Jacobian
    Jacobian jacobian ( const LocalCoordinate& local ) const
    {
      return asImp().jacobianTransposed(local).transposed();
    }

    //! Return inverse of Jacobian
    JacobianInverse jacobianInverse ( const LocalCoordinate &local ) const
    {
      return asImp().jacobianInverseTransposed(local).transposed();
    }

  private:
    // Barton-Nackman trick
    GeometryImp<mydim,cdim,GridImp>& asImp () {return static_cast<GeometryImp<mydim,cdim,GridImp>&>(*this);}
    const GeometryImp<mydim,cdim,GridImp>& asImp () const {return static_cast<const GeometryImp<mydim,cdim,GridImp>&>(*this);}
  }; // end GeometryDefault



  // referenceElement
  // ----------------

  template< int mydim, int cdim, class GridImp, template< int, int, class > class GeometryImp>
  auto referenceElement(const Geometry<mydim,cdim,GridImp,GeometryImp>& geo)
    -> decltype(referenceElement(geo,geo.impl()))
  {
    return referenceElement(geo,geo.impl());
  }

  //! Second-level dispatch to select the correct reference element for a grid geometry.
  /**
   * This function is the default implementation of the second-level reference element dispatch
   * performed by Geometry.
   *
   * \note This function is only important for grid implementors with geometries that do not have
   *       a standard reference element.
   *
   * When referenceElement() is called with a Geometry, it will forward the call to
   * `referenceElement(const Geometry&,const GeometryImplementation&)`. This default implementation
   * will do the right thing as long as your geometry is based on a standard Dune ReferenceElement. If
   * it is not and you want to supply your own reference element implementation, provide an override of
   * this function for your specific geometry implementation.
   *
   * \related Geometry
   */
  template< int mydim, int cdim, class GridImp, template< int, int, class > class GeometryImp, typename Impl>
  auto referenceElement(const Geometry<mydim,cdim,GridImp,GeometryImp>& geo, const Impl&)
    -> decltype(referenceElement<typename GridImp::ctype,mydim>(geo.type()))
  {
    using Geo = Geometry<mydim,cdim,GridImp,GeometryImp>;
    return referenceElement<typename Geo::ctype,Geo::mydimension>(geo.type());
  }

} // namespace Dune

#endif // DUNE_GRID_COMMON_GEOMETRY_HH
