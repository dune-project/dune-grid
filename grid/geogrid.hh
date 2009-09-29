// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/geogrid/grid.hh>

/** \addtogroup GeoGrid
 *
 *  The GeometryGrid is an implementation of the DUNE grid interface that can
 *  wrap any other DUNE grid (called host grid) and replace its geometry.
 *  To this end, the grid also gets a coordinate function that maps the corners
 *  of the host grid into any larger Euklidian space.
 *  The generic geometries are then used to provide a geometry implementation
 *  for the grid, interpolating the corners in a linear (respectively n-linear)
 *  manner.
 *
 *  \image html helix.png
 *  The figure above displays a <tt>GeometryGrid< YaspGrid< 2 >, Helix ></tt>,
 *  where Helix models the following coordinate function:
 *  \f[
 *    \left(\begin{array}{c}r\\\varphi\end{array}\right)
 *    \mapsto
 *    \left(\begin{array}{c}
 *      (r + \frac{1}{5}) \cos( 2 \pi \varphi )\\
 *      (r + \frac{1}{5}) \sin( 2 \pi \varphi )\\
 *      \varphi
 *    \end{array}\right).
 *  \f]
 *  Though YaspGrid can only model plane, Carthesian grids, using GeometryGrid
 *  we have obtained a parallel surface grid with quadrilateral elements.
 *
 *  \section features Features
 *
 *  Features of the GeometryGrid include:
 *  - complete wrapper of the host grid
 *    (i.e., no non-geometric feature of the host grid is lost)
 *  - Only uses the coordinate of the corners of each entity -
 *    no other geometric information of the underlying grid is used.
 *  - provides entities for all codimensions, even if the host grid does not
 *    (though communication is not extended to these codimensions)
 *  .
 *
 *  \section usage Usage
 *
 *  There are three different construction mechanisms for a geometry grid:
 *  - Given a host grid instance and a function mapping
 *    global coordinates from the host grid to some space
 *    with larger or equal dimension.
 *  - Given a vector class assining each index of the host grid a coordinate
 *    vector.
 *  - given an implementation of a local function container, i.e., a class with
 *    a method \c localFunction taking a entity of the host grid and returning
 *    a local function object with a \c evaluate method mapping local
 *    coordinates to global coordinates.
 *    It is required, that the resulting global mapping is continuous.
 *  .
 *  Remark: in the second case no geometry class has to be implemented by the
 *          host grid.
 *          In the first case the host grid must provide an implementation of
 *          the method <tt>corner</tt> on the geometry class for codimension
 *          zero entity.
 */
