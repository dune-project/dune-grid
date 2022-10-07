<!--
SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
-->

# Master (will become release 2.9)

- The `Geometry` interface was extended by methods `jacobian(local)` and `jacobianInverse(local)`
  and corresponding typedefs `Jacobian` and `JacobianInverse`. All grid implementations need to
  provide the new interface. For transition, the methods and typedefs are default-implemented
  in the `Dune::Geometry` interface class which is used for all grid geometries.

- The `Geometry::integrationElement` now needs to return the type `Volume`
  instead of `ctype`. Note that this may be different from `ctype` if the grid
  supports typed dimensions. In such a case, `ctype` is a length, and not
  appropriate for a volume quantity.

- The `FindAlberta.cmake` module only searches for world dimension libraries up to dim 3. This
  can be increased by setting the CMake variable `ALBERTA_MAX_WORLD_DIM`

- A method `GridView::isConforming()` was added
  indicating whether a given grid view is representing a conforming grid.
  For some grid implementations this could vary from the compile time static
  `GridView::conforming` depending on run time parameters.

- The Python bindings now support tensorproduct coordinates of `YaspGrid`

## Deprecations and removals

- Support for PSurface has been removed.

- Support for AmiraMesh has been removed.

- Following the deprecation of `CollectiveCommunication` in dune-common, grids define
  the type `Communication`. The type `CollectiveCommunication` is deprecated and
  will be removed after Dune 2.9.

- The CMake function `add_dune_ug_flags` is deprecated. It should not be used since
  all flags are already set automatically. Also the package flags for the UG grid are
  not registered any more.

- The `dune-uggrid` module does not set the preprocessor flag `HAVE_UG` anymore.
  Use `HAVE_DUNE_UGGRID` instead.

# Release 2.8

- Specialization of `StructuredGridFactory` for `AlbertaGrid` is added to address the special
  numbering requirements of that grid.

- Return type of `GridFactory::createGrid()` changed to `std::unique_ptr`. While still the
  replacement type `ToUniquePtr` from dune-common works, it is marked deprecated and the std
  implementation should be used directly now.

- `UGGrid` removes support for `_2` and `_3` macros.

- `SingleCodimSingleGeomTypeMapper` and `MultipleCodimMultipleGeomTypeMapper` now have an `update(gridView)`
  method to update the stored `GridView` and recalculate the indices after mesh adaptation.

- The `update()` member function of mappers is deprecated. Use the new
  member function `update(gridView)` with a grid view argument when updating
  the mapper after the grid or grid view changes. The interface
  change reflects that grid view has value semantics.

- The "convenience" classes `LeafSingleCodimSingleGeomTypeMapper`,
  `LevelSingleCodimSingleGeomTypeMapper`, `LeafMultipleCodimMultipleGeomTypeMapper`,
  `LevelMultipleCodimMultipleGeomTypeMapper` have been deprecated since they
  don't comply with the new mapper interface.
  Just use `SingleCodimSingleGeomTypeMapper` and `MultipleCodimMultipleGeomTypeMapper`.

- `IdSet` now exports grid `dimension` and `Codim<cd>::Entity`.

- `UGGrid` index sets can now compute the indices of vertices of edges.

- `UGGrid`: Fixed a bug in the numbering of prism edges.

- Various bugs have been fixed in the `UGGrid` subdomain communication implementation.

- `YaspGrid` now supports class template argument deduction
   (https://en.cppreference.com/w/cpp/language/class_template_argument_deduction)
   if you have a C++17-CTAD-compatible compiler (gcc >= 8, clang 5)

- Python bindings have been moved from the `dune-python` module which is now
  obsolete. To activate Python bindings the CMake flag
  `DUNE_ENABLE_PYTHONBINDINGS` needs to be turned on (default is off).
  Furthermore, flags for either shared library or position independent code
  needs to be used.
- Properly implement the `canCommunicate` capability for `UGGrid` and `IdentityGrid`.

- The return type of the `IndexSet::size` methods isn't `IndexType`
  anymore. In general the return type should be an unsigned integral
  type. The actual type is implementation specific. All
  implementations in `dune-grid` now return `std::size_t`, following
  our approach to make all size information be unsigned.

- The grid capability `hasEntityIterator` now defaults to whatever
  the capability `hasEntity` specifies. Most grid manager implementers
  now only need to implement/specialize `hasEntity`.

- Add a new range generator `subEntities` that allows to iterate over the sub-entities
  of a specified codimension for a given element (codim-0 entity)
  ```cpp
  for (const auto& vertex : subEntities(element, Codim<dim>))
    std::cout << vertex.geometry().corner(0) << "\n";
  ```
  The returned range is sized, i.e. `subEntities(element, Codim<dim>).size()`
  returns the number of vertices in the element and is equivalent to
  `referenceElement(element).size(dim)`.

## Deprecations and removals

- Remove `Intersection`'s deprecated enums `dimension` and
  `codimension`. Instead use grid's dimension and 1.

- Remove deprecated `Grid::getRealImplementation`.
  Use the `impl()` member function of the facade classes directly instead.

- Remove GeometryGrid's deprecated constructors accepting raw pointers.

- Remove deprecated `AlbertaGrid::getRealIntersection`.
  Use the `impl()` member function of the intersection class directly instead.

- Remove deprecated `AlbertaGrid::readGridXdr` and `AlbertaGrid::writeGridXdr`.
  Instead use `AlbertaGrid::readGrid` and `AlbertaGrid::writeGrid`.

- Remove deprecated header `dune/common/universalmapper.hh`.

- Support for PSurface is deprecated and will be removed after Dune 2.8.

- Support for AmiraMesh is deprecated and will be removed after Dune 2.8.

# Release 2.7

- The `YaspGrid` class has a new constructor that takes a `Coordinates`
  object as its first argument.  This object can be of type `EquidistantCoordinates`,
  `EquidistantOffsetCoordinates`, or `TensorProductCoordinates`,
  and encapsulates the domain and element sizes.  Previously, this data
  was given directly to different constructors of `YaspGrid`,
  and the constructor had to match the `Coordinates` object used
  as the second template argument of the `YaspGrid` class.
  The new constructor subsumes the previous three, and makes sure
  by design that the correct information is passed at `YaspGrid`
  construction.

- The `GridFactory` exports the type `Communication` and corresponding object
  `comm()` that is used to create the grid. In the `GridFactoryInterface`, this
  defaults to the Communication induced by the process-local communicator
  `MPIHelper::getLocalCommunicator()`, while `UGGrid`'s factory returns the
  Communication of the grid. This can be used to steer the grid creation
  process, see `dune/grid/io/file/gmshreader.hh` for an example.

- The number type used by a `BoundarySegment` is not hard-wired to `double`
  anymore.

- The `Grid::getRealImplementation` member function has been deprecated.
  Use the `impl()` member function of the facade classes directly instead.

- The `AlbertaGrid::getRealIntersection` member function has been deprecated.
  Use the `impl()` member function of the intersection class directly instead.

- The return type of all file reader methods that return a grid has been changed
  from a plain C pointer to the custom pointer class `ToUniquePtr<Grid>`
  (from the `dune-common` module).  Values of this pointer class cast to
  C pointers (but with a deprecation warning), to `std::unique_ptr`, and to
  `std::shared_ptr`.  This marks the beginning of a transition period.
  In the long run, the methods are planned to return objects of type
  `std::unique_ptr`, to make it obvious that the calling code receives
  the ownership of the grid object.  For the time being the calling code
  can still store the return value in a C pointer, but that possibility
  will go away.

- Likewise, the return type of the `GridFactory::createGrid`method has been changed
  from a plain C pointer to the custom pointer class `ToUniquePtr<Grid>`.
  In the long run, the method is planned to return objects of type
  `std::unique_ptr`, to make it obvious that the calling code receives
  the ownership of the grid object.  For the time being the calling code
  can still store the return value in a C pointer, but that possibility
  will go away.  While this procedure allows full backward compatibility
  for code that calls `GridFactory::createGrid`, implementors or third-party
  grid implementations will need to update their implementations of
  `GridFactory::createGrid`.

- The `VTKWriter`s now support custom output data precision for functions
  (via the provided `FieldInfo`), `VTKFunction`s  (via a new virtual interface method `precision()`),
  and coordinates (via an argument to the writer's constructor) that can be
  chosen at runtime. Any field can now choose between VTK's `Float32`, `Float64`,
  `Int32`, `UInt8`, and `UInt32`, represented as `Dune::VTK::Precision::float32`/`float64`/`int32`/`uint8`/`uint32`,
  the `DataArrayWriter` selects the correct type at runtime. The default for
  functions and coordinate is `Float32` as before.

- The `VTKWriter` now supports writing functions that can only be evaluated globally,
  i.e., functions that are not defined with respect to the grid.  Such functions
  will be sampled on the grid vertices (`addVertexData`) or grid element centers (`addCellData`).

- The Capability `hasBackupRestoreFacilities<GeometryGrid<HG, CoordFunction>>`
  now returns `false` in case the `CoordFunction` is not default-constructible.

- The `Geometry` interface now provides the type `Volume` for the return value of the
  method of the same name.  Note that this may be different from `ctype` if you care
  about dimensions.  In that case `ctype` is a length, and not appropriate for
  a quantity that is a volume.

- The `VTKWriter` writer now truncates subnormal floating point values to 0 when writing ASCII files
  (`DUNE::VTK::ascii`). This avoids Paraview crashes on macOS. For this reasons, most VTK files written
  by DUNE 2.7 will differ from the same file written in DUNE 2.6. If you are using VTK files for testing
  results, make sure to use fuzzy float comparisons!

- The `VTKSequenceWriter` now exposes the function `clear()` of the associated `VTKWriter`

- The `VTKSequenceWriter` allows to get and set the time steps storage, which enables serialization
  of a sequence writer. The feature can be used to restart a sequence.

- UG 3.x is no longer supported. Use dune-uggrid instead.


# Release 2.6

- The deprecated `EntityPointer` has been removed completely and `EntityIterator`
  no longer inherits from it.
  As a consequence, the dimension `EntityIterator::dimension`,
  `EntityIterator::codimension`, and `EntityIterator::mydimension` are gone.

- Experimental grid extensions are now always enabled:

    See core/dune-grid!155

  - The method `impl` and the type `Implementation` on the facade classes are
    always public (and documented), now.
    Warning: Implementation details may change without prior notification.
  - The method experimental grid extension `boundaryId` has been removed from the
    intersection interface. Some grid will continue providing it on their
    implementation, i.e., it may still be accessible through
    ```
    intersection.impl().boundaryId()
    ```
  - The DGF block `general` is now always available and
    the DGFWriter will always write a boundary id and can write user-defined
    boundary data, now.

- `MultipleCodimMultipleGeomTypeMapper`: The `Layout` template parameter has
  been deprecated in favor of a function object that indicates which geometry
  types to include in the mapping.  The layout function
  object is passed in the constructor, so instead of
  ```c++
  MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout> mapper1(gv);
  MultipleCodimMultipleGeomTypeMapper<GV, MCMGVertexLayout> mapper2(gv);
  ```
  please write
  ```c++
  MultipleCodimMultipleGeomTypeMapper<GV> mapper1(gv, mcmgElementLayout());
  MultipleCodimMultipleGeomTypeMapper<GV> mapper2(gv, mcmgVertexLayout());
  ```
  See the doxygen documentation for custom layouts and core/dune-grid!177

- The `MCMGMapper` can now be used to attach multiple dofs to each
  entity:

    See core/dune-grid!215

  - the Layout is passed into the constructor and
    returns the number of dofs to attach to the given geometry type
    ```
       MCMGLayout layout = [](GeometryType gt, int griddim) {
         return gt.dim() == griddim? 2:0;
       };
       MCMGMapper mapper(grid,layout);
    ```
    Note: the layout can still return a `bool` with `true` leading to a single dof being attached.
  - The new method `MCMGMapper::indices(entity)` returns an iterable range
    (instance of `IntegralRange<Index>`)
    with the indices of dofs attached to the given entity:
    ```
      for (const auto& i : mapper.indices(entity) )
        dof = vector[i];
    ```

- Two new method were added to the MCMGMapper:
  `size_type size(GeometryType)` and
  `const std::vector< GeometryType >& types (int codim)`
  returning the number of dofs attached to the geometry type and a vector
  with all geometry types on which dofs are attached, respectively.

    See core/dune-grid!215

- The `StructuredGridFactory` now returns a `unique_ptr` instead of a
  `shared_ptr`.  Code that relies on a `shared_ptr`
  needs to explicitly assign the return value to a `shared_ptr`
  variable.

    See core/dune-grid!212

- `SubsamplingVTKWriter` now supports arbitrary refinements, not just powers
  of two.  The old constructor taking a parameter `int levels` has been
  deprecated, you should now pass a parameter `RefinementIntervals intervals`
  instead.  There are convenience functions `refinementIntervals(int
  intervals)` and `refinementLevels(int levels)` to construct parameters of
  type `RefinementIntervals` in dune-geometry.

    See core/dune-grid!193

- `UGGrid` now supports transferring element data during load balancing.

    See core/dune-grid!172

- A convenience method `referenceElement( entity )` was added.

    See core/dune-grid!349
