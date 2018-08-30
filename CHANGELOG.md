# master (will become 2.7)

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
