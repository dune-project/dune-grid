# master (will become 2.7)

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
