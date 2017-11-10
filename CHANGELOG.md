# Release 2.6

- Experimental grid extensions are now always enabled.

- The method `impl` and the type `Implementation` on the facade classes are
  always public (and documented), now.
  Warning: Implementation details may change without prior notification.

- The method experimental grid extension `boundaryId` has been removed from the
  intersection interface. Some grid will continue providing it on their
  implementation, i.e., it may still be accessible through
  ```
  intersection.impl().boundaryId()
  ```

- The DGF block `general` is now always available.

- The DGFWriter will always write a boundary id and can write user-defined
  boundary data, now.

- The `MCMGMapper` can now be used to attach multiple dofs to each
  entity:
  the Layout is passed into the constructor and
  returns the number of dofs to attach to the given geometry type
  ```
     MCMGLayout layout = [](GeometryType gt, int griddim) {
       return gt.dim() == griddim? 2:0;
     };
     MCMGMapper mapper(grid,layout);
  ```
- The new method `MCMGMapper::indices(entity)` returns an iterable range
  (instance of `IntegralRange<Index>`)
  with the indices of dofs attached to the given entity:
  ```
    for (const auto& i : mapper.indices(entity) )
      dof = vector[i];
  ```

  [dune-grid!177]: https://gitlab.dune-project.org/core/dune-grid/merge_requests/177

- Two new method were added to the MCMGMapper:
  `size_type size(GeometryType)` and
  `const std::vector< GeometryType >& types (int codim)`
  returning the number of dofs attached to the geometry type and a vector
  with all geometry types on which dofs are attached, respectively.

- The `StructuredGridFactory` now returns a `unique_ptr` instead of a
  `shared_ptr` ([dune-grid!212][]).  Code that relies on a `shared_ptr`
  needs to explicitly assign the return value to a `shared_ptr`
  variable.

  [dune-grid!212]: https://gitlab.dune-project.org/core/dune-grid/merge_requests/212

- `SubsamplingVTKWriter` now supports arbitrary refinements, not just powers
  of two.  The old constructor taking a parameter `int levels` has been
  deprecated, you should now pass a parameter `RefinementIntervals intervals`
  instead.  There are convenience functions `refinementIntervals(int
  intervals)` and `refinementLevels(int levels)` to construct parameters of
  type `RefinementIntervals` in dune-geometry.

  [dune-grid!193]: https://gitlab.dune-project.org/core/dune-grid/merge_requests/193
