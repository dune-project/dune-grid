# master (will become 2.6)

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
