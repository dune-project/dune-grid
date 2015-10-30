# Add all tests without deprecations
dune_add_test(NAME test-yaspgrid-backuprestore-equidistant-without-deprecations
              SOURCES test-yaspgrid-backuprestore-equidistant.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

dune_add_test(NAME test-yaspgrid-backuprestore-equidistantoffset-without-deprecations
              SOURCES test-yaspgrid-backuprestore-equidistantoffset.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

dune_add_test(NAME test-yaspgrid-backuprestore-tensor-without-deprecations
              SOURCES test-yaspgrid-backuprestore-tensor.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

dune_add_test(NAME test-yaspgrid-tensorgridfactory-without-deprecations
              SOURCES test-yaspgrid-tensorgridfactory.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

dune_add_test(NAME test-yaspgrid-yaspfactory-1d-without-deprecations
              SOURCES test-yaspgrid-yaspfactory-1d.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

dune_add_test(NAME test-yaspgrid-yaspfactory-2d-without-deprecations
              SOURCES test-yaspgrid-yaspfactory-2d.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

dune_add_test(NAME test-yaspgrid-yaspfactory-3d-without-deprecations
              SOURCES test-yaspgrid-yaspfactory-3d.cc
              MPI_RANKS 1 2
              COMPILE_DEFINITIONS DISABLE_DEPRECATED_METHOD_CHECK=1
              COMPILE_FLAGS -Wall -Werror)

# Add all tests with deprecations!
dune_add_test(NAME test-yaspgrid-backuprestore-equidistant-with-deprecations
              SOURCES test-yaspgrid-backuprestore-equidistant.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)

dune_add_test(NAME test-yaspgrid-backuprestore-equidistantoffset-with-deprecations
              SOURCES test-yaspgrid-backuprestore-equidistantoffset.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)

dune_add_test(NAME test-yaspgrid-backuprestore-tensor-with-deprecations
              SOURCES test-yaspgrid-backuprestore-tensor.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)

dune_add_test(NAME test-yaspgrid-tensorgridfactory-with-deprecations
              SOURCES test-yaspgrid-tensorgridfactory.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)

dune_add_test(NAME test-yaspgrid-yaspfactory-1d-with-deprecations
              SOURCES test-yaspgrid-yaspfactory-1d.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)

dune_add_test(NAME test-yaspgrid-yaspfactory-2d-with-deprecations
              SOURCES test-yaspgrid-yaspfactory-2d.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)

dune_add_test(NAME test-yaspgrid-yaspfactory-3d-with-deprecations
              SOURCES test-yaspgrid-yaspfactory-3d.cc
              COMPILE_FLAGS -Wall -Werror -Wno-deprecated-declarations
              COMPILE_ONLY)
