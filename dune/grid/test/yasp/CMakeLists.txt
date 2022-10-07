# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_add_test(NAME test-yaspgrid-backuprestore-equidistant
              SOURCES test-yaspgrid-backuprestore-equidistant.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(NAME test-yaspgrid-backuprestore-equidistantoffset
              SOURCES test-yaspgrid-backuprestore-equidistantoffset.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(NAME test-yaspgrid-backuprestore-tensor
              SOURCES test-yaspgrid-backuprestore-tensor.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(SOURCES test-yaspgrid-entityshifttable.cc)

dune_add_test(NAME test-yaspgrid-tensorgridfactory
              SOURCES test-yaspgrid-tensorgridfactory.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(NAME test-yaspgrid-yaspfactory-1d
              SOURCES test-yaspgrid-yaspfactory-1d.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(NAME test-yaspgrid-yaspfactory-2d
              SOURCES test-yaspgrid-yaspfactory-2d.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(NAME test-yaspgrid-yaspfactory-3d
              SOURCES test-yaspgrid-yaspfactory-3d.cc
              MPI_RANKS 1 2
              TIMEOUT 666
              )

dune_add_test(SOURCES test-yaspgrid-constructor.cc)
