# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

if(Alberta_FOUND)
  add_executable(test-alberta3d-refine test-alberta3d-refine.cc)
  target_link_libraries(test-alberta3d-refine PRIVATE dunegrid)
  add_dune_alberta_flags(test-alberta3d-refine WORLDDIM 3)
  dune_add_test(TARGET test-alberta3d-refine)
endif()
