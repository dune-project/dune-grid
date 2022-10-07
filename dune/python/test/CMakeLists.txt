# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_python_add_test(NAME pytest_gf1
                     SCRIPT test_gf1.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pytest_gf2
                     SCRIPT test_gf2.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pytest_indexset
                     SCRIPT test_indexset.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pyinterpolate
                     SCRIPT interpolate.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pybackrest1
                     SCRIPT backrest1.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
if(dune-uggrid_FOUND)
dune_python_add_test(NAME pytestug
                     SCRIPT test-ug.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pyugfromfile
                     SCRIPT ugfromfile.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
endif()
