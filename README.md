<!--
SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
-->

DUNE-library
============

DUNE, the Distributed and Unified Numerics Environment is a modular toolbox
for solving partial differential equations with grid-based methods.

The main intention is to create slim interfaces allowing an efficient use of
legacy and/or new libraries. Using C++ techniques DUNE allows to use very
different implementation of the same concept (i.e. grid, solver, ...) under
a common interface with a very low overhead.

DUNE was designed with flexibility in mind. It supports easy discretization
using methods, like Finite Elements, Finite Volume and also Finite
Differences. Through separation of data structures DUNE allows fast Linear
Algebra like provided in the ISTL module, or usage of external libraries
like blas.

This package contains the basic DUNE grid classes.

More information
----------------

Check dune-common for more details concerning dependencies, known bugs,
license and installation.
