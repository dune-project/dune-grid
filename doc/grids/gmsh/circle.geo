// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// mesh width associated with points
lc = 1;

// Two antipodal points on the circle, the third point is the circle center
Point(1) = {-1, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 0, 0, lc};

Circle(1) = {1,3,2};
Circle(2) = {2,3,1};


Line Loop(100) = {1,2};  
Plane Surface(200) = {100};  
