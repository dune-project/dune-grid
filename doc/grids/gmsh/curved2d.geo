// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// mesh width associated with points
lc  = 0.4;
lc2 = 0.05;

// points
Point(1) = {-1, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 0, 0, lc};

Point(4) = {2,0.5,0,lc};
Point(5) = {2.5,3,0,lc};
Point(6) = {-1,1,0,lc};

Point(7) = {1.3, 1.5, 0, lc2};
Point(8) = {1.5, 1.5, 0, lc2};
Point(9) = {1.7, 1.5, 0, lc2};

Point(10) = {-0.3, 0.5, 0, lc2};
Point(11) = { 0.1, 0.5, 0, lc2};
Point(12) = { 0.5, 0.5, 0, lc2};

// lines
Circle(1) = {1,3,2};
BSpline(2) = {2,4,5,6,1};

Circle(3) = {7,8,9};
Circle(4) = {9,8,7};

Circle(5) = {10,11,12};
Circle(6) = {12,11,10};

Line Loop(100) = {1,2};  
Line Loop(101) = {3,4};
Line Loop(102) = {5,6};

// surfaces
Plane Surface(200) = {100,101,102};  
