// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// mesh width associated with points
lc = 0.5;

// circle
r  = 1.0;
cx = 0.0;
cy = 0.0;
cz = 0.0;
pbase=0;
lbase=0;
sbase=0;

// points
Point(pbase+1) = {cx, cy, cz, lc};

Point(pbase+2) = {cx+r, cy, cz, lc};
Point(pbase+3) = {cx, cy+r, cz, lc};
Point(pbase+4) = {cx, cy, cz+r, lc};
Point(pbase+5) = {cx-r, cy, cz, lc};
Point(pbase+6) = {cx, cy-r, cz, lc};

Point(pbase+7) = {cx+r, cy, cz-r, lc};
Point(pbase+8) = {cx, cy+r, cz-r, lc};
Point(pbase+9) = {cx-r, cy, cz-r, lc};
Point(pbase+10)= {cx, cy-r, cz-r, lc};
Point(pbase+11)= {cx, cy, cz-r, lc};

// lines
Circle(lbase+1) = {2,1,3};
Circle(lbase+2) = {3,1,4};
Circle(lbase+3) = {4,1,2};
Circle(lbase+4) = {3,1,5};
Circle(lbase+5) = {5,1,4};
Circle(lbase+6) = {5,1,6};
Circle(lbase+7) = {6,1,4};
Circle(lbase+8) = {6,1,2};

Circle(lbase+9) = {8,11,9};
Circle(lbase+10)= {9,11,10};
Circle(lbase+11)= {10,11,7};
Circle(lbase+12)= {7,11,8};

Line(lbase+13)= {8,3};
Line(lbase+14)= {9,5};
Line(lbase+15)= {10,6};
Line(lbase+16)= {7,2};

Line Loop(lbase+101) = {1,2,3};
Line Loop(lbase+102) = {4,5,-2};
Line Loop(lbase+103) = {6,7,-5};
Line Loop(lbase+104) = {8,-3,-7};

Line Loop(lbase+105) = {12,13,-1,-16};
Line Loop(lbase+106) = {9,14,-4,-13};
Line Loop(lbase+107) = {10,15,-6,-14};
Line Loop(lbase+108) = {11,16,-8,-15};

Line Loop(lbase+109) = {9,10,11,12};

// surfaces
Ruled Surface(sbase+1) = {101};  
Ruled Surface(sbase+2) = {102};  
Ruled Surface(sbase+3) = {103};  
Ruled Surface(sbase+4) = {104};  
Ruled Surface(sbase+5) = {105};  
Ruled Surface(sbase+6) = {106};  
Ruled Surface(sbase+7) = {107};  
Ruled Surface(sbase+8) = {108};  
Plane Surface(sbase+9) = {109};

Surface Loop(sbase+400) = {1,2,3,4,5,6,7,8,9};

// volume
Volume(1000) = {400}; 

