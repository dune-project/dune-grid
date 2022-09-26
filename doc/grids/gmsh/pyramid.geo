// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// mesh width associated with points
lc = 1;

// vertices of the pyramid
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0.5, 0.5, 1, lc};

// lines of the pyramid
Line(1) = {1,2};
Line(2) = {3,2};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {2,5};
Line(7) = {3,5};
Line(8) = {4,5};

// faces
Line Loop(100) = {4,1,-2,3};  
Plane Surface(200) = {100};  
Physical Surface(300) = {200};

Line Loop(101) = {5,-8,4};    
Plane Surface(201) = {101};  Physical 
Surface(301) = {201};

Line Loop(102) = {1,6,-5};    
Plane Surface(202) = {102};  
Physical Surface(302) = {202};

Line Loop(103) = {-2,7,-6};   
Plane Surface(203) = {103};  
Physical Surface(303) = {203};

Line Loop(104) = {3,8,-7};    
Plane Surface(204) = {104};  
Physical Surface(304) = {204};

// volume
Surface Loop(400) = {200,-201,-202,-203,-204};
Volume(1000) = {400}; 
Physical Volume(1001) = {1000};