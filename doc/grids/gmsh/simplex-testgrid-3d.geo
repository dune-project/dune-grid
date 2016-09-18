// Gmsh project created on Sun Sep 18 11:55:32 2016
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1.5, 0.8, 0, 1.0};
Delete {
  Point{4};
}
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, 0, 1, 1.0};
Point(6) = {1, 0, 1, 1.0};
Point(7) = {1, 1, 1, 1.0};
Point(8) = {0, 1, 1, 1.0};
Line(1) = {5, 6};
Line(2) = {6, 7};
Line(3) = {7, 8};
Line(4) = {8, 5};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {5, 1};
Line(10) = {6, 2};
Line(11) = {7, 3};
Line(12) = {8, 4};
Line Loop(13) = {4, 1, 2, 3};
Plane Surface(14) = {13};
Line Loop(15) = {8, 5, 6, 7};
Plane Surface(16) = {15};
Line Loop(17) = {2, 11, -6, -10};
Plane Surface(18) = {17};
Line Loop(19) = {4, 9, -8, -12};
Plane Surface(20) = {19};
Line Loop(21) = {12, -7, -11, 3};
Plane Surface(22) = {21};
Line Loop(23) = {9, 5, -10, -1};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 20, 24, 16, 18, 22};
Volume(26) = {25};
