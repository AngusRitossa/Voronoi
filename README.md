# Voronoi

## General Information

Fortune's algorithm has been implemented in C++ to generate 2D Voronoi diagrams from a set of *n* points in O(*n* log *n*). The code can be used to generate a coloured image of the voronoi diagram using [OpenCV](https://opencv.org "OpenCV").

## Usage

### API

All of the API is defined in `voronoi.h`. 

#### Point

This class has two members, `long double x` and `long double y`, which define a point in the 2D plane.

#### VoronoiCell

This class represents a single cell in the Voronoi diagram. It stores two variables
- point, a `Point` which is the point that was used to make this Voronoi cell
- `points`, a vector of `Point` which are the points on the border of this Voronoi cell

The class also has one member function, `sortPoints()`, which sorts the vector `points` in their order on the boundary. 

#### VoronoiDiagram

This is the class that stores the completed Voronoi diagram. Stored in the class is 
- `mnx`, `mxx`, `mny` and `mxy` which define the bounding box
- `points`, a vector of `Point` which are the points which were used to generate the Voronoi diagram
- `edges`, a vector of pairs of `Point` which stores all the edges in the Voronoi diagram (excluding those on the bounding box)
- `cells`, a vector of `VoronoiCell*` which stores all the cells in the Voronoi diagram

#### makeVoronoi()

The function, `VoronoiDiagram* makeVoronoi(int n, long double *x, long double *y, ld mnx, ld mxx, ld mny, ld mxy)`, implemented in `voronoi.cpp`, is what creates the Voronoi diagram. The argument `n` is the number of points that will be used to make the Voronoi diagram, `x` and `y` are arrays of length `n` which store the points, and `mnx`, `mxx`, `mny` and `mxy` define the bounding box. The function returns a pointer to the `VoronoiDiagram` is creates.

#### makeImage()

The function, `void makeImage(std::string imagename, VoronoiDiagram* v, bool colour = true)`, implemented in `construct.cpp` is what creates the image of the Voronoi diagram. The first argument is the name of the image that this function will create, the second argument is the VoronoiDiagram to base the image off and the third argument is whether or not the cells should be coloured. The function creates an image called `imagename.jpg`. This requires [OpenCV](https://opencv.org "OpenCV").

### Sample Code

#### Input Format

Both sample main functions use the same input format.

- The first line contains 4 numbers which define the bounding box. Specifically, the minimum *x*-coordinate, the maximum *x*-coordinate, the minimum *y*-coordinate and the maximum *y*-coordinate
- The second line contains the integer *n* which is the number of points being used to make this Voronoi diagram
- The next *n* lines contain *x*<sub>*i*</sub> and *y*<sub>*i*</sub>, the points that are used to create the Voronoi diagram

[examples/in.txt](https://github.com/AngusRitossa/Voronoi/blob/master/examples/in.txt "examples/in.txt") contains an example of an input file. 

#### main.cpp

`main.cpp`, compiled with `compile.sh` provides a sample implementation of a main function which creates a Voronoi diagram and prints all the edges in it.

#### main2.cpp

`main2.cpp`, compiled with `compile2.sh` provides a sample implementation of a main function which creates a Voronoi diagram and then creates an image named `image.jpg` based on it. 