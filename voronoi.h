#include <cstdio>
#include <algorithm>
#include <utility>
#include <cmath>
#include <vector>
#include <queue>
#include <set>
#include <string>
#include <cassert>
typedef long double ld;
typedef class Point* ppoint;
typedef class GrowingEdge* pedge;
typedef class Edge* pfulledge;
typedef class VoronoiCell* pcell;
class Point // A point in the voronoi diagram
{
public:
	ld x, y; // The coordinates of the point
	Point(ld X, ld Y)
	{
		x = X;
		y = Y;
	}
	ld a(ld d);
	ld b(ld d);
	ld c(ld d);
	bool liesWithinBoundingBox(ld mnx, ld mxx, ld mny, ld mxy);
};

class GrowingEdge // Half an edge in the voronoi diagram; a ray
{
public:
	ld x, y; // The coordinates of the start of the ray
	ld vx, vy; // A vector in the direction of the ray
	pfulledge fullEdge; // The edge that this forms a part of
	// The line in the form ax + by + c = 0
	ld a();
	ld b();
	ld c();
	bool doesContainPoint(ld a, ld b); // Does (a, b) lie on the ray
};
class Edge // A full growing edge
{
public:
	pedge edge1 = nullptr, edge2 = nullptr; // The two growing edges that form this edge
	ppoint point1 = nullptr, point2 = nullptr; // The two ends of this edge
	pcell cell1, cell2; // The two voronoi cells that lie on either side of this edge
	void increaseToInfinity(); // If the edge only has one end, increase the other one 
	bool encaseInBoundingBox(ld mnx, ld mxx, ld mny, ld mxy); // If the edge only has one end, increase the other one 
	void updateBorderingCells(); // Pushes the points to the vectors in the two VoronoiCells that border this edge
	// The line in the form ax + by + c = 0
	ld a();
	ld b();
	ld c();
	bool doesContainPoint(ld a, ld b); // Does (a, b) lie on the line segment
};
typedef class TreapNode* pnode;
class TreapNode
{	
public:
	ppoint focus; // The point that defines this parabola
	pnode left = nullptr, right = nullptr, par = nullptr; // Children in the treap
	pedge leftEdge = nullptr, rightEdge = nullptr; // The growing edge on either side of this
	pnode leftInBeachLine = nullptr, rightInBeachLine = nullptr; // In the actual beach line, left and right nodes 
	pcell cell; // The voronoi cell that contains the focus
	ppoint leftHyperbola() { return leftInBeachLine ? leftInBeachLine->focus : nullptr; };
	ppoint rightHyperbola() { return rightInBeachLine ? rightInBeachLine->focus : nullptr; };
	int sz = 1; // Number of elements in its subtree
	int leftSubtreeSize(); // Number of elements in the left subtree
	int prior; // Random priority assigned to the node
	// The range for which this hyperbola is on the beachline
	ld leftBoundary(ld d);
	ld rightBoundary(ld d);
	// Determine left and right hyperbola's as well as subtree size
	void updateDetails();
	// Do the growing edges on either size intersect, if so, where?
	std::pair<bool, ld> edgeIntersection();

	TreapNode() // Set random priority
	{
		prior = rand();
	}
};
std::pair<pedge, pedge> findGrowingEdges(ld startx, ld starty, ppoint a, ppoint b); // Find the growing edges between two points
class PqNode // An element in the priority queue
{
public:
	ld y;
	bool siteEvent; // Is it a site event or edge-intersection event
	ld compy() const { return siteEvent ? y : y + 1e-3; }
	ppoint point; // If its a site event, which point are we inserting
	pnode node; // If its an edge-intersection event, what node in the treap is it
	bool operator<(const PqNode &o) const
	{
		return compy() < o.compy();
	}
	bool operator>(const PqNode &o) const
	{
		return compy() > o.compy();
	}
};
// Represents a single voronoi cell defined by one point
class VoronoiCell
{
public:
	Point point; // The 'centre' of the cell
	std::vector<Point> points; // The edges that make up this VoronoiCell
	VoronoiCell(ppoint p) : point (*p) {}
	void sortPoints(); // Sorts the points around the origin point
};
bool comp(Point a, Point b); // Computes whether a point is anticlockwise to another relative to the centre of the cell
// Return type
class VoronoiDiagram // Class that stores the completed voronoi diagram
{
public:
	std::vector<std::pair<Point, Point> > edges; // Stores all the edges in the voronoi diagram
	std::vector<VoronoiCell*> cells; // All the voronoi cells in this voronoi diagram
	std::vector<Point> points; // The points that were used to make the diagram
	ld mnx, mxx, mny, mxy;
	VoronoiDiagram(ld MNX, ld MXX, ld MNY, ld MXY)
	{
		mnx = MNX, mxx = MXX, mny = MNY, mxy = MXY;
	}
	void addPoint(ld x, ld y); // Add this point as a corner to the appropriate cell (used for corners of the bounding box)
};
VoronoiDiagram* makeVoronoi(int n, ld *x, ld *y, ld mnx, ld mxx, ld mny, ld mxy);
bool edgeIntersection(pedge edge1, pedge edge2); // Intersection of two complete edges
// Treap functions
void split(pnode a, pnode &left, pnode &right, int am);
void actualsplit(pnode a, pnode &left, pnode &right, int am);
void merge(pnode &a, pnode left, pnode right);
void actualmerge(pnode &a, pnode left, pnode right);
// Creating the image
void makeImage(std::string imagename, VoronoiDiagram* v, bool colour = true);
bool cross(ld x, ld y, ld a, ld b); // Sign of the cross product of two vectors