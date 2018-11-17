#include "voronoi.h"
#ifdef DEBUG
	#define D(x...) printf(x)
#else 
	#define D(x...) 
#endif
/* Functions to generate a quadratic equation given a focus point and directrix */
ld Point::a(ld d)
{
	return 1/(2*(y-d));
}
ld Point::b(ld d)
{
	return (2*x)/(2*(d-y));
}
ld Point::c(ld d)
{
	return (x*x + y*y - d*d)/(2*(y-d));
}
bool Point::liesWithinBoundingBox(ld mnx, ld mxx, ld mny, ld mxy)
{
	return mnx <= x && x <= mxx && mny <= y && y <= mxy;
}
/* Functions to generate a the equation for the line */
ld GrowingEdge::a() // y2 - y1
{
	return vy;
}
ld GrowingEdge::b() // x1 - x2 
{
	return -vx;
}
ld GrowingEdge::c() // -Ax1 - By1
{
	return -a()*x - b()*y;
}
ld Edge::a() // y2 - y1
{
	return point2->y - point1->y;
}
ld Edge::b() // x1 - x2 
{
	return point1->x - point2->x;
}
ld Edge::c() // -Ax1 - By1
{
	return -a()*point1->x - b()*point1->y;
}
/* Does a point lie on a growing edge (the ray) 
   Assumes that the point lies on the line formed by extending the ray both directions */
bool GrowingEdge::doesContainPoint(ld a, ld b)
{
	if (std::fabs(vy) < 1e-7) // Line is horizontal, compare x
	{
		if (vx > 0) return a + 1e-7 >= x;
		else return a - 1e-7 <= x;
	}
	// Otherwise, compare y
	if (vy > 0) return b + 1e-7 >= y;
	else return b - 1e-7 <= y;
}
bool Edge::doesContainPoint(ld a, ld b)
{
	if (point1->x != point2->x && point1->y != point2->y)
		return std::min(point1->x, point2->x) <= a && a <= std::max(point1->x, point2->x)
		&& std::min(point1->y, point2->y) <= b && b <= std::max(point1->y, point2->y);
	else if (point1->x != point2->x)
		return std::min(point1->x, point2->x) <= a && a <= std::max(point1->x, point2->x);
	else
		return std::min(point1->y, point2->y) <= b && b <= std::max(point1->y, point2->y);
}
/* Find the growing edges from two points */
std::pair<pedge, pedge> findGrowingEdges(ld startx, ld starty, ppoint a, ppoint b)
{
	pedge left = new GrowingEdge();
	pedge right = new GrowingEdge();
	left->x = right->x = startx;
	left->y = right->y = starty;
	// Calculate vectors in the correct direction
	left->vx = a->y - b->y;
	left->vy = b->x - a->x;
	right->vx = -left->vx;
	right->vy = -left->vy;
	if (left->vx > right->vx) std::swap(left, right);
	return { left, right };
}  
/* Functions to find the left and right boundaries of a hyperbola */
std::pair<ld, ld> quadraticSolution(ld a, ld b, ld c) // Returns the two real roots of a quadratic, assumes they exist
{
	// Special case time: a == 0
	if (fabs(a) < 1e-6)
	{
		return { -c/b, -c/b };
	}
	ld d = b*b - 4*a*c;
	assert(d >= 0);
	d = sqrt(d);
	d /= 2*a;
	if (d < 0) d = -d;
	ld ans = -b/(2*a);
	return { ans-d, ans+d };
}
ld TreapNode::leftBoundary(ld d)
{
	if (!leftHyperbola()) return -1e18; // Is the leftmost, return negative infinity
	// Special case time: vertical lines
	if (fabs(d-focus->y) < 1e-6) { return focus->x; }
	if (fabs(d-leftHyperbola()->y) < 1e-6) { return leftHyperbola()->x; }
	std::pair<ld, ld> inters = quadraticSolution(focus->a(d) - leftHyperbola()->a(d), focus->b(d) - leftHyperbola()->b(d), focus->c(d) - leftHyperbola()->c(d));
	// If we are lower than our left neighbour, take the left (lower) intersection, otherwise its the right one
	if (focus->y < leftHyperbola()->y) return inters.first;
	else return inters.second;
}
ld TreapNode::rightBoundary(ld d)
{
	if (!rightHyperbola()) return 1e18; // Is the rightmost, return infinity
	// Special case time: vertical lines
	if (fabs(d-focus->y) < 1e-6) { return focus->x; }
	if (fabs(d-rightHyperbola()->y) < 1e-6) { return rightHyperbola()->x; }
	std::pair<ld, ld> inters = quadraticSolution(focus->a(d) - rightHyperbola()->a(d), focus->b(d) - rightHyperbola()->b(d), focus->c(d) - rightHyperbola()->c(d));
	// If we are lower than our right neighbour, take the right (higher) intersection, otherwise its the left one
	if (focus->y < rightHyperbola()->y) return inters.second;
	else return inters.first;
}
void TreapNode::updateDetails()
{
	sz = 1; // Size of the subtree
	if (left) sz += left->sz, left->par = this;
	if (right) sz += right->sz, right->par = this;
}
int TreapNode::leftSubtreeSize()
{
	if (!left) return 0;
	return left->sz;
}
/* The intersection of the two growing edges around a point */
ld xCoordOfEdgeIntersection; // Stores the x coordinate of the last edge intersection
ld yCoordOfEdgeIntersection; // Stores the y coordinate of the last edge intersection
std::pair<bool, ld> TreapNode::edgeIntersection()
{
	if (!leftEdge || !rightEdge) return { 0, 0.0 }; // No intersection because at least one doesn't exist
	// Find the point of intersection
	ld det = leftEdge->a()*rightEdge->b() - rightEdge->a()*leftEdge->b();
	if (std::fabs(det) < 1e-7) return { 0, 0.0 }; // No intersection because they are parallel
	ld x = (leftEdge->b()*rightEdge->c() - rightEdge->b()*leftEdge->c()) / det;
	ld y = (rightEdge->a()*leftEdge->c() - leftEdge->a()*rightEdge->c()) / det;
	if (!leftEdge->doesContainPoint(x, y) || !rightEdge->doesContainPoint(x, y)) return { 0, 0.0 }; // This point does not lie on both rays
	// Find the height of the sweepline where they would intersect on the beachline
	xCoordOfEdgeIntersection = x; // Store globally so it can be accessed later
	yCoordOfEdgeIntersection = y;
	y -= hypot(x-focus->x, y-focus->y);
	return { 1, y };
}
/* Intersection of two complete edges */
bool edgeIntersection(pfulledge edge1, pfulledge edge2)
{
	if (!edge1 || !edge2) return 0; // No intersection because at least one doesn't exist
	// Find the point of intersection
	ld det = edge1->a()*edge2->b() - edge2->a()*edge1->b();
	if (fabs(det) < 1e-7) return 0; // No intersection because they are parallel
	ld x = (edge1->b()*edge2->c() - edge2->b()*edge1->c()) / det;
	ld y = (edge2->a()*edge1->c() - edge1->a()*edge2->c()) / det;
	if (!edge1->doesContainPoint(x, y) || !edge2->doesContainPoint(x, y)) return 0; // This point does not lie on both rays
	// Find the height of the sweepline where they would intersect on the beachline
	xCoordOfEdgeIntersection = x; // Store globally so it can be accessed later
	yCoordOfEdgeIntersection = y;
	return 1;
}
/* Treap Code */
pnode leftmost(pnode a) // Returns the leftmost node in the treap
{
	if (a->left) return leftmost(a->left);
	else return a;
}
pnode rightmost(pnode a) // Returns the rightmost node in the treap
{
	if (a->right) return rightmost(a->right);
	else return a;
}
void split(pnode a, pnode &left, pnode &right, int am)
{
	actualsplit(a, left, right, am); // Perform the split
	// Update leftInBeachLine and rightInBeachLine for the bordering nodes
	if (left && right)
	{
		pnode r = rightmost(left);
		pnode l = leftmost(right);
		r->rightInBeachLine = l->leftInBeachLine = nullptr;
	}
}
void actualsplit(pnode a, pnode &left, pnode &right, int am) // Splits the treap so that left contains the lowest am nodes and right has the rest
{
	if (!a)
	{
		left = right = nullptr;
		return;
	}
	if (a->leftSubtreeSize()+1 <= am) // Put a in the left subtree
	{
		am -= a->leftSubtreeSize()+1;
		actualsplit(a->right, a->right, right, am);
		left = a;
		left->updateDetails();
	}
	else // Put a in the right subtree
	{
		actualsplit(a->left, left, a->left, am);
		right = a;
		right->updateDetails();
	}
}
void merge(pnode &a, pnode left, pnode right)
{
	// Update leftInBeachLine and rightInBeachLine for the bordering nodes
	if (left && right)
	{
		pnode r = rightmost(left);
		pnode l = leftmost(right);
		r->rightInBeachLine = l;
		l->leftInBeachLine = r;
	}
	actualmerge(a, left, right); // Perform the split
}
void actualmerge(pnode &a, pnode left, pnode right) // Merges left and right into one treap
{
	if (!left || !right)
	{
		a = left ? left : right;
		return;
	}
	if (left->prior > right->prior) // Put left above right
	{
		actualmerge(left->right, left->right, right);
		a = left;
	}
	else
	{
		actualmerge(right->left, left, right->left);
		a = right;
	}
	a->updateDetails();
}
std::pair<bool, int> findOnBeachLine(pnode a, ld x, ld d) // Returns the 1-indexed location of the parabola which contains the x-coordinate at height d
{
	ld l = a->leftBoundary(d);
	ld r = a->rightBoundary(d);
	D("l: %.2Lf r: %.2Lf - target %.2Lf - point (%.2Lf, %.2Lf) - d %.2Lf\n", l, r, x, a->focus->x, a->focus->y, d);
	if (fabs(x-l) < 1e-6) return { 0, a->leftSubtreeSize() }; // Is between two, on the left
	if (fabs(x-r) < 1e-6) return { 0, a->leftSubtreeSize() + 1 }; // Is between two, on the right
	if (l <= x && x <= r) return { 1, a->leftSubtreeSize() + 1 };
	else if (x < l) // Is in the left subtree
	{
		return findOnBeachLine(a->left, x, d);
	}
	else // Is in right subtree
	{
		std::pair<bool, int> ans = findOnBeachLine(a->right, x, d);
		ans.second += a->leftSubtreeSize() + 1;
		return ans;
	}
}
int findIndexOnBeachLine(pnode a) // Returns the 1-indexed location of this node
{
	int am = a->leftSubtreeSize()+1;
	while (a->par)
	{
		pnode par = a->par;
		if (par->right == a) am += par->leftSubtreeSize()+1;
		a = par;
	}
	return am;
}
int findAmongEqualHeight(pnode a, ld x) // Finds the location in the beach line when it is ordered by x-coord of focus
{
	if (!a) return 0;
	if (a->focus->x > x) return findAmongEqualHeight(a->left, x);
	else return findAmongEqualHeight(a->right, x) + a->leftSubtreeSize()+1;
}
/* Creates the edge intersection event (if applicable) and pushes to the pq */
void pushEdgeIntersectionEvent(std::priority_queue<PqNode> &pq, pnode node)
{
	// Find the intersection
	std::pair<bool, ld> inter = node->edgeIntersection();
	if (!inter.first) return; // There is no intersection
	PqNode* event = new PqNode(); // Create the new event
	event->y = inter.second;
	event->siteEvent = 0;
	event->node = node;
	D("new event\n");
	pq.push(*event);
}
/* Processing edges */
std::pair<ld, ld> extendRayToInfinity(Point p, pedge e)
{
	ld x, y;
	if (fabs(e->vx) > fabs(e->vy)) // Increase vx
	{
		ld am;
		if (e->vx > 0) am = 1e10 - p.x;
		else am = -1e10 + p.x;
		ld multam = am / e->vx;
		x = p.x + multam*e->vx; // = 1e18
		y = p.y + multam*e->vy;
	}
	else
	{
		ld am;
		if (e->vy > 0) am = 1e10 - p.y;
		else am = -1e10 + p.y;
		ld multam = am / e->vy; 
		y = p.y + multam*e->vy; // = 1e18
		x = p.x + multam*e->vx;
	}
	return { x, y };
}
void Edge::increaseToInfinity()
{
	if (point1 && point2) return; // Aleady has two ends
	if (!point1 && !point2) 
	{
		// Extend both rays individually
		std::pair<ld, ld> a, b;
		a = extendRayToInfinity(Point(edge1->x, edge1->y), edge1);
		b = extendRayToInfinity(Point(edge2->x, edge2->y), edge2);
		// Then get the end of each to make the line
		point1 = new Point(a.first, a.second);
		point2 = new Point(b.first, b.second);
		return;
	}
	if (point2) std::swap(point1, point2), std::swap(edge1, edge2);
	point2 = new Point(0, 0);
	std::pair coords = extendRayToInfinity(*point1, edge2);
	point2->x = coords.first;
	point2->y = coords.second;
}
bool Edge::encaseInBoundingBox(ld mnx, ld mxx, ld mny, ld mxy)
{
	increaseToInfinity();
	bool point1within = point1->liesWithinBoundingBox(mnx, mxx, mny, mxy);
	bool point2within = point2->liesWithinBoundingBox(mnx, mxx, mny, mxy);
	if (point1within && point2within) return 1;
	if (point2within) std::swap(point1, point2), std::swap(edge1, edge2), std::swap(point1within, point2within);
	if (point1within) // One point lies in the bounding box
	{
		// Look for one intersect
		// Top edge
		pfulledge side = new Edge();
		side->point1 = new Point(0, 0);
		side->point2 = new Point(0, 0);
		side->point1->x = mxx, side->point2->x = mxx, side->point1->y = mny, side->point2->y = mxy;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			point2->x = xCoordOfEdgeIntersection, point2->y = yCoordOfEdgeIntersection; return 1;
		}
		// Bottom edge
		side->point1->x = mnx, side->point2->x = mnx, side->point1->y = mny, side->point2->y = mxy;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			point2->x = xCoordOfEdgeIntersection, point2->y = yCoordOfEdgeIntersection; return 1;
		}
		// Left edge
		side->point1->x = mnx, side->point2->x = mxx, side->point1->y = mny, side->point2->y = mny;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			point2->x = xCoordOfEdgeIntersection, point2->y = yCoordOfEdgeIntersection; return 1;
		}
		// Right edge
		side->point1->x = mnx, side->point2->x = mxx, side->point1->y = mxy, side->point2->y = mxy;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			point2->x = xCoordOfEdgeIntersection, point2->y = yCoordOfEdgeIntersection; return 1;
		}
		printf("Failed to find intersection - Point (%.2Lf, %.2Lf) - Edge (%.2Lf, %.2Lf)\n", point1->x, point1->y, point2->x, point2->y);
		assert(0);
		return 0;
	}
	else
	{
		// Look for exactly two intersections
		std::vector<std::pair<ld, ld> > intersections;
		// Top edge
		pfulledge side = new Edge();
		side->point1 = new Point(0, 0);
		side->point2 = new Point(0, 0);
		side->point1->x = mxx, side->point2->x = mxx, side->point1->y = mny, side->point2->y = mxy;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			intersections.emplace_back(xCoordOfEdgeIntersection, yCoordOfEdgeIntersection);
		}
		// Bottom edge
		side->point1->x = mnx, side->point2->x = mnx, side->point1->y = mny, side->point2->y = mxy;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			intersections.emplace_back(xCoordOfEdgeIntersection, yCoordOfEdgeIntersection);
		}
		// Left edge
		side->point1->x = mnx, side->point2->x = mxx, side->point1->y = mny, side->point2->y = mny;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			intersections.emplace_back(xCoordOfEdgeIntersection, yCoordOfEdgeIntersection);
		}
		// Right edge
		side->point1->x = mnx, side->point2->x = mxx, side->point1->y = mxy, side->point2->y = mxy;
		if (edgeIntersection(this, side)) // Found the end location
		{ 
			intersections.emplace_back(xCoordOfEdgeIntersection, yCoordOfEdgeIntersection);
		}
		std::sort(intersections.begin(), intersections.end());
		intersections.erase(std::unique(intersections.begin(), intersections.end()), intersections.end());
		assert(intersections.size() < 3);
		if (intersections.size() == 2)
		{
			point1->x = intersections[0].first;
			point1->y = intersections[0].second;
			point2->x = intersections[1].first;
			point2->y = intersections[1].second;
			return 1;
		}
		else return 0;
	}
}
void Edge::updateBorderingCells()
{
	cell1->points.push_back(*point1);
	cell1->points.push_back(*point2);
	cell2->points.push_back(*point1);
	cell2->points.push_back(*point2);
}
/* Sort all the points around a specific cell*/
Point sortingpoint(0, 0);
void VoronoiCell::sortPoints()
{
	// Split the points into two vectors based on whether they are left or right of the centre point
	std::vector<Point> left, right;
	for (auto a : points)
	{
		if (a.x < point.x) left.push_back(a);
		else right.push_back(a);
	}
	// Sort individually
	sortingpoint = point;
	std::sort(left.begin(), left.end(), comp);
	std::sort(right.begin(), right.end(), comp);
	// Push all the points back into the original vector
	points.clear();
	for (auto a : right) points.push_back(a);
	for (auto a : left) points.push_back(a);
}
bool comp(Point a, Point b) // Computes whether a point is anticlockwise to another relative to the centre of the cell
{
	return cross(a.x - sortingpoint.x, a.y - sortingpoint.y, b.x - sortingpoint.x, b.y - sortingpoint.y);
}
void printdfs(pnode a)
{
	if (!a) return;
	printdfs(a->left);
	D("(%.2Lf, %.2Lf), ", a->focus->x, a->focus->y);
	printdfs(a->right);
}
void printBeachLine(pnode a)
{
	D("Beach line:\n");
	printdfs(a);
	D("\n");
}
bool cross(ld x, ld y, ld a, ld b) // Sign of the cross product
{
	return x*b - y*a > 0;
}
void updateFullEdge(pedge a, ld x, ld y) // Ends the growing edge at the coordinates (x, y)
{
	if (a->fullEdge->edge1 == a)
	{
		a->fullEdge->point1 = new Point(x, y);
	}
	else
	{
		a->fullEdge->point2 = new Point(x, y);
	}
}
void sameHeightCreateEdges(pnode a, pnode b, ld mxhei, std::vector<pfulledge> &edges) // For two nodes at the same height (the highest nodes), create the growing edges between them 
{
	// Find the growing edges, one will be pointing directly upwards, one downwards
	pedge growingEdge = new GrowingEdge();
	growingEdge->x = (a->focus->x + b->focus->x)/2;
	growingEdge->y = 1e10;
	growingEdge->vx = 0;
	growingEdge->vy = -1;
	// Create the full edge
	pfulledge fullEdge = new Edge();
	fullEdge->edge1 = growingEdge;
	fullEdge->point2 = new Point(growingEdge->x, 1e10);
	fullEdge->cell1 = a->cell;
	fullEdge->cell2 = b->cell;
	growingEdge->fullEdge = fullEdge;
	edges.push_back(fullEdge);
	// Set their left and right edges to the downwards edge
	a->rightEdge = b->leftEdge = growingEdge;
}
void VoronoiDiagram::addPoint(ld x, ld y) // Adds this point as a corner to the appropriate cell (used for corners of the bounding box)
{
	pcell closest = nullptr;
	ld dis = 1e18;
	for (auto cell : cells)
	{
		ld d = hypot(cell->point.x - x, cell->point.y - y);
		if (d < dis) dis = d, closest = cell;
	}
	closest->points.push_back(Point(x, y));
}
/* Code that runs Fortune's algorithm */
VoronoiDiagram* makeVoronoi(int n, ld *x, ld *y, ld mnx, ld mxx, ld mny, ld mxy)
{
	VoronoiDiagram* diagram = new VoronoiDiagram(mnx, mxx, mny, mxy);
	std::vector<pfulledge> edges;

	std::priority_queue<PqNode> pq; // The event queue
	for (int i = 0; i < n; i++)
	{
		PqNode* a = new PqNode();
		a->siteEvent = 1;
		a->point = new Point(x[i], y[i]);
		a->y = y[i];
		pq.push(*a);
		diagram->points.push_back(*a->point);
	}
	pnode treap = nullptr;
	ld mxhei;
	while (!pq.empty()) // Handling the event queue
	{
		PqNode a = pq.top();
		pq.pop();
		D("\nPopped off pq: %.2Lf %d\n", a.y, a.siteEvent);
		if (a.siteEvent)
		{
			if (!treap) // Is the first node
			{
				mxhei = a.y; // This is the highest value.
				treap = new TreapNode();
				treap->focus = a.point;
				treap->updateDetails();
				treap->cell = new VoronoiCell(a.point);
				diagram->cells.push_back(treap->cell);
				continue;
			}
			else if (mxhei-a.y < 1e-3) // This is the same height as the highest value (or close) - special case
			{
				D("Has the same height as the maximum\n");
				int loc = findAmongEqualHeight(treap, a.point->x);
				pnode newmid = new TreapNode();
				newmid->focus = a.point;
				newmid->cell = new VoronoiCell(a.point);
				diagram->cells.push_back(newmid->cell);
				pnode left, right;
				split(treap, left, right, loc);
				if (left) // Find the growing edges
				{
					sameHeightCreateEdges(rightmost(left), newmid, mxhei, edges);
				}
				if (right) // Find the growing edges
				{
					sameHeightCreateEdges(newmid, leftmost(right), mxhei, edges);
				}
				// Insert newmid into the treap
				merge(treap, left, newmid);
				merge(treap, treap, right);
				continue;
			}
			printBeachLine(treap);
			// Find the arc on the beachline that this splits
			std::pair<bool, int> locOnBeachLine = findOnBeachLine(treap, a.point->x, a.y);
			int loc = locOnBeachLine.second;
			if (!locOnBeachLine.first) // Doesn't split one, instead goes between two
			{
				D("Middle of two beach line things\n");
				pnode left, right;
				split(treap, left, right, loc);
				pnode last = rightmost(left);
				pnode fir = leftmost(right);
				ld interx = a.point->x;
				ld intery = last->focus->a(a.y)*interx*interx + last->focus->b(a.y)*interx + last->focus->c(a.y);
				D("interx intery %Lf %Lf\n", interx, intery);
				// Get update the half edge that is between last & fir
				updateFullEdge(last->rightEdge, interx, intery);
				// Create new node
				pnode newmid = new TreapNode();
				newmid->focus = a.point;
				newmid->cell = new VoronoiCell(a.point); // Create the new cell for the point being processed
				diagram->cells.push_back(newmid->cell);
				newmid->leftEdge = last->rightEdge = findGrowingEdges(interx, intery, a.point, last->focus).first;
				newmid->rightEdge = fir->leftEdge = findGrowingEdges(interx, intery, a.point, fir->focus).second;
				// Create full edges
				pfulledge fullEdge = new Edge();
				fullEdge->edge2 = newmid->leftEdge;
				fullEdge->point1 = new Point(interx, intery);
				fullEdge->cell1 = newmid->cell;
				fullEdge->cell2 = last->cell;
				newmid->leftEdge->fullEdge = fullEdge;
				edges.push_back(fullEdge); // Save this edge to be processed later

				fullEdge = new Edge();
				fullEdge->edge2 = newmid->rightEdge;
				fullEdge->point1 = new Point(interx, intery);
				fullEdge->cell1 = newmid->cell;
				fullEdge->cell2 = fir->cell;
				newmid->rightEdge->fullEdge = fullEdge;
				edges.push_back(fullEdge); // Save this edge to be processed later

				pushEdgeIntersectionEvent(pq, fir);
				pushEdgeIntersectionEvent(pq, last);

				merge(treap, left, newmid);
				merge(treap, treap, right);
				continue;
			}
			D("Found intersection: %d out of %d\n", loc, treap->sz);
			// Split treap into (elements before loc), (loc), (elements after loc)
			pnode left, mid, right;
			split(treap, left, right, loc);
			split(left, left, mid, loc-1);
			assert(mid->sz == 1);
			// Find the coordinates where the new arc splits the old one by subsituting into quadratic equation
			ld interx = a.point->x;
			ld intery = mid->focus->a(a.y)*interx*interx + mid->focus->b(a.y)*interx + mid->focus->c(a.y);
			D("interx intery %Lf %Lf\n", interx, intery);
			std::pair<pedge, pedge> growingEdges = findGrowingEdges(interx, intery, a.point, mid->focus);
			pfulledge fullEdge = new Edge(); // The full edge that consists of these half edges
			fullEdge->edge1 = growingEdges.first;
			fullEdge->edge2 = growingEdges.second;
			growingEdges.first->fullEdge = growingEdges.second->fullEdge = fullEdge;
			fullEdge->cell1 = mid->cell; 
			fullEdge->cell2 = new VoronoiCell(a.point); // Create the new cell for the point being processed
			diagram->cells.push_back(fullEdge->cell2);
			edges.push_back(fullEdge); // Save this edge to be processed later

			// Create three new treap nodes to be added to the beachline (the outer two being the old mid, the inner being the new arc)
			pnode newleft = new TreapNode();
			pnode newmid = new TreapNode();
			pnode newright = new TreapNode();

			newleft->focus = newright->focus = mid->focus;
			newleft->cell = newright->cell = mid->cell;

			newmid->focus = a.point;
			newmid->cell = fullEdge->cell2;

			newleft->leftEdge = mid->leftEdge;
			newright->rightEdge = mid->rightEdge;
			newleft->rightEdge = newmid->leftEdge = growingEdges.first;
			newright->leftEdge = newmid->rightEdge = growingEdges.second;
			// Discount the old edge-intersection event and create three new ones
			mid->leftEdge = mid->rightEdge = nullptr;
			pushEdgeIntersectionEvent(pq, newleft);
			pushEdgeIntersectionEvent(pq, newmid);
			pushEdgeIntersectionEvent(pq, newright);
			// Merge the treaps back together
			merge(newleft, newleft, newmid);
			merge(newleft, newleft, newright);
			merge(treap, left, newleft);
			merge(treap, treap, right);
		}
		else
		{
			std::pair<bool, ld> inter = a.node->edgeIntersection();
			if (inter != std::make_pair(true, a.y)) continue; // This edge intersection no longer exists
			// Fix up the half edges, since they now have an end point
			updateFullEdge(a.node->leftEdge, xCoordOfEdgeIntersection, yCoordOfEdgeIntersection);
			updateFullEdge(a.node->rightEdge, xCoordOfEdgeIntersection, yCoordOfEdgeIntersection);
			// Find the nodes to the left and right of a
			pnode left = a.node->leftInBeachLine;
			pnode right = a.node->rightInBeachLine;
			D("Two points: (%.2Lf, %.2Lf), (%.2Lf, %.2Lf)\n", a.node->leftEdge->x, a.node->leftEdge->y, a.node->rightEdge->x, a.node->rightEdge->y);
			// Create the new growing edge
			ld startx = xCoordOfEdgeIntersection; // Just about the dodgiest 2 lines in this
			ld starty = yCoordOfEdgeIntersection;
			D("Node at (%Lf, %Lf)\n", startx, starty);
			std::pair<pedge, pedge> growingEdges = findGrowingEdges(startx, starty, left->focus, right->focus);
			// Find which of the two is pointing out of the new cell
			pedge downwardsEdge;
			if (cross(a.node->leftEdge->vx, a.node->leftEdge->vy, growingEdges.first->vx, growingEdges.first->vy)
			 == cross(a.node->leftEdge->vx, a.node->leftEdge->vy, a.node->rightEdge->vx, a.node->rightEdge->vy)) 
				downwardsEdge = growingEdges.first;
			else downwardsEdge = growingEdges.second;
			D("growing edge %Lf %Lf\n", downwardsEdge->vx, downwardsEdge->vy);
			// Create the fulledge for this growing edge
			pfulledge fullEdge = new Edge(); 
			fullEdge->edge2 = downwardsEdge;
			fullEdge->point1 = new Point(startx, starty);
			fullEdge->cell1 = left->cell; // Save the two cells that 
			fullEdge->cell2 = right->cell; // this edge is between
			downwardsEdge->fullEdge = fullEdge;
			edges.push_back(fullEdge); // Save this edge to be processed later
			// Update the growing edges & intersection points of left & right
			left->rightEdge = right->leftEdge = downwardsEdge;
			pushEdgeIntersectionEvent(pq, left);
			pushEdgeIntersectionEvent(pq, right);
			// Locate and remove a.node from the beachline
			int loc = findIndexOnBeachLine(a.node);
			pnode treapLeft, treapMid, treapRight;
			split(treap, treapLeft, treapRight, loc);
			split(treapLeft, treapLeft, treapMid, loc-1);
			assert(treapMid->sz == 1);
			assert(treapMid == a.node);
			// Remerge, without mid
			merge(treap, treapLeft, treapRight);
			a.node->leftEdge = a.node->rightEdge = nullptr;
		}	
	}
	// Add corners to the appropriate cells
	diagram->addPoint(mnx, mny);
	diagram->addPoint(mnx, mxy);
	diagram->addPoint(mxx, mny);
	diagram->addPoint(mxx, mxy);
	// Now, to deal with the bounding box
	for (auto a : edges) // Process all the edges
	{
		if (a->encaseInBoundingBox(mnx, mxx, mny, mxy)) diagram->edges.push_back({ *a->point1, *a->point2 }), a->updateBorderingCells();
		else D("Edge ignored\n");
	}
	return diagram;
}
