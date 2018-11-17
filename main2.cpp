#include "voronoi.h"
#include <cstdio>
int n;
ld x[1000010], y[1000010];
int main()
{
	ld mnx, mxx, mny, mxy;
	scanf("%Lf%Lf%Lf%Lf", &mnx, &mxx, &mny, &mxy);
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
	{
		scanf("%Lf%Lf", &x[i], &y[i]);
	}
	auto diagram = makeVoronoi(n, x, y, mnx, mxx, mny, mxy);
	makeImage("image.jpg", diagram);
}