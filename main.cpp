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
	printf("There are %lu edges\n", diagram->edges.size());
	for (auto a : diagram->edges)
	{
		printf("(%.2Lf, %.2Lf) - (%.2Lf, %.2Lf)\n", a.first.x, a.first.y, a.second.x, a.second.y);
	}
}