#include <bits/stdc++.h>
using namespace std;
double genpoint()
{
	int a = rand()%100000;
	return a/100.0;
}
int main()
{
	srand(time(nullptr));
	int n;
	scanf("%d", &n);
	printf("%d %d %d %d\n", -10, 1010, -10, 1010);
	printf("%d\n", n);
	for (int i = 0; i < n; i++)
	{
		printf("%.2lf %.2lf\n", genpoint(), genpoint());
	}
}