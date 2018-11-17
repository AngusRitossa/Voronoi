/* Constructs an image based on a voronoi diagram */
#include "voronoi.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp> 
int WIDTH = 5000;
int HEIGHT = 5000;
/*int main()
{    
	// Create white empty image
	cv::Mat image = cv::Mat(1000, 1000, CV_8UC3, cv::Scalar(255,255,255));
	cv::line(image, cv::Point(0,0), cv::Point(511,511), cv::Scalar(255,0,0), 1);
	cv::imwrite("image.jpg", image);
}*/
std::pair<ld, ld> convertToPixels(ld x, ld y, ld width, ld height, VoronoiDiagram* v)
{
	ld a = (x - v->mnx)/(v->mxx - v->mnx);
	a *= width;
	ld b = (y - v->mny)/(v->mxy - v->mny);
	b *= height;
	return { a, b };
}
void makeImage(std::string imagename, VoronoiDiagram* v, bool colour)
{
	srand(time(nullptr));
	cv::Mat image = cv::Mat(WIDTH, HEIGHT, CV_8UC3, cv::Scalar(255,255,255));
	if (colour) // Colour the cells
	{
		for (auto cell : v->cells)
		{
			cell->sortPoints();
			cv::Point*
			ppt[1];
			ppt[0] = new cv::Point[cell->points.size()];
			for (unsigned int i = 0; i < cell->points.size(); i++)
			{
				std::pair<ld, ld> p = convertToPixels(cell->points[i].x, cell->points[i].y, WIDTH, HEIGHT, v);
				ppt[0][i] = cv::Point(p.first, p.second);
			}
			int npt[1];
			npt[0] = cell->points.size();
			const cv::Point* PPT[1] = { ppt[0] };
			cv::fillPoly(image, PPT, npt, 1, cv::Scalar(rand()%256, rand()%256, rand()%256), 8);
			delete[] ppt[0];
		}
	}
	// Draw points
	for (auto a : v->points)
	{
		std::pair<ld, ld> p = convertToPixels(a.x, a.y, WIDTH, HEIGHT, v);
		cv::circle(image, cv::Point(p.first, p.second), 20.0, cv::Scalar( 0, 0, 0 ), -1, 8);
	}
	// Draw lines
	for (auto a : v->edges)
	{
		std::pair<ld, ld> fir = convertToPixels(a.first.x, a.first.y, WIDTH, HEIGHT, v);
		std::pair<ld, ld> sec = convertToPixels(a.second.x, a.second.y, WIDTH, HEIGHT, v);
		cv::line(image, cv::Point(fir.first, fir.second), cv::Point(sec.first, sec.second), cv::Scalar(0,0,0), 4);
	}
	cv::imwrite(imagename, image);
}

