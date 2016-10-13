#ifndef Distance_H_
#define Distance_H_

#include "math.h"

double Distance_3D (double *p1, double *p2, double bx, double by, double bz)
{
	double delta_x; double delta_y; double delta_z;
	delta_x = fabs (p2[0] - p1[0]);
	delta_y = fabs (p2[1] - p1[1]);
	delta_z = fabs (p2[2] - p1[2]);
	if (fabs (delta_x) > 0.5 * bx)
		delta_x = delta_x - bx;
	if (fabs (delta_y) > 0.5 * by)
		delta_y = delta_y - by;
	if (fabs (delta_z) > 0.5 * bz)
		delta_z = delta_z - bz;

	return sqrt (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}

double Distance_3D_1 (double x1, double y1, double z1, double x2, double y2, double z2, double bx, double by, double bz)
{
	double delta_x; double delta_y; double delta_z;
	delta_x = fabs (x2 - x1);
	delta_y = fabs (y2 - y1);
	delta_z = fabs (z2 - z1);
	if (fabs (delta_x) > 0.5 * bx)
		delta_x = delta_x - bx;
	if (fabs (delta_y) > 0.5 * by)
		delta_y = delta_y - by;
	if (fabs (delta_z) > 0.5 * bz)
		delta_z = delta_z - bz;

	return sqrt (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}

#endif
