/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "light_source.h"

const double PI = acos(-1.0);

Point3D PointLight::get_random_sample()
{
	double a1 = (double)rand() / RAND_MAX * 2 * PI;
	double a2 = (double)rand() / RAND_MAX * 2 * PI;
	double r = (double)rand() / RAND_MAX * _radius;
	return Point3D(
		_pos[0] + r * sin(a1) * cos(a2),
		_pos[1] + r * sin(a1) * sin(a2),
		_pos[0] + r * cos(a1));
}

void PointLight::ambient( Ray3D& ray ) {
	// returns the ambient part of phong shading

	ray.col = ray.intersection.mat->ambient * _col_ambient;
}

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	Vector3D s = _pos - ray.intersection.point; // light direction
	Vector3D n = ray.intersection.normal; // normal
	Vector3D d = -1.0 * ray.dir; // ray direction

	s.normalize();
	d.normalize();

	printf("xxx\n");
	printf("%lf %lf %lf\n", n[0], n[1], n[2]);
	printf("%lf %lf %lf\n", ray.intersection.point[0], ray.intersection.point[1], ray.intersection.point[2]);
	printf("%lf\n", n.dot(s));
	printf("zzz\n");

	// ambient
	ray.col = ray.intersection.mat->ambient * _col_ambient;

	// diffuse
	ray.col = ray.col + std::max(0.0, n.dot(s)) * ray.intersection.mat->diffuse * _col_diffuse;

	// specular
	Vector3D m = 2 * n.dot(s) * n - s;
	ray.col = ray.col + pow(std::max(0.0, m.dot(d)), ray.intersection.mat->specular_exp) * ray.intersection.mat->specular * _col_specular;

	ray.col.clamp();
}

