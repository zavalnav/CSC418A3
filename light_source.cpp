/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <stdlib.h>
#include <cstdio>
#include "light_source.h"
#define TEXTURE_MAPPING 0
Point3D PointLight::get_random_sample()
{
	for (; ; )
	{
		double x = (double)rand() / RAND_MAX * _radius * 2 - _radius;
		double y = (double)rand() / RAND_MAX * _radius * 2 - _radius;
		double z = (double)rand() / RAND_MAX * _radius * 2 - _radius;
		if (x * x + y * y + z * z < _radius)
			return Point3D(x + _pos[0], y + _pos[1], z + _pos[2]);
	}
}

void PointLight::ambient( Ray3D& ray ) {
	// returns the ambient part of phong shading

	ray.col = ray.intersection.mat->ambient * _col_ambient;
}

void PointLight::shade( Ray3D& ray) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	Colour kd = ray.intersection.mat->diffuse;

	if (TEXTURE_MAPPING && ray.intersection.texture_enabled)
	{
		int row, col;
        //get texel coords
        SphereMapping sphereMap;
        sphereMap.getTexel(&row, &col, ray.intersection.texture.getHeight(), ray.intersection.texture.getWidth(), ray.intersection.pointOS);
            //change diffuse coefficient to that of texture
            unsigned char * red = ray.intersection.texture.getRed(row,col);
            unsigned char * green = ray.intersection.texture.getGreen(row,col);
            unsigned char * blue = ray.intersection.texture.getBlue(row,col);
            //get values for the color components in proper type
            double RedValue = (*red)/255.0;
            double BlueValue = (*blue)/255.0;
            double GreenValue = (*green)/255.0;
            kd = ray.intersection.mat->diffuse*Colour(RedValue, GreenValue, BlueValue); 
	}



	Vector3D s = _pos - ray.intersection.point; // light direction
	Vector3D n = ray.intersection.normal; // normal
	Vector3D d = -1.0 * ray.dir; // ray direction

	s.normalize();
	d.normalize();

	// ambient
	ray.col = ray.intersection.mat->ambient * _col_ambient;

	// diffuse
	// ray.col = ray.col + std::max(0.0, n.dot(s)) * ray.intersection.mat->diffuse * _col_diffuse;
	ray.col = ray.col + std::max(0.0, n.dot(s)) * kd * _col_diffuse;

	// specular
	Vector3D m = 2 * n.dot(s) * n - s;
	ray.col = ray.col + pow(std::max(0.0, m.dot(d)), ray.intersection.mat->specular_exp) * ray.intersection.mat->specular * _col_specular;

	ray.col.clamp();
}

