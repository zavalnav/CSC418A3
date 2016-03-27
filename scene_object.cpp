/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {

	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	//Transform ray into object space
	//Since z = 0 and normal = (0,0,1)
	//R(t) = (X0 + Xd*t, Y0 + Yd*t, Z0 + Zd*t) = (X0, Y0, Z0)
	//printf("ray orig in view: (%d, %d, %d) \n", ray.origin[0], ray.origin[1], ray.origin[2]);
	//printf("ray dir vecotor in view: (%d, %d, %d) \n", ray.dir[0], ray.dir[1], ray.dir[2]);

	Point3D originInObject = worldToModel * ray.origin; 
	Vector3D dirInObject = worldToModel * ray.dir;
	//printf("ray orig in world: (%d, %d, %d) \n", originInObject[0], originInObject[1], originInObject[2]);
	//printf("ray dir vecotor: (%d, %d, %d) \n", dirInObject[0], dirInObject[1], dirInObject[2]);

	Vector3D normal = Vector3D(0, 0, 1);
	Vector3D P1 = Vector3D(0.5, 0.5, 0);
	float dist = -originInObject[2] / dirInObject[2];
	Point3D intersect_point = originInObject + dist * dirInObject;
/*
	Point3D intersect_point = Point3D(
		(originInObject[0] - 0.5)/dirInObject[0], 
		(originInObject[1] - 0.5)/dirInObject[1], 
		originInObject[2]/dirInObject[2]);
*/
	if (intersect_point[0] >= -0.5 && intersect_point[0] <= 0.5 &&
		intersect_point[1] >= -0.5 && intersect_point[1] <= 0.5)
	{
		Point3D new_intersection = modelToWorld * intersect_point;
		float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
		if (new_t_value < ray.intersection.t_value)
		{
			ray.intersection.none = false;
			ray.intersection.t_value = new_t_value;
			ray.intersection.point = new_intersection;
			ray.intersection.normal = transNorm(worldToModel, normal);
			return true;
		}
		else
			return false;
		//ray.intersection.point = originInObject;
		//ray.intersection.normal = normal;
		//ray.intersection.none = false;
	}
	else
		return false;

	//Find intesection with object. Our ray = dirInObject

	//return !ray.intersection.none;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	// The intersection point & normal are in world space!

	Point3D x1 = worldToModel * ray.origin; // o
	Point3D x2 = worldToModel * (ray.origin + ray.dir);
	Point3D x0 = Point3D(0, 0, 0); // c

	Vector3D v0 = x0 - x1, v1 = x0 - x2, v2 = x2 - x1;
	Vector3D v = v0.cross(v1);
	float dist = v.length() / v2.length();

	if (!(dist > 1))
	{
		Vector3D I = x2 - x1;
		I.normalize();
		Vector3D J = x1 - x0;
		float d0 = -I.dot(J) - sqrt(I.dot(J) * I.dot(J) - J.length() * J.length() + 1);
		float d1 = -I.dot(J) + sqrt(I.dot(J) * I.dot(J) - J.length() * J.length() + 1);

		Point3D i_model; // intersection point in object space
		Point3D new_intersection; // intersection point in world space
		if (d0 > 0)
		{
			i_model = x1 + d0 * I;
			new_intersection = modelToWorld * i_model;
		}
		else if (d1 > 0)
		{
			i_model = x1 + d1 * I;
			new_intersection = modelToWorld * i_model;
		}

		float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
		if (new_t_value < ray.intersection.t_value)
		{
			ray.intersection.none = false;
			ray.intersection.point = new_intersection;
			ray.intersection.t_value = new_t_value;
			ray.intersection.normal = transNorm(worldToModel, i_model - x0);
			return true;
		}
		else
			return false;
	}
	else
		return false;

	//return !ray.intersection.none;
}

