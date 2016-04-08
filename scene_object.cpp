/***********************************************************
     Starter code for Assignment 3
     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005
		implements scene_object.h
***********************************************************/

#include <cmath>
#include <iostream>
#include <stdio.h>
#include "scene_object.h"

const double eps = 1e-4;

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
	if (std::abs(dirInObject[2]) < eps) return false;
	float dist = -originInObject[2] / dirInObject[2];
	if (dist < eps) return false; // the ray shoots away or the origin is on the plane
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
		if (ray.intersection.none || new_t_value + eps < ray.intersection.t_value)
		{
			ray.intersection.none = false;
			ray.intersection.t_value = new_t_value;
			ray.intersection.point = new_intersection;
			ray.intersection.normal = transNorm(worldToModel, normal);
			ray.intersection.normal.normalize();
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
		if (d0 > eps)
		{
			i_model = x1 + d0 * I;
			new_intersection = modelToWorld * i_model;
		}
		else if (d1 > eps)
		{
			i_model = x1 + d1 * I;
			new_intersection = modelToWorld * i_model;
		}
		else
			return false;

		float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
		if (ray.intersection.none || new_t_value + eps < ray.intersection.t_value)
		{
			ray.intersection.pointOS = i_model;
			ray.intersection.none = false;
			ray.intersection.point = new_intersection;
			ray.intersection.t_value = new_t_value;
			ray.intersection.normal = transNorm(worldToModel, i_model - x0);
			ray.intersection.normal.normalize();
			return true;
		}
		else
			return false;
	}
	else
		return false;
}


bool UnitCone::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ){

	Point3D O(0, 0, 0);
        Point3D origin = worldToModel * ray.origin;
        Vector3D dir = worldToModel * ray.dir;

        float twoPI = (float)(M_PI * 2.0);
        float A=1;
        float B=1;
        float C=1;
        float A2, B2, C2;
        A2 = 2.0f * A;
        B2 = 2.0f * B;
        C2 = 2.0f * C;
        float a2_1;
        float a, b, c, discrim;
        float dMagOS;
        float zmin = 0;
        float zmax = 1;
        float miny, maxy;
        float tWS;
        float theta;
        float thetamax = 45;
        float t_valueWorld;

        miny = -zmin - eps; 
        maxy = zmax + eps;
        theta = (float)(thetamax * M_PI / 180.0f - M_PI);

        dMagOS = 1.0f / dir.length();  
        bool flag = false;

        a = A * dir[0] * dir[0] + B * dir[1] * dir[1] - C * dir[2] * dir[2];
        b = A2 * origin[0] * dir[0] + B2 * origin[1] * dir[1] - C2 * origin[2] * dir[2];
        c = A * origin[0] * origin[0] + B * origin[1] * origin[1] - C * origin[2] * origin[2];

	// check if intersect with the bottum round
	if (std::abs(dir[2]) > eps)
	{
		float t_value = (1 - origin[2]) / dir[2];
		Point3D I = origin + t_value * dir;

		// check if the ray shoots away or the origin is on the plane
		// check if the ray is in range
		if (t_value > eps && I[0] * I[0] + I[1] * I[1] < 1)
		{
			Point3D new_intersection = modelToWorld * I;
			float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
			if (ray.intersection.none || new_t_value + eps < ray.intersection.t_value)
			{
				ray.intersection.none = false;
				ray.intersection.point = new_intersection;
				ray.intersection.t_value = new_t_value;
				ray.intersection.normal = transNorm(worldToModel, Vector3D(0, 0, 1));
				ray.intersection.normal.normalize();
				flag = true;
			}
		}
	}

        if (std::abs(a) < eps)
	{
		float t_value = -c / b;
		Point3D intersection_model = origin + t_value * dir;
		if (t_value > eps && 0 <= intersection_model[2] && intersection_model[2] <= 1)
		{
			Point3D new_intersection = modelToWorld * intersection_model;
			float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
			if (ray.intersection.none || new_t_value + eps < ray.intersection.t_value)
			{
				ray.intersection.none = false;
				ray.intersection.point = new_intersection;
				ray.intersection.t_value = new_t_value;
				Vector3D normal_model = intersection_model - O;
				normal_model[2] = -normal_model[2];
				ray.intersection.normal = transNorm(worldToModel, normal_model);
				ray.intersection.normal.normalize();
				flag = true;
			}
		}
	}
	else
	{
		discrim = b * b - 4.0 * a * c;

		if (discrim < 0) return flag; // no intersection

		discrim = (float)sqrt(discrim);

		float t_value = (-b - discrim) / (2 * a);            // near intersection
		Point3D intersection_model = origin + t_value * dir;

		if (t_value > eps && 0 <= intersection_model[2] && intersection_model[2] <= 1)
		{
			Point3D new_intersection = modelToWorld * intersection_model;
			float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
			if (ray.intersection.none || new_t_value + eps < ray.intersection.t_value)
			{
				ray.intersection.none = false;
				ray.intersection.point = new_intersection;
				ray.intersection.t_value = new_t_value;
				Vector3D normal_model = intersection_model - O;
				normal_model[2] = -normal_model[2];
				ray.intersection.normal = transNorm(worldToModel, normal_model);
				ray.intersection.normal.normalize();
				flag = true;
			}
		}

		t_value = (-b + discrim) / (2 * a);            // far intersection
		intersection_model = origin + t_value * dir;
		if (t_value > eps && 0 <= intersection_model[2] && intersection_model[2] <= 1)
		{
			Point3D new_intersection = modelToWorld * intersection_model;
			float new_t_value = (new_intersection[0] - ray.origin[0]) / ray.dir[0];
			if (ray.intersection.none || new_t_value + eps < ray.intersection.t_value)
			{
				ray.intersection.none = false;
				ray.intersection.point = new_intersection;
				ray.intersection.t_value = new_t_value;
				Vector3D normal_model = intersection_model - O;
				normal_model[2] = -normal_model[2];
				ray.intersection.normal = transNorm(worldToModel, normal_model);
				ray.intersection.normal.normalize();
				flag = true;
			}
		}
        }

	return flag;
}


