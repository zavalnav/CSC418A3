/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>

#define NUMANTIALIASE 5
#define ANTIALIAS 1
#define REFLECTION 0
#define GLOSSY_REFLECTION 0 // only available when REFLECTION = 1
#define REFRACTION 0
#define SOFT_SHADOWS 0
#define FOCALLENGTH 5
#define DOF 1

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
	_maxDepth = 2;
	_sample_num = 100;
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			//printf("intersected an object\n");
 			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows here if needed.

		if (SOFT_SHADOWS)
		{
			Colour final_colour(0, 0, 0);
			curLight->light->ambient(ray);
			Colour ambient_colour = ray.col;
			curLight->light->shade(ray);
			Colour phong_colour = ray.col;

			srand((unsigned)time(NULL));

			for (int i = 0; i < _sample_num; ++ i)
			{
				Vector3D dir_shadow = curLight->light->get_random_sample() - ray.intersection.point;
				Ray3D ray_shadow(ray.intersection.point, dir_shadow);
				traverseScene(_root, ray_shadow);
				if (!ray_shadow.intersection.none)
					final_colour = final_colour + ambient_colour;
				else
					final_colour = final_colour + phong_colour;
			}
			ray.col = 1.0 / _sample_num * final_colour;
		}
		else
		{
			Vector3D dir_shadow = curLight->light->get_position() - ray.intersection.point;
			Ray3D ray_shadow(ray.intersection.point, dir_shadow);
			traverseScene(_root, ray_shadow);
			if (!ray_shadow.intersection.none)
			{
				curLight->light->ambient(ray);
				printf("mm\n");
			}
			else
				curLight->light->shade(ray);
		}
		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray, int depth, bool air ) {
	Colour col(0.0, 0.0, 0.0); 
	if (depth == _maxDepth) return col;
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray); 
		float refl = ray.intersection.mat->refl_coef;
		float refr = ray.intersection.mat->refr_coef;
//		col = (1 - refl - refr) * ray.col;
		col = ray.col;

		if (REFLECTION && air && refl > 0)
		{
			if (GLOSSY_REFLECTION && ray.intersection.mat->refl_sigma > 0)
			{
				// handle blurred reflection
				float nlength = -2 * ray.dir.dot(ray.intersection.normal);
				Vector3D refl_dir = nlength * ray.intersection.normal + ray.dir;
				refl_dir.normalize();

				Colour refl_col(0.0, 0.0, 0.0);

				for (int i = 0; i < _sample_num; ++ i)
				{
					Vector3D u(-refl_dir[1], refl_dir[0], 0);
					u.normalize();
					Vector3D v = u.cross(refl_dir);
					double a = (double)rand() / RAND_MAX * 2 * M_PI;
					double r = (double)rand() / RAND_MAX * ray.intersection.mat->refl_sigma;
					Vector3D random_dir = refl_dir + r * (sin(a) * u + cos(a) * v);

					Ray3D ray_refl(ray.intersection.point, random_dir);
					Colour this_col = shadeRay( ray_refl, depth + 1, air );
					refl_col = refl_col + refl * this_col;
				}
				col = col + refl / _sample_num * refl_col;
			}
			else
			{
				// handle reflection
				float nlength = -2 * ray.dir.dot(ray.intersection.normal);
				Ray3D ray_refl(ray.intersection.point, nlength * ray.intersection.normal + ray.dir);

				Colour refl_col = shadeRay( ray_refl, depth + 1, air );
				col = col + refl * refl_col;
			}
		}

		if (REFRACTION && refr > 0)
		{
			// handle refraction
			Ray3D ray_refr;
			ray_refr.origin = ray.intersection.point;

			if (air) // shoot ray from air to material
			{
				double nlength = -ray.dir.dot(ray.intersection.normal);
				Vector3D d1 = ray.dir + nlength * ray.intersection.normal;
				double cos1 = nlength / ray.dir.length();
				double sin1 = sqrt(1 - cos1 * cos1);
				double n2 = ray.intersection.mat->refr_index;
				double sin2 = sin1 / n2;
				double cos2 = sqrt(1 - sin2 * sin2);
				double rlength = nlength / cos2 * sin2;
				Vector3D d2 = rlength / d1.length() * d1;
				ray_refr.dir = ray.dir - d1 + d2;
				//printf("%lf\n", ray_refr.dir.cross(ray.intersection.normal).dot(ray.dir));
				Colour refr_col = shadeRay( ray_refr, depth + 1, !air );
				col = col + refr * refr_col;
			}
			else // shoot ray from material to air
			{
				double nlength = ray.dir.dot(ray.intersection.normal);
				Vector3D d1 = ray.dir - nlength * ray.intersection.normal;
				double cos1 = nlength / ray.dir.length();
				double sin1 = sqrt(1 - cos1 * cos1);
				double n1 = ray.intersection.mat->refr_index;
				double sin2 = sin1 * n1;
				if (sin2 < 1)
				{
					double cos2 = sqrt(1 - sin2 * sin2);
					double rlength = nlength / cos2 * sin2;
					Vector3D d2 = rlength / d1.length() * d1;
					ray_refr.dir = ray.dir - d1 + d2;
					Colour refr_col = shadeRay( ray_refr, depth + 1, !air );
					//printf("%lf\n", ray_refr.dir.cross(ray.intersection.normal).dot(ray.dir));
					//printf("%lf %lf %lf\n", refr_col[0], refr_col[1], refr_col[2]);
					col = col + refr_col; // assume no reflection inside material
				}
			}
		}

		col.clamp();
	}
	
	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	freopen("scene.txt", "w", stdout);

	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) 
	{
		for (int j = 0; j < _scrWidth; j++) 
		{
			// Sets up ray origin and direction in view space, image plane is at z = -1.
			Point3D origin(0, 0, 0); // Create a point called origin 
			Point3D imagePlane;	// the image plane where it is some distance away from the camera 

			imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor; // x position of image plane (factor influences the final width of the image plane) 
			imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor; // y position of image plane 
			Colour col(0,0,0); 
			srand(time(NULL));  

			for (float parti = i; parti < i + 1.0f; parti += 0.5f){
				for (float partj = j; partj < j + 1.0f; partj += 0.5f){
					imagePlane[0] = (-double(width)/2 + ((double)rand()/(double)RAND_MAX) + 0.25 + parti)/factor;
					imagePlane[1] = (-double(height)/2 + ((double)rand()/(double)RAND_MAX) + 0.25 + partj)/factor;
					imagePlane[2] = -1; 
					// Construct the ray direction based on the points on the image plane computed 
					Vector3D dir = Vector3D(imagePlane[0], imagePlane[1], imagePlane[2]);
					// point Aimed is the position of the focal plane in specified direction 
					Point3D pointAimed; 
					Point3D new_origin = eye; 
					
					if(DOF == 1) 
					{
						dir = viewToWorld * dir;
						dir.normalize();
						pointAimed = eye + FOCALLENGTH*dir; 
						double radius = 1; // radius of the lens 
		    			double rand_x = rand()/double(RAND_MAX+1);
		    			double rand_y = rand()/double(RAND_MAX+1);
						Vector3D u = up; 
						Vector3D v = up.cross(view); 
						new_origin[0]= eye[0] -(radius/2)*u[0]-(radius/2)*v[0]+radius*(rand_x)*u[0]+radius*(rand_y)*v[0]; 
						new_origin[1]= eye[1] -(radius/2)*u[1]-(radius/2)*v[1]+radius*(rand_x)*u[1]+radius*(rand_y)*v[1]; 
						new_origin[2]= eye[2] -(radius/2)*u[2]-(radius/2)*v[2]+radius*(rand_x)*u[2]+radius*(rand_y)*v[2]; 
						//getting the new direction of ray
						dir = pointAimed - new_origin;
						Ray3D ray = Ray3D(eye, dir);
						dir.normalize();
					}
					else
					{
						dir = viewToWorld * dir;
						dir.normalize();
					}
					Ray3D ray = Ray3D(new_origin, dir); 	
					Colour colTwo = shadeRay(ray, 0, true); 
					col[0] += colTwo[0]; 
					col[1] += colTwo[1]; 
					col[2] += colTwo[2]; 
					}
				}
				col[0] /= (double) NUMANTIALIASE;  
				col[1] /= (double) NUMANTIALIASE; 
				col[2] /= (double) NUMANTIALIASE;  
				// Finally, color the pixel based on the color that was given 
				// Add less and less each sweep 
				_rbuffer[i*width+j] += int(col[0]*255);
				_gbuffer[i*width+j] += int(col[1]*255);
				_bbuffer[i*width+j] += int(col[2]*255);
		}
	}
	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	raytracer.loadScene(width, height, 0);

	return 0;
}


void Raytracer::loadScene(int width, int height, int scene)
{
	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
			Colour(0.628281, 0.555802, 0.366065), 51.2,
			0.0, 0.0, 0.0, 0.0);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
			Colour(0.316228, 0.316228, 0.316228), 12.8,
			0.5, 1.0, 0.0, 0.0);
	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.588235, 0.670588, 0.729412),
			Colour(0.9, 0.9, 0.9), 1.5,
			1.0, 1.0, 1.0, 1.5);

	if (scene == 0)
	{
		// Camera parameters.
		Point3D eye(0, 0, 1);
		Vector3D view(0, 0, -1);
		Vector3D up(0, 1, 0);
		double fov = 60;

		// Defines a point light source.
		addLightSource( new PointLight(Point3D(0, 0, 5), 1, Colour(0.9, 0.9, 0.9)));
		// Add a unit sphere nto the scene with material mat.
		SceneDagNode* sphere = addObject( new UnitSphere(), &gold );
		SceneDagNode* sphere1 = addObject( new UnitSphere(), &gold );
		SceneDagNode* sphere2 = addObject( new UnitSphere(), &glass );

//		SceneDagNode* cone = addObject( new UnitCone(), &glass );

		//double f[3] = {0.4, 0.4, 0.4};
		//scale(cone, Point3D(0, 0, 0), f);

		//SceneDagNode* sphere3 = raytracer.addObject( new UnitSphere(), &gold );
		//SceneDagNode* sphere4 = raytracer.addObject( new UnitSphere(), &gold );

		SceneDagNode* plane = addObject( new UnitSquare(), &jade );
	
		// Apply some transformations to the unit sphere.
		//Sphere1
		double factor1[3] = { 1.0, 2.0, 1.0 };
		translate(sphere, Vector3D(0, 0, -5));
		rotate(sphere, 'x', -45); 
		rotate(sphere, 'z', 45); 
		scale(sphere, Point3D(0, 0, 0), factor1);

		//Sphere2
		double factor3[3] = { 0.5, 0.5, 0.5 };
		translate(sphere1, Vector3D(1, 1, -7));	
		rotate(sphere1, 'x', 30); 
		rotate(sphere1, 'z', 45); 
		scale(sphere1, Point3D(0, 0, 0), factor3);

		//Sphere3
		//raytracer.rotate(sphere2, 'x', 30); 
		//raytracer.rotate(sphere2, 'y', 45); 
		translate(sphere2, Vector3D(-3, 1, -5));	

		//Plane
		double factor2[3] = { 6.0, 6.0, 6.0 };
		translate(plane, Vector3D(0, 0, -7));	
		rotate(plane, 'z', 45); 
		scale(plane, Point3D(0, 0, 0), factor2);

		// Render the scene, feel free to make the image smaller for
		// testing purposes.	
		render(width, height, eye, view, up, fov, "view1.bmp");

		// Render it from a different point of view.
		Point3D eye2(4, 2, 1);
		Vector3D view2(-4, -2, -6);
		render(width, height, eye2, view2, up, fov, "view2.bmp");
	}
	if (scene == 1)
	{
		Point3D eye(0, 1, 4);
		Vector3D view(0, 0, -1);
		Vector3D up(0, 1, 0);
		double fov = 60;

		addLightSource( new PointLight(Point3D(0, -5, -1), 1, Colour(0.9, 0.9, 0.9)));

		for (int i = 0; i < 6; ++ i)
			for (int j = 0; j < 6; ++ j)
			{
				SceneDagNode* plane;
				if ((i + j) & 1)
					plane = addObject( new UnitSquare(), &jade );
				else
					plane = addObject( new UnitSquare(), &gold );
				rotate(plane, 'x', 90);
				translate(plane, Vector3D(i - 2.5, j - 2.5, 0.0));
			}

		SceneDagNode* sphere = addObject( new UnitSphere(), &glass );
		translate(sphere, Vector3D(0, 2, 0));

		render(width, height, eye, view, up, fov, "refraction.bmp");
	}
	if (scene == 2)
	{
		SceneDagNode* cone = addObject( new UnitCone(), &gold );
		rotate(cone, 'y', 90);
		addLightSource( new PointLight(Point3D(2, -1, 5), 1, Colour(0.9, 0.9, 0.9)));

		Point3D eye(2, 1, 4);
		Vector3D view(-2, 0, -4);
		Vector3D up(0, 1, 0);
		double fov = 60;
		render(width, height, eye, view, up, fov, "cone.bmp");		
	}
}