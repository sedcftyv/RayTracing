#pragma once

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

#include "pbrt.h"
#include "geometry.h"
#include "transform.h"

struct CameraSample {
	Point2f pFilm;
	Point2f pLens;
	Float time;
};

class Camera {
public:
	// Camera Interface
	Camera(const Transform &CameraToWorld);
		//Float shutterOpen,Float shutterClose, Film *film, const Medium *medium);
	virtual ~Camera();
	virtual Float GenerateRay(const CameraSample &sample, Ray *ray) const = 0;
	//virtual Float GenerateRayDifferential(const CameraSample &sample,RayDifferential *rd) const;
	//virtual Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
	//virtual void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
	//virtual Spectrum Sample_Wi(const Interaction &ref, const Point2f &u,Vector3f *wi, Float *pdf, Point2f *pRaster,VisibilityTester *vis) const;

	// Camera Public Data
	Transform CameraToWorld;
	//const Float shutterOpen, shutterClose;
	//Film *film;
	//const Medium *medium;
};

class ProjectiveCamera : public Camera {
public:
	// ProjectiveCamera Public Methods
	ProjectiveCamera(const Transform &CameraToWorld,
		const Transform &CameraToScreen,
		const Bounds2f &screenWindow,
		Float lensr, Float focald,int image_width=300,int image_height=300)
		: Camera(CameraToWorld),
		CameraToScreen(CameraToScreen) {
		// Initialize depth of field parameters
		lensRadius = lensr;
		focalDistance = focald;

		// Compute projective camera transformations

		// Compute projective camera screen transformations
		ScreenToRaster = Scale(image_width, image_height, 1) *
			Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x),
				1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
			Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
		RasterToScreen = Inverse(ScreenToRaster);
		RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
	}

protected:
	// ProjectiveCamera Protected Data
	Transform CameraToScreen, RasterToCamera;
	Transform ScreenToRaster, RasterToScreen;
	Float lensRadius, focalDistance;
};

class OrthographicCamera : public ProjectiveCamera {
public:
	// OrthographicCamera Public Methods
	OrthographicCamera(const Transform &CameraToWorld,
		const Bounds2f &screenWindow, Float lensRadius,Float focalDistance)
		: ProjectiveCamera(CameraToWorld, Orthographic(0, 1), screenWindow,lensRadius, focalDistance) {
		// Compute differential changes in origin for orthographic camera rays
		//dxCamera = RasterToCamera(Vector3f(1, 0, 0));
		//dyCamera = RasterToCamera(Vector3f(0, 1, 0));
	}
	Float GenerateRay(const CameraSample &sample, Ray *) const;

private:
	// OrthographicCamera Private Data
	//Vector3f dxCamera, dyCamera;
};

OrthographicCamera *CreateOrthographicCamera(const Transform &cam2world);

class PerspectiveCamera : public ProjectiveCamera {
public:
	// PerspectiveCamera Public Methods
	PerspectiveCamera(const Transform &CameraToWorld, const Bounds2f &screenWindow, Float lensRadius, Float focalDistance, Float fov,int image_width, int image_height)
		:ProjectiveCamera(CameraToWorld, Perspective(fov, focalDistance, 1000.f),
			screenWindow,lensRadius,focalDistance, image_width, image_height)
	{
		// Compute image plane bounds at z=1 for _PerspectiveCamera_
		//Point2i res = Point2i(300, 300);
		//Point3f pMin = RasterToCamera(Point3f(0, 0, 0));
		//Point3f pMax = RasterToCamera(Point3f(res.x, res.y, 0));
		//pMin /= pMin.z;
		//pMax /= pMax.z;
		//A = std::abs((pMax.x - pMin.x) * (pMax.y - pMin.y));
	}
	Float GenerateRay(const CameraSample &sample, Ray *) const;
private:
	// PerspectiveCamera Private Data
	//Vector3f dxCamera, dyCamera;
	//Float A;
};

PerspectiveCamera *CreatePerspectiveCamera(const Transform &cam2world, Float focaldistance, Float fov, int image_width, int image_height);





#endif