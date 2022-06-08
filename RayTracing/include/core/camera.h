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
		Camera(const Transform &CameraToWorld);
		virtual ~Camera();
		virtual Float GenerateRay(const CameraSample &sample, Ray *ray) const = 0;				
		Transform CameraToWorld;
};

class ProjectiveCamera : public Camera {
public:
		ProjectiveCamera(const Transform &CameraToWorld,
		const Transform &CameraToScreen,
		const Bounds2f &screenWindow,
		Float lensr, Float focald,int image_width=300,int image_height=300)
		: Camera(CameraToWorld),
		CameraToScreen(CameraToScreen) {
				lensRadius = lensr;
		focalDistance = focald;
		ScreenToRaster = Scale(image_width, image_height, 1) *
		Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x),1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
		Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
		RasterToScreen = Inverse(ScreenToRaster);
		RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
		}

protected:
		Transform CameraToScreen, RasterToCamera;
		Transform ScreenToRaster, RasterToScreen;
		Float lensRadius, focalDistance;
};

class OrthographicCamera : public ProjectiveCamera {
public:
		OrthographicCamera(const Transform &CameraToWorld,
		const Bounds2f &screenWindow, Float lensRadius,Float focalDistance)
		: ProjectiveCamera(CameraToWorld, Orthographic(0, 1), screenWindow,lensRadius, focalDistance) {}
		Float GenerateRay(const CameraSample &sample, Ray *) const;
};

OrthographicCamera *CreateOrthographicCamera(const Transform &cam2world);

class PerspectiveCamera : public ProjectiveCamera {
public:
		PerspectiveCamera(const Transform &CameraToWorld, const Bounds2f &screenWindow, Float lensRadius, Float focalDistance, Float fov,int image_width, int image_height)
		:ProjectiveCamera(CameraToWorld, Perspective(fov, focalDistance, 1000.f),
			screenWindow,lensRadius,focalDistance, image_width, image_height){}
		Float GenerateRay(const CameraSample &sample, Ray *) const;
};

PerspectiveCamera *CreatePerspectiveCamera(const Transform &cam2world, Float focaldistance, Float fov, int image_width, int image_height);





#endif