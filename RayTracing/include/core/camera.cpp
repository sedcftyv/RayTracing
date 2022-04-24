#include "camera.h"
#include "assert.h"
Camera::~Camera() { }

Camera::Camera(const Transform &CameraToWorld): CameraToWorld(CameraToWorld){
	if (CameraToWorld.HasScale())
		assert(
			"Scaling detected in world-to-camera transformation!\n"
			"The system has numerous assumptions, implicit and explicit,\n"
			"that this transform will have no scale factors in it.\n"
			"Proceed at your own risk; your image may have errors or\n"
			"the system may crash as a result of this.");

}


Float OrthographicCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	//ProfilePhase prof(Prof::GenerateCameraRay);
	// Compute raster and camera sample positions
	Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
	Point3f pCamera = RasterToCamera(pFilm);
	*ray = Ray(pCamera, Vector3f(0, 0, 1));
	// Modify ray for depth of field
	*ray = CameraToWorld(*ray);
	return 1;
}

OrthographicCamera *CreateOrthographicCamera(const Transform &cam2world) {
	Float lensradius = 0.0f;
	Float focaldistance = 0.0f;
	Float frame = 300.f / 300.f;
	Bounds2f screen;
	if (frame > 1.f) {
		screen.pMin.x = -frame;
		screen.pMax.x = frame;
		screen.pMin.y = -1.f;
		screen.pMax.y = 1.f;
	}
	else {
		screen.pMin.x = -1.f;
		screen.pMax.x = 1.f;
		screen.pMin.y = -1.f / frame;
		screen.pMax.y = 1.f / frame;
	}

	float ScreenScale = 2.0f;
	screen.pMin.x *= ScreenScale;
	screen.pMax.x *= ScreenScale;
	screen.pMin.y *= ScreenScale;
	screen.pMax.y *= ScreenScale;

	return new OrthographicCamera(cam2world, screen,lensradius, focaldistance);
}

Float PerspectiveCamera::GenerateRay(const CameraSample &sample,Ray *ray) const {
	//ProfilePhase prof(Prof::GenerateCameraRay);
	// Compute raster and camera sample positions
	Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
	Point3f pCamera = RasterToCamera(pFilm);
	*ray = Ray(Point3f(0, 0, 0), Normalize(Vector3f(pCamera)));
	*ray = CameraToWorld(*ray);
	return 1;
}

PerspectiveCamera *CreatePerspectiveCamera(const Transform &cam2world,int image_width,int image_height)
{
	Float lensradius = 0.0f;
	Float focaldistance = 0.0f;
	Float frame = (float)image_width / (float)image_height;
	Bounds2f screen;
	if (frame > 1.f) {
		screen.pMin.x = -frame;
		screen.pMax.x = frame;
		screen.pMin.y = -1.f;
		screen.pMax.y = 1.f;
	}
	else {
		screen.pMin.x = -1.f;
		screen.pMax.x = 1.f;
		screen.pMin.y = -1.f / frame;
		screen.pMax.y = 1.f / frame;
	}
	Float fov = 90.0f;
	Float halffov = 45.0f;
	return new PerspectiveCamera(cam2world, screen, lensradius, focaldistance, fov, image_width, image_height);
}
