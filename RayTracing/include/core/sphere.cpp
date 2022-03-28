#include "sphere.h"


Bounds3f Sphere::ObjectBound() const {
	return Bounds3f(Point3f(-radius, -radius, -radius),
		Point3f(radius, radius, radius));
}


bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
	bool testAlphaTexture) const {
	
	Point3f pHit;
	Ray ray = (*WorldToObject)(r);
	Vector3f oc = ray.o - Point3f(0.0f, 0.0f, 0.0f);
	// Compute quadratic sphere coefficients

	// Initialize _EFloat_ ray coordinate values
	Float a = Dot(ray.d,ray.d);
	Float b = 2.0*Dot(oc, ray.d);
	Float c = Dot(oc, oc) - radius * radius;
	Float discriminant = b * b - 4 * a*c;
	return discriminant > 0;
}

bool Sphere::IntersectP(const Ray &r, bool testAlphaTexture) const {

	Point3f pHit;
	Ray ray = (*WorldToObject)(r);
	Vector3f oc = ray.o - Point3f(0.0f, 0.0f, 0.0f);
	// Compute quadratic sphere coefficients

	// Initialize _EFloat_ ray coordinate values
	Float a = Dot(ray.d, ray.d);
	Float b = 2.0*Dot(oc, ray.d);
	Float c = Dot(oc, oc) - radius * radius;
	Float discriminant = b * b - 4 * a*c;
	return discriminant > 0;
}