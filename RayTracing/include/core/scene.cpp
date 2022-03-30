#include "scene.h"

static long long nIntersectionTests = 0;
static long long nShadowTests = 0;

bool Scene::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
	++nIntersectionTests;
	//DCHECK_NE(ray.d, Vector3f(0, 0, 0));
	return aggregate->Intersect(ray, isect);
}

bool Scene::IntersectP(const Ray &ray) const {
	++nShadowTests;
	//DCHECK_NE(ray.d, Vector3f(0, 0, 0));
	return aggregate->IntersectP(ray);
}