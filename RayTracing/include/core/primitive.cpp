#include "primitive.h"
#include "interaction.h"

static long long primitiveMemory = 0;

Primitive::~Primitive() {}

const AreaLight *Aggregate::GetAreaLight() const {
	return nullptr;
}

const Material *Aggregate::GetMaterial() const {
	return nullptr;
}

const AreaLight *GeometricPrimitive::GetAreaLight() const {
	return areaLight.get();
}

const Material *GeometricPrimitive::GetMaterial() const {
	return material.get();
}

GeometricPrimitive::GeometricPrimitive(const std::shared_ptr<Shape> &shape, const std::shared_ptr<Material> &material, const std::shared_ptr<AreaLight> &areaLight)
	: shape(shape),material(material), areaLight(areaLight){
	primitiveMemory += sizeof(*this);
}

Bounds3f GeometricPrimitive::WorldBound() const { return shape->WorldBound(); }

bool GeometricPrimitive::IntersectP(const Ray &r) const {
	return shape->IntersectP(r);
}

bool GeometricPrimitive::Intersect(const Ray &r,
	SurfaceInteraction *isect) const {
	Float tHit;
	if (!shape->Intersect(r, &tHit, isect)) return false;
	r.tMax = tHit;
	isect->primitive = this;
	return true;
}

void GeometricPrimitive::ComputeScatteringFunctions(SurfaceInteraction *isect, TransportMode mode,bool allowMultipleLobes) const {
	if (material)
		material->ComputeScatteringFunctions(isect, mode,allowMultipleLobes);
}

SimpleAccel::SimpleAccel(std::vector<std::shared_ptr<Primitive>> p) {
	if (p.empty()) return;
	primitives = p;
	for (int i = 0; i < primitives.size(); ++i)
	{
		bound=Union(primitives[i]->WorldBound(),bound);
	}
}

SimpleAccel::~SimpleAccel() {
}

Bounds3f SimpleAccel::WorldBound() const{
	return bound;
}
bool SimpleAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
	bool hit = false;
	for (int i = 0; i <primitives.size(); ++i)
		if (primitives[i]->Intersect(ray, isect))
			hit = true;
	return hit;
}

bool SimpleAccel::IntersectP(const Ray &ray) const {
	for (int i = 0; i < primitives.size(); ++i)
		if (primitives[i]->IntersectP(ray))
			return true;
	return false;
}