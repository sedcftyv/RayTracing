#include "lightdistrib.h"
#include "scene.h"

LightDistribution::~LightDistribution() {}

std::unique_ptr<LightDistribution> CreateLightSampleDistribution(
	const std::string &name, const Scene &scene) {
		return std::unique_ptr<LightDistribution>{
		new UniformLightDistribution(scene)};
}

UniformLightDistribution::UniformLightDistribution(const Scene &scene) {
	std::vector<Float> prob(scene.lights.size(), Float(1));
	distrib.reset(new Distribution1D(&prob[0], int(prob.size())));
}

const Distribution1D *UniformLightDistribution::Lookup(const Point3f &p) const {
	return distrib.get();
}