#pragma once

#ifndef PBRT_SAMPLERS_STRATIFIED_H
#define PBRT_SAMPLERS_STRATIFIED_H

#include "sampler.h"
#include "rng.h"

class StratifiedSampler : public PixelSampler {
public:
	// StratifiedSampler Public Methods
	StratifiedSampler(int xPixelSamples, int yPixelSamples, bool jitterSamples,
		int nSampledDimensions)
		: PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions),
		xPixelSamples(xPixelSamples),
		yPixelSamples(yPixelSamples),
		jitterSamples(jitterSamples) {}
	void StartPixel(const Point2i &);
	std::unique_ptr<Sampler> Clone(int seed);

private:
	// StratifiedSampler Private Data
	const int xPixelSamples, yPixelSamples;
	const bool jitterSamples;
};

StratifiedSampler *CreateStratifiedSampler();










#endif