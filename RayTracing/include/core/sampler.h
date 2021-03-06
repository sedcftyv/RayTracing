#pragma once

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include "camera.h"



class Sampler {
public:
	virtual ~Sampler();
	Sampler(int64_t samplesPerPixel);
	virtual void StartPixel(const Point2i &p);
	virtual Float Get1D() = 0;
	virtual Point2f Get2D() = 0;
	CameraSample GetCameraSample(const Point2i &pRaster);
	void Request1DArray(int n);
	void Request2DArray(int n);
	virtual int RoundCount(int n) const { return n; }
	const Float *Get1DArray(int n);
	const Point2f *Get2DArray(int n);
	virtual bool StartNextSample();
	virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
	virtual bool SetSampleNumber(int64_t sampleNum);
	int64_t CurrentSampleNumber() const { return currentPixelSampleIndex; }
	const int64_t samplesPerPixel;

protected:
	Point2i currentPixel;
	int64_t currentPixelSampleIndex;
	std::vector<int> samples1DArraySizes, samples2DArraySizes;
	std::vector<std::vector<Float>> sampleArray1D;
	std::vector<std::vector<Point2f>> sampleArray2D;

private:
	size_t array1DOffset, array2DOffset;
};

class PixelSampler : public Sampler {
public:
	PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);
	bool StartNextSample();
	bool SetSampleNumber(int64_t);
	Float Get1D();
	Point2f Get2D();

protected:
	std::vector<std::vector<Float>> samples1D;
	std::vector<std::vector<Point2f>> samples2D;
	int current1DDimension = 0, current2DDimension = 0;
	RNG rng;
};





#endif