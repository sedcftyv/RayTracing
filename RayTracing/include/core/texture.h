#pragma once

#ifndef PBRT_CORE_TEXTURE_H
#define PBRT_CORE_TEXTURE_H

#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"

template <typename T>
class Texture {
public:
	// Texture Interface
	virtual T Evaluate(const SurfaceInteraction &) const = 0;
	virtual ~Texture() {}
};

template <typename T>
class ConstantTexture : public Texture<T> {
public:
	// ConstantTexture Public Methods
	ConstantTexture(const T &value) : value(value) {}
	T Evaluate(const SurfaceInteraction &) const { return value; }
private:
	T value;
};





#endif