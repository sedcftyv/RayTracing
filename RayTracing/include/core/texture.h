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

class TextureMapping2D {
public:
	// TextureMapping2D Interface
	virtual ~TextureMapping2D();
	virtual Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
		Vector2f *dstdy) const = 0;
};

class UVMapping2D : public TextureMapping2D {
public:
	// UVMapping2D Public Methods
	UVMapping2D(Float su = 1, Float sv = 1, Float du = 0, Float dv = 0);
	Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
		Vector2f *dstdy) const;

private:
	const Float su, sv, du, dv;
};

Float Lanczos(Float x, Float tau = 2);
#endif