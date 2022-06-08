#pragma once

#ifndef PBRT_H
#define PBRT_H

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <assert.h>
#include <malloc.h>
#include <map>
#include <cctype>

using std::vector;
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;
using std::make_unique;
using std::string;


//#define CHECK_NE(a)
//#define CHECK_EQ(a)
//#define CHECK_GE(a)
//#define CHECK_LT(a)
//#define DCHECK(a)
//#define DCHECK_NE(a)
//#define CHECK(a)
//#define CHECK_GT(a)

#ifndef PBRT_L1_CACHE_LINE_SIZE
#define PBRT_L1_CACHE_LINE_SIZE 64
#endif



#define PBRT_CONSTEXPR constexpr
#ifdef PBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 
template <typename T>
class Vector2;
template <typename T>
class Vector3;
template <typename T>
class Point3;
template <typename T>
class Point2;
template <typename T>
class Normal3;
class Ray;
class RayDifferential;
template <typename T>
class Bounds2;
template <typename T>
class Bounds3;
class Transform;
struct Interaction;
class SurfaceInteraction;
class Shape;
class Primitive;
class GeometricPrimitive;
class TransformedPrimitive;
class MediumInterface;
class BxDF;
class BSDF;
class RGBSpectrum;
template <typename T>
class Texture;
class Material;
typedef RGBSpectrum Spectrum;
class Scene;
class Light;
class VisibilityTester;
class AreaLight;
class Sampler;
class TriangleMesh;
class Distribution2D;
class MicrofacetDistribution;
//inline void *AllocAligned(size_t size) {
//	return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
//}


template <typename T, int logBlockSize=2>
class BlockedArray {
public:
	BlockedArray(int uRes, int vRes, const T *d = nullptr)
	: uRes(uRes), vRes(vRes), uBlocks(RoundUp(uRes) >> logBlockSize) {
	int nAlloc = RoundUp(uRes) * RoundUp(vRes);
	data = new T[nAlloc];
	for (int i = 0; i < nAlloc; ++i) new (&data[i]) T();
	if (d)
		for (int v = 0; v < vRes; ++v)
			for (int u = 0; u < uRes; ++u) 
				(*this)(u, v) = d[v * uRes + u];
	}
	PBRT_CONSTEXPR int BlockSize() const { return 1 << logBlockSize; }
	int RoundUp(int x) const {
		return (x + BlockSize() - 1) & ~(BlockSize() - 1);
	}
	int uSize() const { return uRes; }
	int vSize() const { return vRes; }
	~BlockedArray() {
		for (int i = 0; i < uRes * vRes; ++i) data[i].~T();
		delete[]data;
	}
	int Block(int a) const { return a >> logBlockSize; }
	int Offset(int a) const { return (a & (BlockSize() - 1)); }
	T &operator()(int u, int v) {
		int bu = Block(u), bv = Block(v);
		int ou = Offset(u), ov = Offset(v);
		int offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
		offset += BlockSize() * ov + ou;
		return data[offset];
	}
	const T &operator()(int u, int v) const {
		int bu = Block(u), bv = Block(v);
		int ou = Offset(u), ov = Offset(v);
		int offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
		offset += BlockSize() * ov + ou;
		return data[offset];
	}
	void GetLinearArray(T *a) const {
		for (int v = 0; v < vRes; ++v)
			for (int u = 0; u < uRes; ++u) 
				*a++ = (*this)(u, v);
	}

private:
	T *data;
	const int uRes, vRes, uBlocks;
};

inline bool HasExtension(const std::string &value, const std::string &ending) {
	if (ending.size() > value.size()) return false;
	return std::equal(ending.rbegin(), ending.rend(), value.rbegin(),
		[](char a, char b) { return std::tolower(a) == std::tolower(b); });
}

template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}

template <typename T> inline T Mod(T a, T b) {
	T result = a - (a / b) * b;
	return (T)((result < 0) ? result + b : result);
}

template <> inline Float Mod(Float a, Float b) {
	return std::fmod(a, b);
}

static const Float Pi = 3.14159265358979323846;
static const Float InvPi = 0.31830988618379067154;
static const Float Inv2Pi = 0.15915494309189533577;
static const Float Inv4Pi = 0.07957747154594766788;
static const Float PiOver2 = 1.57079632679489661923;
static const Float PiOver4 = 0.78539816339744830961;
static const Float Sqrt2 = 1.41421356237309504880;
const Float ShadowEpsilon = 0.0001f;
const Float eps = 1e-4;
static constexpr Float MachineEpsilon = std::numeric_limits<Float>::epsilon() * 0.5;
static constexpr Float Infinity = std::numeric_limits<Float>::infinity();

inline Float Radians(Float deg) {
	return (Pi / 180) * deg;
}
inline Float Degrees(Float rad) {
	return (180 / Pi) * rad;
}

inline Float Log2(Float x) {
	const Float invLog2 = 1.442695040888963387004650940071;
	return std::log(x) * invLog2;
}

inline int Log2Int(uint32_t v) {
	unsigned long lz = 0;
	if (_BitScanReverse(&lz, v)) return lz;
	return 0;
}
inline int CountTrailingZeros(uint32_t v) {
	unsigned long index;
	if (_BitScanForward(&index, v))
		return index;
	else
		return 32;
}

template <typename T> inline bool IsPowerOf2(T v) {
	return v && !(v & (v - 1));
}

inline int32_t RoundUpPow2(int32_t v) {
	v--;
	v |= v >> 1;    v |= v >> 2;
	v |= v >> 4;    v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

template <typename Predicate> int FindInterval(int size,
	const Predicate &pred) {
	int first = 0, len = size;
	while (len > 0) {
		int half = len >> 1, middle = first + half;
		if (pred(middle)) {
			first = middle + 1;
			len -= half + 1;
		}
		else
			len = half;
	}
	return Clamp(first - 1, 0, size - 2);
}

inline constexpr Float gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline Float GammaCorrect(Float value) {
	if (value <= 0.0031308f) return 12.92f * value;
	return 1.055f * std::pow(value, (Float)(1.f / 2.4f)) - 0.055f;
}

inline Float InverseGammaCorrect(Float value) {
	if (value <= 0.04045f) return value * 1.f / 12.92f;
	return std::pow((value + 0.055f) * 1.f / 1.055f, (Float)2.4f);
}

inline Float Lerp(Float t, Float v1, Float v2) 
{ return (1 - t) * v1 + t * v2; }

inline Float ErfInv(Float x) {
	Float w, p;
	x = Clamp(x, -.99999f, .99999f);
	w = -std::log((1 - x) * (1 + x));
	if (w < 5) {
		w = w - 2.5f;
		p = 2.81022636e-08f;
		p = 3.43273939e-07f + p * w;
		p = -3.5233877e-06f + p * w;
		p = -4.39150654e-06f + p * w;
		p = 0.00021858087f + p * w;
		p = -0.00125372503f + p * w;
		p = -0.00417768164f + p * w;
		p = 0.246640727f + p * w;
		p = 1.50140941f + p * w;
	}
	else {
		w = std::sqrt(w) - 3;
		p = -0.000200214257f;
		p = 0.000100950558f + p * w;
		p = 0.00134934322f + p * w;
		p = -0.00367342844f + p * w;
		p = 0.00573950773f + p * w;
		p = -0.0076224613f + p * w;
		p = 0.00943887047f + p * w;
		p = 1.00167406f + p * w;
		p = 2.83297682f + p * w;
	}
	return p * x;
}

inline Float Erf(Float x) {
	Float a1 = 0.254829592f;
	Float a2 = -0.284496736f;
	Float a3 = 1.421413741f;
	Float a4 = -1.453152027f;
	Float a5 = 1.061405429f;
	Float p = 0.3275911f;

	int sign = 1;
	if (x < 0) sign = -1;
	x = std::abs(x);

	Float t = 1 / (1 + p * x);
	Float y =1 -(((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

	return sign * y;
}



#endif