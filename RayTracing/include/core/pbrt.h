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

#define CHECK_NE(a)
#define CHECK_EQ(a)
#define CHECK_GE(a)
#define CHECK_LT(a)
#define DCHECK(a)
#define CHECK(a)

#define PBRT_CONSTEXPR constexpr
#ifdef PBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif // PBRT_FLOAT_AS_DOUBLE
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
class BSDF;
class BSSRDF;
class MemoryArena;
class RGBSpectrum;
class TransportMode;
typedef RGBSpectrum Spectrum;



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
		//<< Bisect range based on value of pred at middle >>
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

inline void FreeAligned(void *ptr) {
	if (!ptr) return;
#if defined(PBRT_HAVE__ALIGNED_MALLOC)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

#endif