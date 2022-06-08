#include "sampling.h"
#include "geometry.h"
#include "shape.h"

void StratifiedSample1D(Float *samp, int nSamples, RNG &rng, bool jitter) {
	Float invNSamples = (Float)1 / nSamples;
	for (int i = 0; i < nSamples; ++i) {
		Float delta = jitter ? rng.UniformFloat() : 0.5f;
		samp[i] = std::min((i + delta) * invNSamples, OneMinusEpsilon);
	}
}

void StratifiedSample2D(Point2f *samp, int nx, int ny, RNG &rng, bool jitter) {
	Float dx = (Float)1 / nx, dy = (Float)1 / ny;
	for (int y = 0; y < ny; ++y)
		for (int x = 0; x < nx; ++x) {
			Float jx = jitter ? rng.UniformFloat() : 0.5f;
			Float jy = jitter ? rng.UniformFloat() : 0.5f;
			samp->x = std::min((x + jx) * dx, OneMinusEpsilon);
			samp->y = std::min((y + jy) * dy, OneMinusEpsilon);
			++samp;
		}
}

void LatinHypercube(Float *samples, int nSamples, int nDim, RNG &rng) {
	Float invNSamples = (Float)1 / nSamples;
	for (int i = 0; i < nSamples; ++i)
		for (int j = 0; j < nDim; ++j) {
			Float sj = (i + (rng.UniformFloat())) * invNSamples;
			samples[nDim * i + j] = std::min(sj, OneMinusEpsilon);
		}
		for (int i = 0; i < nDim; ++i) {
		for (int j = 0; j < nSamples; ++j) {
			int other = j + rng.UniformUInt32(nSamples - j);
			std::swap(samples[nDim * j + i], samples[nDim * other + i]);
		}
	}
}

Point2f RejectionSampleDisk(RNG &rng) {
	Point2f p;
	do {
		p.x = 1 - 2 * rng.UniformFloat();
		p.y = 1 - 2 * rng.UniformFloat();
	} while (p.x * p.x + p.y * p.y > 1);
	return p;
}

Vector3f UniformSampleHemisphere(const Point2f &u) {
	Float z = u[0];
	Float r = std::sqrt(std::max((Float)0, (Float)1. - z * z));
	Float phi = 2 * Pi * u[1];
	return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformHemispherePdf() { return Inv2Pi; }

Vector3f UniformSampleSphere(const Point2f &u) {
	Float z = 1 - 2 * u[0];
	Float r = std::sqrt(std::max((Float)0, (Float)1 - z * z));
	Float phi = 2 * Pi * u[1];
	return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformSpherePdf() { return Inv4Pi; }

Point2f UniformSampleDisk(const Point2f &u) {
	Float r = std::sqrt(u[0]);
	Float theta = 2 * Pi * u[1];
	return Point2f(r * std::cos(theta), r * std::sin(theta));
}

Point2f ConcentricSampleDisk(const Point2f &u) {
	Point2f uOffset = 2.f * u - Vector2f(1, 1);
	if (uOffset.x == 0 && uOffset.y == 0) return Point2f(0, 0);
	Float theta, r;
	if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
		r = uOffset.x;
		theta = PiOver4 * (uOffset.y / uOffset.x);
	}
	else {
		r = uOffset.y;
		theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
	}
	return r * Point2f(std::cos(theta), std::sin(theta));
}

Float UniformConePdf(Float cosThetaMax) {
	return 1 / (2 * Pi * (1 - cosThetaMax));
}

Vector3f UniformSampleCone(const Point2f &u, Float cosThetaMax) {
	Float cosTheta = ((Float)1 - u[0]) + u[0] * cosThetaMax;
	Float sinTheta = std::sqrt((Float)1 - cosTheta * cosTheta);
	Float phi = u[1] * 2 * Pi;
	return Vector3f(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta,cosTheta);
}

Vector3f UniformSampleCone(const Point2f &u, Float cosThetaMax,
	const Vector3f &x, const Vector3f &y,
	const Vector3f &z) {
	Float cosTheta = Lerp(u[0], cosThetaMax, 1.f);
	Float sinTheta = std::sqrt((Float)1. - cosTheta * cosTheta);
	Float phi = u[1] * 2 * Pi;
	return std::cos(phi) * sinTheta * x + std::sin(phi) * sinTheta * y +cosTheta * z;
}

Point2f UniformSampleTriangle(const Point2f &u) {
	Float su0 = std::sqrt(u[0]);
	return Point2f(1 - su0, u[1] * su0);
}

Distribution2D::Distribution2D(const Float *func, int nu, int nv) {
	pConditionalV.reserve(nv);
	for (int v = 0; v < nv; ++v) {
		pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
	}
	std::vector<Float> marginalFunc;
	marginalFunc.reserve(nv);
	for (int v = 0; v < nv; ++v)
		marginalFunc.push_back(pConditionalV[v]->funcInt);
	pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
}
