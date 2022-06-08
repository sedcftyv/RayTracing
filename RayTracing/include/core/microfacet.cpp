#include "core/microfacet.h"
#include "reflection.h"

static void BeckmannSample11(Float cosThetaI, Float U1, Float U2,
	Float *slope_x, Float *slope_y) {
	/* Special case (normal incidence) */
	if (cosThetaI > .9999) {
		Float r = std::sqrt(-std::log(1.0f - U1));
		Float sinPhi = std::sin(2 * Pi * U2);
		Float cosPhi = std::cos(2 * Pi * U2);
		*slope_x = r * cosPhi;
		*slope_y = r * sinPhi;
		return;
	}
	/* The original inversion routine from the paper contained
	   discontinuities, which causes issues for QMC integration
	   and techniques like Kelemen-style MLT. The following code
	   performs a numerical inversion with better behavior */
	Float sinThetaI =
		std::sqrt(std::max((Float)0, (Float)1 - cosThetaI * cosThetaI));
	Float tanThetaI = sinThetaI / cosThetaI;
	Float cotThetaI = 1 / tanThetaI;
	/* Search interval -- everything is parameterized
	   in the Erf() domain */
	Float a = -1, c = Erf(cotThetaI);
	Float sample_x = std::max(U1, (Float)1e-6f);
	/* Start with a good initial guess */
	/* We can do better (inverse of an approximation computed in
	 * Mathematica) */
	Float thetaI = std::acos(cosThetaI);
	Float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
	Float b = c - (1 + c) * std::pow(1 - sample_x, fit);
	/* Normalization factor for the CDF */
	static const Float SQRT_PI_INV = 1.f / std::sqrt(Pi);
	Float normalization =1 /(1 + c + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));
	int it = 0;
	while (++it < 10) {
		/* Bisection criterion -- the oddly-looking
		   Boolean expression are intentional to check
		   for NaNs at little additional cost */
		if (!(b >= a && b <= c)) b = 0.5f * (a + c);
		/* Evaluate the CDF and its derivative
		   (i.e. the density function) */
		Float invErf = ErfInv(b);
		Float value =normalization *(1 + b + SQRT_PI_INV * tanThetaI * std::exp(-invErf * invErf)) -sample_x;
		Float derivative = normalization * (1 - invErf * tanThetaI);
		if (std::abs(value) < 1e-5f) break;
		/* Update bisection intervals */
		if (value > 0)
			c = b;
		else
			a = b;
		b -= value / derivative;
	}
	/* Now convert back into a slope value */
	*slope_x = ErfInv(b);
	/* Simulate Y component */
	*slope_y = ErfInv(2.0f * std::max(U2, (Float)1e-6f) - 1.0f);

}
static Vector3f BeckmannSample(const Vector3f &wi, Float alpha_x, Float alpha_y,
	Float U1, Float U2) {
	Vector3f wiStretched =
	Normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));
	Float slope_x, slope_y;
	BeckmannSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);
	Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
	slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
	slope_x = tmp;
	slope_x = alpha_x * slope_x;
	slope_y = alpha_y * slope_y;
	return Normalize(Vector3f(-slope_x, -slope_y, 1.f));
}

MicrofacetDistribution::~MicrofacetDistribution() {}

Float BeckmannDistribution::D(const Vector3f &wh) const {
	Float tan2Theta = Tan2Theta(wh);
	if (std::isinf(tan2Theta)) return 0.;
	Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
	return std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +Sin2Phi(wh) / (alphay * alphay))) /(Pi * alphax * alphay * cos4Theta);
}

Float TrowbridgeReitzDistribution::D(const Vector3f &wh) const {
	Float tan2Theta = Tan2Theta(wh);
	if (std::isinf(tan2Theta)) return 0.;
	const Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
	Float e =(Cos2Phi(wh) / (alphax * alphax) + Sin2Phi(wh) / (alphay * alphay)) *tan2Theta;
	return 1 / (Pi * alphax * alphay * cos4Theta * (1 + e) * (1 + e));
}

Float BeckmannDistribution::Lambda(const Vector3f &w) const {
	Float absTanTheta = std::abs(TanTheta(w));
	if (std::isinf(absTanTheta)) return 0.;
	Float alpha =std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
	Float a = 1 / (alpha * absTanTheta);
	if (a >= 1.6f) return 0;
	return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

Float TrowbridgeReitzDistribution::Lambda(const Vector3f &w) const {
	Float absTanTheta = std::abs(TanTheta(w));
	if (std::isinf(absTanTheta)) return 0.;
	Float alpha =std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
	Float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
	return (-1 + std::sqrt(1.f + alpha2Tan2Theta)) / 2;
}

std::string BeckmannDistribution::ToString() const {
	return StringPrintf("[ BeckmannDistribution alphax: %f alphay: %f ]",alphax, alphay);
}

std::string TrowbridgeReitzDistribution::ToString() const {
	return StringPrintf("[ TrowbridgeReitzDistribution alphax: %f alphay: %f ]",alphax, alphay);
}

Vector3f BeckmannDistribution::Sample_wh(const Vector3f &wo,
	const Point2f &u) const {
	if (!sampleVisibleArea) {	
		Float tan2Theta, phi;
		if (alphax == alphay) {
			Float logSample = std::log(1 - u[0]);
			tan2Theta = -alphax * alphax * logSample;
			phi = u[1] * 2 * Pi;
		}
		else {
			Float logSample = std::log(1 - u[0]);
			phi = std::atan(alphay / alphax *std::tan(2 * Pi * u[1] + 0.5f * Pi));
			if (u[1] > 0.5f) phi += Pi;
			Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
			Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
			tan2Theta = -logSample /(cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
		}
		Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
		Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
		Vector3f wh = SphericalDirection(sinTheta, cosTheta, phi);
		if (!SameHemisphere(wo, wh)) wh = -wh;
		return wh;
	}
	else {
		Vector3f wh;
		bool flip = wo.z < 0;
		wh = BeckmannSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) wh = -wh;
		return wh;
	}
}

static void TrowbridgeReitzSample11(Float cosTheta, Float U1, Float U2,
	Float *slope_x, Float *slope_y) {
		if (cosTheta > .9999) {
		Float r = sqrt(U1 / (1 - U1));
		Float phi = 6.28318530718 * U2;
		*slope_x = r * cos(phi);
		*slope_y = r * sin(phi);
		return;
	}

	Float sinTheta =std::sqrt(std::max((Float)0, (Float)1 - cosTheta * cosTheta));
	Float tanTheta = sinTheta / cosTheta;
	Float a = 1 / tanTheta;
	Float G1 = 2 / (1 + std::sqrt(1.f + 1.f / (a * a)));
	Float A = 2 * U1 / G1 - 1;
	Float tmp = 1.f / (A * A - 1.f);
	if (tmp > 1e10) tmp = 1e10;
	Float B = tanTheta;
	Float D = std::sqrt(std::max(Float(B * B * tmp * tmp - (A * A - B * B) * tmp), Float(0)));
	Float slope_x_1 = B * tmp - D;
	Float slope_x_2 = B * tmp + D;
	*slope_x = (A < 0 || slope_x_2 > 1.f / tanTheta) ? slope_x_1 : slope_x_2;
	Float S;
	if (U2 > 0.5f) {
		S = 1.f;
		U2 = 2.f * (U2 - .5f);
	}
	else {
		S = -1.f;
		U2 = 2.f * (.5f - U2);
	}
	Float z =(U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
		(U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
		*slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);

}

static Vector3f TrowbridgeReitzSample(const Vector3f &wi, Float alpha_x,
	Float alpha_y, Float U1, Float U2) {
	Vector3f wiStretched =
	Normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));
	Float slope_x, slope_y;
	TrowbridgeReitzSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);
	Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
	slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
	slope_x = tmp;
	slope_x = alpha_x * slope_x;
	slope_y = alpha_y * slope_y;
	return Normalize(Vector3f(-slope_x, -slope_y, 1.));
}

Vector3f TrowbridgeReitzDistribution::Sample_wh(const Vector3f &wo,
	const Point2f &u) const {
	Vector3f wh;
	if (!sampleVisibleArea) {
		Float cosTheta = 0, phi = (2 * Pi) * u[1];
		if (alphax == alphay) {
			Float tanTheta2 = alphax * alphax * u[0] / (1.0f - u[0]);
			cosTheta = 1 / std::sqrt(1 + tanTheta2);
		}
		else {
			phi =std::atan(alphay / alphax * std::tan(2 * Pi * u[1] + .5f * Pi));
			if (u[1] > .5f) phi += Pi;
			Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
			const Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
			const Float alpha2 =1 / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
			Float tanTheta2 = alpha2 * u[0] / (1 - u[0]);
			cosTheta = 1 / std::sqrt(1 + tanTheta2);
		}
		Float sinTheta =std::sqrt(std::max((Float)0., (Float)1. - cosTheta * cosTheta));
		wh = SphericalDirection(sinTheta, cosTheta, phi);
		if (!SameHemisphere(wo, wh)) wh = -wh;
	}
	else {
		bool flip = wo.z < 0;
		wh = TrowbridgeReitzSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) wh = -wh;
	}
	return wh;
}

Float MicrofacetDistribution::Pdf(const Vector3f &wo,
	const Vector3f &wh) const {
	if (sampleVisibleArea)
		return D(wh) * G1(wo) * AbsDot(wo, wh) / AbsCosTheta(wo);
	else
		return D(wh) * AbsCosTheta(wh);
}

std::string PBRDistribution::ToString() const {
	return StringPrintf("[ PBRDistribution roughness: %f ]",roughness);
}

Float PBRDistribution::Lambda(const Vector3f &w) const {
	return 0.f;
}

Vector3f PBRDistribution::Sample_wh(const Vector3f &wo,
	const Point2f &u) const {
	Vector3f wh;
	bool flip = wo.z < 0;
	wh = TrowbridgeReitzSample(flip ? -wo : wo, roughness, roughness, u[0], u[1]);
	if (flip) wh = -wh;
	return wh;
}

Float PBRDistribution::D(const Vector3f &wh) const {
	Float a = roughness * roughness;
	Float a2 = a * a;
	Float NdotH = std::max(Dot(N, wh), 0.0f);
	Float NdotH2 = NdotH * NdotH;
	Float nom = a2;
	Float denom = (NdotH2 * (a2 - 1.0) + 1.0);
	denom = Pi * denom * denom;
	return nom / denom;
}

inline Float GeometrySchlickGGX(Float NdotV, Float roughness)
{
	Float r = (roughness + 1.0);
	Float k = (r*r) / 8.0;
	Float nom = NdotV;
	Float denom = NdotV * (1.0 - k) + k;
	return nom / denom;
}

Float PBRDistribution::G(const Vector3f &wo, const Vector3f &wi) const {
	Float NdotV = std::max(Dot(N, wo), 0.0f);
	Float NdotL = std::max(Dot(N, wi), 0.0f);
	Float ggx2 = GeometrySchlickGGX(NdotV, roughness);
	Float ggx1 = GeometrySchlickGGX(NdotL, roughness);
	return ggx1 * ggx2;
}