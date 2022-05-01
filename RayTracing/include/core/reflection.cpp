#include "reflection.h"
#include "sampling.h"
#include "microfacet.h"
Fresnel::~Fresnel() {}

Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) {
	cosThetaI = Clamp(cosThetaI, -1, 1);
	// Potentially swap indices of refraction
	bool entering = cosThetaI > 0.f;
	if (!entering) {
		std::swap(etaI, etaT);
		cosThetaI = std::abs(cosThetaI);
	}

	// Compute _cosThetaT_ using Snell's law
	Float sinThetaI = std::sqrt(std::max((Float)0, 1 - cosThetaI * cosThetaI));
	Float sinThetaT = etaI / etaT * sinThetaI;

	// Handle total internal reflection
	if (sinThetaT >= 1) return 1;
	Float cosThetaT = std::sqrt(std::max((Float)0, 1 - sinThetaT * sinThetaT));
	Float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
		((etaT * cosThetaI) + (etaI * cosThetaT));
	Float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
		((etaI * cosThetaI) + (etaT * cosThetaT));
	return (Rparl * Rparl + Rperp * Rperp) / 2;
}

Spectrum FrConductor(Float cosThetaI, const Spectrum &etai,
	const Spectrum &etat, const Spectrum &k) {
	cosThetaI = Clamp(cosThetaI, -1, 1);
	Spectrum eta = etat / etai;
	Spectrum etak = k / etai;

	Float cosThetaI2 = cosThetaI * cosThetaI;
	Float sinThetaI2 = 1. - cosThetaI2;
	Spectrum eta2 = eta * eta;
	Spectrum etak2 = etak * etak;

	Spectrum t0 = eta2 - etak2 - sinThetaI2;
	Spectrum a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etak2);
	Spectrum t1 = a2plusb2 + cosThetaI2;
	Spectrum a = Sqrt(0.5f * (a2plusb2 + t0));
	Spectrum t2 = (Float)2 * cosThetaI * a;
	Spectrum Rs = (t1 - t2) / (t1 + t2);

	Spectrum t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
	Spectrum t4 = t2 * sinThetaI2;
	Spectrum Rp = Rs * (t3 - t4) / (t3 + t4);

	return 0.5 * (Rp + Rs);
}

Spectrum FresnelConductor::Evaluate(Float cosThetaI) const {
	return FrConductor(std::abs(cosThetaI), etaI, etaT, k);
}

std::string FresnelConductor::ToString() const {
	return std::string("[ FresnelConductor etaI: ") + etaI.ToString() +
		std::string(" etaT: ") + etaT.ToString() + std::string(" k: ") +
		k.ToString() + std::string(" ]");
}

Spectrum FresnelDielectric::Evaluate(Float cosThetaI) const {
	return FrDielectric(cosThetaI, etaI, etaT);
}

std::string FresnelDielectric::ToString() const {
	return StringPrintf("[ FrenselDielectric etaI: %f etaT: %f ]", etaI, etaT);
}

Spectrum SpecularReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
	const Point2f &sample, Float *pdf,
	BxDFType *sampledType) const {
	// Compute perfect specular reflection direction
	*wi = Vector3f(-wo.x, -wo.y, wo.z);
	*pdf = 1;
	return fresnel->Evaluate(CosTheta(*wi)) * R / AbsCosTheta(*wi);
}

std::string SpecularReflection::ToString() const {
	return std::string("[ SpecularReflection R: ") + R.ToString() +
		std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

Spectrum SpecularTransmission::Sample_f(const Vector3f &wo, Vector3f *wi,
	const Point2f &sample, Float *pdf,
	BxDFType *sampledType) const {
	// Figure out which �� is incident and which is transmitted
	bool entering = CosTheta(wo) > 0;
	Float etaI = entering ? etaA : etaB;
	Float etaT = entering ? etaB : etaA;

	// Compute ray direction for specular transmission
	if (!Refract(wo, Faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi))
		return 0;
	*pdf = 1;
	Spectrum ft = T * (Spectrum(1.) - fresnel.Evaluate(CosTheta(*wi)));
	// Account for non-symmetry with transmission to different medium
	if (mode == TransportMode::Radiance) ft *= (etaI * etaI) / (etaT * etaT);
	return ft / AbsCosTheta(*wi);
}

std::string SpecularTransmission::ToString() const {
	return std::string("[ SpecularTransmission: T: ") + T.ToString() +
		StringPrintf(" etaA: %f etaB: %f ", etaA, etaB) +
		std::string(" fresnel: ") + fresnel.ToString() +
		std::string(" mode : ") +
		(mode == TransportMode::Radiance ? std::string("RADIANCE")
			: std::string("IMPORTANCE")) +
		std::string(" ]");
}


Spectrum BxDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
	Float *pdf, BxDFType *sampledType) const {
	// Cosine-sample the hemisphere, flipping the direction if necessary
	*wi = CosineSampleHemisphere(u);
	if (wo.z < 0) wi->z *= -1;
	*pdf = Pdf(wo, *wi);
	return f(wo, *wi);
}

Float BxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
	return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}

Spectrum BxDF::rho(const Vector3f &w, int nSamples, const Point2f *u) const {
	Spectrum r(0.);
	for (int i = 0; i < nSamples; ++i) {
		// Estimate one term of Undefined control sequence \roman
		Vector3f wi;
		Float pdf = 0;
		Spectrum f = Sample_f(w, &wi, u[i], &pdf);
		if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
	}
	return r / nSamples;
}

Spectrum BxDF::rho(int nSamples, const Point2f *u1, const Point2f *u2) const {
	Spectrum r(0.f);
	for (int i = 0; i < nSamples; ++i) {
		// Estimate one term of Undefined control sequence \roman
		Vector3f wo, wi;
		wo = UniformSampleHemisphere(u1[i]);
		Float pdfo = UniformHemispherePdf(), pdfi = 0;
		Spectrum f = Sample_f(wo, &wi, u2[i], &pdfi);
		if (pdfi > 0)
			r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
	}
	return r / (Pi * nSamples);
}

Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const {
	return R * InvPi;
}


// BSDF Method Definitions
Spectrum BSDF::f(const Vector3f &woW, const Vector3f &wiW,
	BxDFType flags) const {
	//ProfilePhase pp(Prof::BSDFEvaluation);
	Vector3f wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
	if (wo.z == 0) return 0.;
	bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;
	Spectrum f(0.f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags) &&
			((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
			(!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
			f += bxdfs[i]->f(wo, wi);
	return f;
}

Spectrum BSDF::rho(int nSamples, const Point2f *samples1,
	const Point2f *samples2, BxDFType flags) const {
	Spectrum ret(0.f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags))
			ret += bxdfs[i]->rho(nSamples, samples1, samples2);
	return ret;
}

Spectrum BSDF::rho(const Vector3f &woWorld, int nSamples, const Point2f *samples,
	BxDFType flags) const {
	Vector3f wo = WorldToLocal(woWorld);
	Spectrum ret(0.f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags))
			ret += bxdfs[i]->rho(wo, nSamples, samples);
	return ret;
}

Spectrum BSDF::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,
	const Point2f &u, Float *pdf, BxDFType type,
	BxDFType *sampledType) const {
	//ProfilePhase pp(Prof::BSDFSampling);
	// Choose which _BxDF_ to sample
	int matchingComps = NumComponents(type);
	if (matchingComps == 0) {
		*pdf = 0;
		if (sampledType) *sampledType = BxDFType(0);
		return Spectrum(0);
	}
	int comp =
		std::min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

	// Get _BxDF_ pointer for chosen component
	BxDF *bxdf = nullptr;
	int count = comp;
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(type) && count-- == 0) {
			bxdf = bxdfs[i];
			break;
		}
	CHECK(bxdf != nullptr);
	//VLOG(2) << "BSDF::Sample_f chose comp = " << comp << " / matching = " <<
	//	matchingComps << ", bxdf: " << bxdf->ToString();

	// Remap _BxDF_ sample _u_ to [0,1)2
	Point2f uRemapped(std::min(u[0] * matchingComps - comp, OneMinusEpsilon),
		u[1]);

	// Sample chosen _BxDF_
	Vector3f wi, wo = WorldToLocal(woWorld);
	if (wo.z == 0) return 0.;
	*pdf = 0;
	if (sampledType) *sampledType = bxdf->type;
	Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);
	//VLOG(2) << "For wo = " << wo << ", sampled f = " << f << ", pdf = "
	//	<< *pdf << ", ratio = " << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.))
	//	<< ", wi = " << wi;
	if (*pdf == 0) {
		if (sampledType) *sampledType = BxDFType(0);
		return 0;
	}
	*wiWorld = LocalToWorld(wi);

	// Compute overall PDF with all matching _BxDF_s
	if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1)
		for (int i = 0; i < nBxDFs; ++i)
			if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type))
				*pdf += bxdfs[i]->Pdf(wo, wi);
	if (matchingComps > 1) *pdf /= matchingComps;

	// Compute value of BSDF for sampled direction
	if (!(bxdf->type & BSDF_SPECULAR)) {
		bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
		f = 0.;
		for (int i = 0; i < nBxDFs; ++i)
			if (bxdfs[i]->MatchesFlags(type) &&
				((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
				(!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
				f += bxdfs[i]->f(wo, wi);
	}
	//VLOG(2) << "Overall f = " << f << ", pdf = " << *pdf << ", ratio = "
	//	<< ((*pdf > 0) ? (f / *pdf) : Spectrum(0.));
	return f;
}

Float BSDF::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
	BxDFType flags) const {
	//ProfilePhase pp(Prof::BSDFPdf);
	if (nBxDFs == 0.f) return 0.f;
	Vector3f wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
	if (wo.z == 0) return 0.;
	Float pdf = 0.f;
	int matchingComps = 0;
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags)) {
			++matchingComps;
			pdf += bxdfs[i]->Pdf(wo, wi);
		}
	Float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
	return v;
}

//std::string BSDF::ToString() const {
//	std::string s = StringPrintf("[ BSDF eta: %f nBxDFs: %d", eta, nBxDFs);
//	for (int i = 0; i < nBxDFs; ++i)
//		s += StringPrintf("\n  bxdfs[%d]: ", i) + bxdfs[i]->ToString();
//	return s + std::string(" ]");
//}

std::string MicrofacetReflection::ToString() const {
	return std::string("[ MicrofacetReflection R: ") + R.ToString() +
		std::string(" distribution: ") + distribution->ToString() +
		std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

Spectrum MicrofacetReflection::f(const Vector3f &wo, const Vector3f &wi) const {
	Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
	Vector3f wh = wi + wo;
	// Handle degenerate cases for microfacet reflection
	if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
	wh = Normalize(wh);
	// For the Fresnel call, make sure that wh is in the same hemisphere
	// as the surface normal, so that TIR is handled correctly.
	Spectrum F = fresnel->Evaluate(Dot(wi, Faceforward(wh, Vector3f(0, 0, 1))));
	return R * distribution->D(wh) * distribution->G(wo, wi) * F /
		(4 * cosThetaI * cosThetaO);
}

Spectrum MicrofacetReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
	const Point2f &u, Float *pdf,
	BxDFType *sampledType) const {
	// Sample microfacet orientation Undefined control sequence \wh and reflected direction Undefined control sequence \wi
	if (wo.z == 0) return 0.;
	Vector3f wh = distribution->Sample_wh(wo, u);
	if (Dot(wo, wh) < 0) return 0.;   // Should be rare
	*wi = Reflect(wo, wh);
	if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

	// Compute PDF of _wi_ for microfacet reflection
	*pdf = distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
	return f(wo, *wi);
}

Float MicrofacetReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
	if (!SameHemisphere(wo, wi)) return 0;
	Vector3f wh = Normalize(wo + wi);
	return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
}

Spectrum MetalRoughnessReflection::f(const Vector3f &wo, const Vector3f &wi) const {
	Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
	Vector3f wh = wi + wo;
	// Handle degenerate cases for microfacet reflection
	if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
	if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
	wh = Normalize(wh);
	// For the Fresnel call, make sure that wh is in the same hemisphere
	// as the surface normal, so that TIR is handled correctly.
	//std::cout << wh << std::endl;
	Float cos = abs(Dot(wi, wh));
	//Float cosThetaI = std::max(Dot(wi, Faceforward(wh, Vector3f(0, 0, 1))),0.f);
	Spectrum F0(0.04);
	F0 = (1-metallic) * F0 + metallic*c;
	Spectrum F= F0 + (Spectrum(1.0) - F0) * pow(1.0 - cos, 5.0);
	Spectrum specular = distribution->D(wh) * distribution->G(wo, wi) * F /
		(4 * cosThetaI * cosThetaO+0.001);
	//Float D = DistributionGGX(N, wh, roughness);
	//Float G = GeometrySmith(N, wo, wi, roughness);
	//std::cout << 's' << ' ' << distribution->D(wh) << ' ' 	<< distribution->G(wo, wi) << ' '<<F<<std::endl;
	Spectrum kS = F;
	Spectrum kD = Spectrum(1.0) - kS;
	kD *= 1.0 - metallic;
	//Spectrum F = fresnel->Evaluate(Dot(wi, Faceforward(wh, Vector3f(0, 0, 1))));
	//std::cout << 's'<<' '<<specular << ' ' << kD * c*InvPi << std::endl;
	return specular+kD*c*InvPi;
}

Spectrum MetalRoughnessReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
	const Point2f &u, Float *pdf,
	BxDFType *sampledType) const {
	// Sample microfacet orientation Undefined control sequence \wh and reflected direction Undefined control sequence \wi
	if (wo.z == 0) return 0.;
	Vector3f wh = distribution->Sample_wh(wo, u);
	if (Dot(wo, wh) < 0) return 0.;   // Should be rare
	*wi = Reflect(wo, wh);
	if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

	// Compute PDF of _wi_ for microfacet reflection
	*pdf = distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
	return f(wo, *wi);
}

Float MetalRoughnessReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
	if (!SameHemisphere(wo, wi)) return 0;
	Vector3f wh = Normalize(wo + wi);
	return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
}

std::string MicrofacetTransmission::ToString() const {
	return std::string("[ MicrofacetTransmission T: ") + T.ToString() +
		std::string(" distribution: ") + distribution->ToString() +
		StringPrintf(" etaA: %f etaB: %f", etaA, etaB) +
		std::string(" fresnel: ") + fresnel.ToString() +
		std::string(" mode : ") +
		(mode == TransportMode::Radiance ? std::string("RADIANCE")
			: std::string("IMPORTANCE")) +
		std::string(" ]");
}

Spectrum MicrofacetTransmission::f(const Vector3f &wo,
	const Vector3f &wi) const {
	if (SameHemisphere(wo, wi)) return 0;  // transmission only

	Float cosThetaO = CosTheta(wo);
	Float cosThetaI = CosTheta(wi);
	if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0);

	// Compute Undefined control sequence \wh from Undefined control sequence \wo and Undefined control sequence \wi for microfacet transmission
	Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
	Vector3f wh = Normalize(wo + wi * eta);
	if (wh.z < 0) wh = -wh;

	// Same side?
	if (Dot(wo, wh) * Dot(wi, wh) > 0) return Spectrum(0);

	Spectrum F = fresnel.Evaluate(Dot(wo, wh));

	Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
	Float factor = (mode == TransportMode::Radiance) ? (1 / eta) : 1;

	return (Spectrum(1.f) - F) * T *
		std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta *
			AbsDot(wi, wh) * AbsDot(wo, wh) * factor * factor /
			(cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
}

Spectrum MicrofacetTransmission::Sample_f(const Vector3f &wo, Vector3f *wi,
	const Point2f &u, Float *pdf,
	BxDFType *sampledType) const {
	if (wo.z == 0) return 0.;
	Vector3f wh = distribution->Sample_wh(wo, u);
	if (Dot(wo, wh) < 0) return 0.;  // Should be rare

	Float eta = CosTheta(wo) > 0 ? (etaA / etaB) : (etaB / etaA);
	if (!Refract(wo, (Normal3f)wh, eta, wi)) return 0;
	*pdf = Pdf(wo, *wi);
	return f(wo, *wi);
}

Float MicrofacetTransmission::Pdf(const Vector3f &wo,
	const Vector3f &wi) const {
	if (SameHemisphere(wo, wi)) return 0;
	// Compute Undefined control sequence \wh from Undefined control sequence \wo and Undefined control sequence \wi for microfacet transmission
	Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
	Vector3f wh = Normalize(wo + wi * eta);

	if (Dot(wo, wh) * Dot(wi, wh) > 0) return 0;

	// Compute change of variables _dwh\_dwi_ for microfacet transmission
	Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
	Float dwh_dwi =
		std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));
	return distribution->Pdf(wo, wh) * dwh_dwi;
}

//FresnelBlend::FresnelBlend(const Spectrum &Rd, const Spectrum &Rs,
//	MicrofacetDistribution *distribution)
//	: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
//	Rd(Rd),
//	Rs(Rs),
//	distribution(distribution) {}
//Spectrum FresnelBlend::f(const Vector3f &wo, const Vector3f &wi) const {
//	auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
//	Spectrum diffuse = (28.f / (23.f * Pi)) * Rd * (Spectrum(1.f) - Rs) *
//		(1 - pow5(1 - .5f * AbsCosTheta(wi))) *
//		(1 - pow5(1 - .5f * AbsCosTheta(wo)));
//	Vector3f wh = wi + wo;
//	if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0);
//	wh = Normalize(wh);
//	Spectrum specular =
//		distribution->D(wh) /
//		(4 * AbsDot(wi, wh) * std::max(AbsCosTheta(wi), AbsCosTheta(wo))) *
//		SchlickFresnel(Dot(wi, wh));
//	return diffuse + specular;
//}
//
//std::string FresnelBlend::ToString() const {
//	return std::string("[ FresnelBlend Rd: ") + Rd.ToString() +
//		std::string(" Rs: ") + Rs.ToString() +
//		std::string(" distribution: ") + distribution->ToString() +
//		std::string(" ]");
//}

//Spectrum FresnelBlend::Sample_f(const Vector3f &wo, Vector3f *wi,
//	const Point2f &uOrig, Float *pdf,
//	BxDFType *sampledType) const {
//	Point2f u = uOrig;
//	if (u[0] < .5) {
//		u[0] = std::min(2 * u[0], OneMinusEpsilon);
//		// Cosine-sample the hemisphere, flipping the direction if necessary
//		*wi = CosineSampleHemisphere(u);
//		if (wo.z < 0) wi->z *= -1;
//	}
//	else {
//		u[0] = std::min(2 * (u[0] - .5f), OneMinusEpsilon);
//		// Sample microfacet orientation Undefined control sequence \wh and reflected direction Undefined control sequence \wi
//		Vector3f wh = distribution->Sample_wh(wo, u);
//		*wi = Reflect(wo, wh);
//		if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
//	}
//	*pdf = Pdf(wo, *wi);
//	return f(wo, *wi);
//}
//
//Float FresnelBlend::Pdf(const Vector3f &wo, const Vector3f &wi) const {
//	if (!SameHemisphere(wo, wi)) return 0;
//	Vector3f wh = Normalize(wo + wi);
//	Float pdf_wh = distribution->Pdf(wo, wh);
//	return .5f * (AbsCosTheta(wi) * InvPi + pdf_wh / (4 * Dot(wo, wh)));
//}
