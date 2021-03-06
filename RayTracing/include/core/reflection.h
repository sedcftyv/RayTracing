#pragma once
#ifndef PBRT_CORE_REFLECTION_H
#define PBRT_CORE_REFLECTION_H

#include "pbrt.h"
#include "geometry.h"
#include "shape.h"
#include "spectrum.h"
#include "interaction.h"

inline Float CosTheta(const Vector3f &w) { return w.z; }
inline Float Cos2Theta(const Vector3f &w) { return w.z * w.z; }
inline Float AbsCosTheta(const Vector3f &w) { return std::abs(w.z); }
inline Float Sin2Theta(const Vector3f &w) {
	return std::max((Float)0, (Float)1 - Cos2Theta(w));
}

inline Float SinTheta(const Vector3f &w) { return std::sqrt(Sin2Theta(w)); }

inline Float TanTheta(const Vector3f &w) { return SinTheta(w) / CosTheta(w); }

inline Float Tan2Theta(const Vector3f &w) {
	return Sin2Theta(w) / Cos2Theta(w);
}

inline Float CosPhi(const Vector3f &w) {
	Float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
}

inline Float SinPhi(const Vector3f &w) {
	Float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
}

inline Float Cos2Phi(const Vector3f &w) { return CosPhi(w) * CosPhi(w); }

inline Float Sin2Phi(const Vector3f &w) { return SinPhi(w) * SinPhi(w); }

inline Float CosDPhi(const Vector3f &wa, const Vector3f &wb) {
	Float waxy = wa.x * wa.x + wa.y * wa.y;
	Float wbxy = wb.x * wb.x + wb.y * wb.y;
	if (waxy == 0 || wbxy == 0)
		return 1;
	return Clamp((wa.x * wb.x + wa.y * wb.y) / std::sqrt(waxy * wbxy), -1, 1);
}

inline Vector3f Reflect(const Vector3f &wo, const Vector3f &n) {
	return -wo + 2 * Dot(wo, n) * n;
}

inline bool Refract(const Vector3f &wi, const Normal3f &n, Float eta,
	Vector3f *wt) {
	Float cosThetaI = Dot(n, wi);
	Float sin2ThetaI = std::max(Float(0), Float(1 - cosThetaI * cosThetaI));
	Float sin2ThetaT = eta * eta * sin2ThetaI;
	if (sin2ThetaT >= 1) return false;
	Float cosThetaT = std::sqrt(1 - sin2ThetaT);
	*wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
	return true;
}

inline bool SameHemisphere(const Vector3f &w, const Vector3f &wp) {
	return w.z * wp.z > 0;
}

inline bool SameHemisphere(const Vector3f &w, const Normal3f &wp) {
	return w.z * wp.z > 0;
}

enum BxDFType {
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE = 1 << 2,
	BSDF_GLOSSY = 1 << 3,
	BSDF_SPECULAR = 1 << 4,
	BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
	BSDF_TRANSMISSION,
};

class BxDF {
public:
	virtual ~BxDF() {}
	BxDF(BxDFType type) : type(type) {}
	bool MatchesFlags(BxDFType t) const { return (type & t) == type; }
	virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;
	virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi,
		const Point2f &sample, Float *pdf,
		BxDFType *sampledType = nullptr) const;
	virtual Spectrum rho(const Vector3f &wo, int nSamples,
		const Point2f *samples) const;
	virtual Spectrum rho(int nSamples, const Point2f *samples1,
		const Point2f *samples2) const;
	virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
	
	const BxDFType type;
};

class Fresnel {
public:
	virtual ~Fresnel();
	virtual Spectrum Evaluate(Float cosI) const = 0;
	virtual std::string ToString() const = 0;
};

inline std::ostream &operator<<(std::ostream &os, const Fresnel &f) {
	os << f.ToString();
	return os;
}

class FresnelConductor : public Fresnel {
public:
	Spectrum Evaluate(Float cosThetaI) const;
	FresnelConductor(const Spectrum &etaI, const Spectrum &etaT,const Spectrum &k)
		: etaI(etaI), etaT(etaT), k(k) {}
	std::string ToString() const;
private:
	Spectrum etaI, etaT, k;
};

class FresnelDielectric : public Fresnel {
public:
	Spectrum Evaluate(Float cosThetaI) const;
	FresnelDielectric(Float etaI, Float etaT) : etaI(etaI), etaT(etaT) {}
	std::string ToString() const;

private:
	Float etaI, etaT;
};

class FresnelNoOp : public Fresnel {
public:
	Spectrum Evaluate(Float) const { return Spectrum(1.); }
	std::string ToString() const { return "[ FresnelNoOp ]"; }
};

class BSDF {
public:
	BSDF(const SurfaceInteraction &si, Float eta = 1)
	: eta(eta),ns(si.shading.n),ng(si.n),ss(Normalize(si.shading.dpdu)),
		ts(Cross(ns, ss)) {}
	~BSDF() {
		for (int i = 0; i < nBxDFs; ++i){
			delete bxdfs[i];
		}
	}
	void Add(BxDF *b) {
		bxdfs[nBxDFs++] = b;
	}
	int NumComponents(BxDFType flags = BSDF_ALL) const;
	Vector3f WorldToLocal(const Vector3f &v) const {
		return Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
	}
	Vector3f LocalToWorld(const Vector3f &v) const {
		return Vector3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
			ss.y * v.x + ts.y * v.y + ns.y * v.z,
			ss.z * v.x + ts.z * v.y + ns.z * v.z);
	}
	Spectrum f(const Vector3f &woW, const Vector3f &wiW,BxDFType flags = BSDF_ALL) const;
	Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2,BxDFType flags = BSDF_ALL) const;
	Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples,BxDFType flags = BSDF_ALL) const;
	Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,Float *pdf, BxDFType type = BSDF_ALL,BxDFType *sampledType = nullptr) const;
	Float Pdf(const Vector3f &wo, const Vector3f &wi,BxDFType flags = BSDF_ALL) const;const Float eta;

private:
	
	const Normal3f ns, ng;
	const Vector3f ss, ts;
	int nBxDFs = 0;
	static PBRT_CONSTEXPR int MaxBxDFs = 8;
	BxDF *bxdfs[MaxBxDFs];
	friend class MixMaterial;
};

inline int BSDF::NumComponents(BxDFType flags) const {
	int num = 0;
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags)) ++num;
	return num;
}


class LambertianReflection : public BxDF {
public:
	LambertianReflection(const Spectrum &R)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
	Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
	Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
	Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
	
private:
	const Spectrum R;
};

class SpecularReflection : public BxDF {
public:
	SpecularReflection(const Spectrum &R, Fresnel *fresnel)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),R(R),fresnel(fresnel) {}
	~SpecularReflection(){
		if (fresnel!=nullptr)
			delete fresnel;
	}
	Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
		return Spectrum(0.f);
	}
	Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
		Float *pdf, BxDFType *sampledType) const;
	Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
	std::string ToString() const;

private:
	const Spectrum R;
	const Fresnel *fresnel;
};

class SpecularTransmission : public BxDF {
public:
	SpecularTransmission(const Spectrum &T, Float etaA, Float etaB,
	TransportMode mode)
	: BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
	T(T),
	etaA(etaA),
	etaB(etaB),
	fresnel(etaA, etaB),
	mode(mode) {}
	Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
		return Spectrum(0.f);
	}
	Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,Float *pdf, BxDFType *sampledType) const;
	Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
	std::string ToString() const;

private:
	const Spectrum T;
	const Float etaA, etaB;
	const FresnelDielectric fresnel;
	const TransportMode mode;
};

class MicrofacetReflection : public BxDF {
public:
	MicrofacetReflection(const Spectrum &R,
	MicrofacetDistribution *distribution, Fresnel *fresnel)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),R(R),distribution(distribution),
	fresnel(fresnel) {}
	~MicrofacetReflection()
	{
		if (fresnel != nullptr)
			delete fresnel;
		if (distribution != nullptr)
			delete distribution;
	}
	Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
	Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
		Float *pdf, BxDFType *sampledType) const;
	Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
	std::string ToString() const;

private:
	const Spectrum R;
	const MicrofacetDistribution *distribution;
	const Fresnel *fresnel;
};

class MetalRoughnessReflection : public BxDF {
public:
	MetalRoughnessReflection(const Spectrum &c,const Float &metallic, const Float &roughness,MicrofacetDistribution *distribution)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
		c(c), metallic(metallic), roughness(roughness),
		distribution(distribution){}
	~MetalRoughnessReflection(){
		if (distribution != nullptr)
			delete distribution;
	}
	Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
	Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
		Float *pdf, BxDFType *sampledType) const;
	Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
	
private:
	const Spectrum c;
	const Float metallic,roughness;
	const MicrofacetDistribution *distribution;
};

class MicrofacetTransmission : public BxDF {
public:
	MicrofacetTransmission(const Spectrum &T,
	MicrofacetDistribution *distribution, Float etaA,Float etaB, TransportMode mode)
		: BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
		T(T),distribution(distribution),etaA(etaA),etaB(etaB),fresnel(etaA, etaB),mode(mode) {}
	Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
	Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
		Float *pdf, BxDFType *sampledType) const;
	Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
	std::string ToString() const;

private:
	const Spectrum T;
	const MicrofacetDistribution *distribution;
	const Float etaA, etaB;
	const FresnelDielectric fresnel;
	const TransportMode mode;
};




#endif