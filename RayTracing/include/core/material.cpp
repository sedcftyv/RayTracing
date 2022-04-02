#include "material.h"
#include "reflection.h"
#include "core/texture.h"
void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,
	bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);

	// Evaluate textures for _MatteMaterial_ material and allocate BRDF
	si->bsdf = new BSDF(*si);
	Spectrum r = Kd->Evaluate(*si).Clamp();
	Float sig = Clamp(sigma->Evaluate(*si), 0, 90);
	if (!r.IsBlack()) {
		if (sig == 0)
			si->bsdf->Add(new LambertianReflection(r));
		//else
		//	si->bsdf->Add(new OrenNayar(r));
	}
}

void MirrorMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,
	bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);
	si->bsdf = new BSDF(*si);
	Spectrum R = Kr->Evaluate(*si).Clamp();
	if (!R.IsBlack())
		si->bsdf->Add(new SpecularReflection(
			R, new FresnelNoOp()));
}

void GlassMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);
	Float eta = index->Evaluate(*si);
	//Float urough = uRoughness->Evaluate(*si);
	//Float vrough = vRoughness->Evaluate(*si);
	Spectrum R = Kr->Evaluate(*si).Clamp();
	Spectrum T = Kt->Evaluate(*si).Clamp();
	// Initialize _bsdf_ for smooth or rough dielectric
	si->bsdf = new BSDF(*si, eta);
	if (!R.IsBlack() && !T.IsBlack())
		si->bsdf->Add(new SpecularTransmission(T,1.f,eta,mode));
	//if (R.IsBlack() && T.IsBlack()) return;

	//bool isSpecular = urough == 0 && vrough == 0;
	//if (isSpecular && allowMultipleLobes) {
	//	si->bsdf->Add(
	//		ARENA_ALLOC(arena, FresnelSpecular)(R, T, 1.f, eta, mode));
	//}
	//else {
	//	if (remapRoughness) {
	//		urough = TrowbridgeReitzDistribution::RoughnessToAlpha(urough);
	//		vrough = TrowbridgeReitzDistribution::RoughnessToAlpha(vrough);
	//	}
	//	MicrofacetDistribution *distrib =
	//		isSpecular ? nullptr
	//		: ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(
	//			urough, vrough);
	//	if (!R.IsBlack()) {
	//		Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);
	//		if (isSpecular)
	//			si->bsdf->Add(
	//				ARENA_ALLOC(arena, SpecularReflection)(R, fresnel));
	//		else
	//			si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetReflection)(
	//				R, distrib, fresnel));
	//	}
	//	if (!T.IsBlack()) {
	//		if (isSpecular)
	//			si->bsdf->Add(ARENA_ALLOC(arena, SpecularTransmission)(
	//				T, 1.f, eta, mode));
	//		else
	//			si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
	//				T, distrib, 1.f, eta, mode));
	//	}
	//}
}
