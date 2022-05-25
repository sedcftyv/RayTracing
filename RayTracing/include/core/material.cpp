#include "material.h"
#include "reflection.h"
#include "core/texture.h"
#include "microfacet.h"
void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,
	bool allowMultipleLobes) const {
	
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);

	// Evaluate textures for _MatteMaterial_ material and allocate BRDF
	si->bsdf = new BSDF(*si);
	Spectrum r = Kd->Evaluate(*si).Clamp();
	//Float sig = Clamp(sigma->Evaluate(*si), 0, 90);
	//if (!r.IsBlack()) {
	//	if (sig == 0)
	si->bsdf->Add(new LambertianReflection(r));
		//else
		//	si->bsdf->Add(new OrenNayar(r));
	//}
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

void PlasticMaterial::ComputeScatteringFunctions(
	SurfaceInteraction *si,  TransportMode mode,
	bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);
	si->bsdf = new BSDF(*si);
	// Initialize diffuse component of plastic material
	Spectrum kd = Kd->Evaluate(*si).Clamp();
	if (!kd.IsBlack())
		si->bsdf->Add(new LambertianReflection(kd));

	// Initialize specular component of plastic material
	Spectrum ks = Ks->Evaluate(*si).Clamp();
	if (!ks.IsBlack()) {
		Fresnel *fresnel = new FresnelDielectric(1.5f, 1.f);
		// Create microfacet distribution _distrib_ for plastic material
		Float rough = roughness->Evaluate(*si);
		//Spectrum tmp = roughness->Evaluate(*si).Clamp();
		//Float rough = 0.f;
		if (remapRoughness)
			rough = TrowbridgeReitzDistribution::RoughnessToAlpha(rough);
		MicrofacetDistribution *distrib =
			new TrowbridgeReitzDistribution(rough, rough);
		BxDF *spec =
			//new MicrofacetReflection(1., distrib, fresnel);
			new MicrofacetReflection(ks, distrib, fresnel);
		si->bsdf->Add(spec);
	}
}

void MetalMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,
	bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);
	si->bsdf = new BSDF(*si);

	Float uRough =
		uRoughness ? uRoughness->Evaluate(*si) : roughness->Evaluate(*si);
	Float vRough =
		vRoughness ? vRoughness->Evaluate(*si) : roughness->Evaluate(*si);
	if (remapRoughness) {
		uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(uRough);
		vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(vRough);
	}
	Fresnel *frMf = new FresnelConductor(1., eta->Evaluate(*si),
		k->Evaluate(*si));
	MicrofacetDistribution *distrib =
		new TrowbridgeReitzDistribution(uRough, vRough);
	si->bsdf->Add(new MicrofacetReflection(1., distrib, frMf));
}

void MetalRoughnessMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,bool allowMultipleLobes) const {
	// Perform bump mapping with _bumpMap_, if present
	//if (bumpMap) Bump(bumpMap, si);
	si->bsdf = new BSDF(*si);

	Spectrum c = base_color->Evaluate(*si).Clamp();
	if (metalroughness == nullptr)
	{
		si->bsdf->Add(new LambertianReflection(c));
	}
	else
	{
		Spectrum tmp = metalroughness->Evaluate(*si).Clamp();
		Float metallic = tmp[2], roughness = std::max(0.05f,tmp[1]);
		//metallic = 1.0; roughness = 0.05;
		MicrofacetDistribution *distrib =
			new TrowbridgeReitzDistribution(roughness, roughness);
			//new BeckmannDistribution(roughness, roughness);
			//new PBRDistribution(roughness,si->n);
		si->bsdf->Add(new MetalRoughnessReflection(c, metallic, roughness, distrib));
	}

}