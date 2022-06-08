#include "material.h"
#include "reflection.h"
#include "core/texture.h"
#include "microfacet.h"
void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,TransportMode mode,
	bool allowMultipleLobes) const {
	si->bsdf = new BSDF(*si);
	Spectrum r = Kd->Evaluate(*si).Clamp();
	si->bsdf->Add(new LambertianReflection(r));
}

void MirrorMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,
	bool allowMultipleLobes) const {
	si->bsdf = new BSDF(*si);
	Spectrum R = Kr->Evaluate(*si).Clamp();
	if (!R.IsBlack())
		si->bsdf->Add(new SpecularReflection(
	R, new FresnelNoOp()));
}

void GlassMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,bool allowMultipleLobes) const {
	Float eta = index->Evaluate(*si);
	Spectrum R = Kr->Evaluate(*si).Clamp();
	Spectrum T = Kt->Evaluate(*si).Clamp();
	si->bsdf = new BSDF(*si, eta);
	if (!R.IsBlack() && !T.IsBlack())
		si->bsdf->Add(new SpecularTransmission(T,1.f,eta,mode));	
}

void PlasticMaterial::ComputeScatteringFunctions(
	SurfaceInteraction *si,  TransportMode mode,
	bool allowMultipleLobes) const {
		si->bsdf = new BSDF(*si);
		Spectrum kd = Kd->Evaluate(*si).Clamp();
		if (!kd.IsBlack())
			si->bsdf->Add(new LambertianReflection(kd));
		Spectrum ks = Ks->Evaluate(*si).Clamp();
		if (!ks.IsBlack()) {
			Fresnel *fresnel = new FresnelDielectric(1.5f, 1.f);
			Float rough = roughness->Evaluate(*si);
			if (remapRoughness)
				rough = TrowbridgeReitzDistribution::RoughnessToAlpha(rough);
			MicrofacetDistribution *distrib =new TrowbridgeReitzDistribution(rough, rough);
			BxDF *spec =new MicrofacetReflection(ks, distrib, fresnel);
			si->bsdf->Add(spec);
		}
}

void MetalMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,
	bool allowMultipleLobes) const {
			si->bsdf = new BSDF(*si);
	Float uRough =
		uRoughness ? uRoughness->Evaluate(*si) : roughness->Evaluate(*si);
	Float vRough =
		vRoughness ? vRoughness->Evaluate(*si) : roughness->Evaluate(*si);
	if (remapRoughness) {
		uRough = TrowbridgeReitzDistribution::RoughnessToAlpha(uRough);
		vRough = TrowbridgeReitzDistribution::RoughnessToAlpha(vRough);
	}
	Fresnel *frMf = new FresnelConductor(1., eta->Evaluate(*si),k->Evaluate(*si));
	MicrofacetDistribution *distrib =
		new TrowbridgeReitzDistribution(uRough, vRough);
	si->bsdf->Add(new MicrofacetReflection(1., distrib, frMf));
}

void MetalRoughnessMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
	TransportMode mode,bool allowMultipleLobes) const {
	si->bsdf = new BSDF(*si);
	Spectrum c = base_color->Evaluate(*si).Clamp();
	if (metalroughness == nullptr)
		si->bsdf->Add(new LambertianReflection(c));
	else
	{
		Spectrum tmp = metalroughness->Evaluate(*si).Clamp();
		Float metallic = tmp[2], roughness = std::max(0.05f,tmp[1]);
		MicrofacetDistribution *distrib =new TrowbridgeReitzDistribution(roughness, roughness);
		si->bsdf->Add(new MetalRoughnessReflection(c, metallic, roughness, distrib));
	}
}