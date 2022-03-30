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