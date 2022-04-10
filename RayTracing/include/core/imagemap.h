#pragma once

#ifndef PBRT_TEXTURES_IMAGEMAP_H
#define PBRT_TEXTURES_IMAGEMAP_H

#include "pbrt.h"
#include "core/texture.h"
#include "core/mipmap.h"
#include "core/imageio.h"

struct TexInfo {
	TexInfo(const std::string &f, bool dt, Float ma, ImageWrap wm, Float sc,
		bool gamma)
		: filename(f),
		doTrilinear(dt),
		maxAniso(ma),
		wrapMode(wm),
		scale(sc),
		gamma(gamma) {}
	std::string filename;
	bool doTrilinear;
	Float maxAniso;
	ImageWrap wrapMode;
	Float scale;
	bool gamma;
	bool operator<(const TexInfo &t2) const {
		if (filename != t2.filename) return filename < t2.filename;
		if (doTrilinear != t2.doTrilinear) return doTrilinear < t2.doTrilinear;
		if (maxAniso != t2.maxAniso) return maxAniso < t2.maxAniso;
		if (scale != t2.scale) return scale < t2.scale;
		if (gamma != t2.gamma) return !gamma;
		return wrapMode < t2.wrapMode;
	}
};



template <typename Tmemory, typename Treturn>
class ImageTexture : public Texture<Treturn> {
public:
	// ImageTexture Public Methods
	ImageTexture(std::unique_ptr<TextureMapping2D> mapping,const std::string &filename, bool doTrilinear, Float maxAniso,
		ImageWrap wrapMode, Float scale, bool gamma) : mapping(std::move(mapping))
	{
		mipmap = GetTexture(filename, doTrilinear, maxAniso, wrapMode, scale, gamma);
	}
	static void ClearCache() {
		textures.erase(textures.begin(), textures.end());
	}
	Treturn Evaluate(const SurfaceInteraction &si) const {
		Vector2f dstdx, dstdy;
		Point2f st = mapping->Map(si, &dstdx, &dstdy);
		Tmemory mem = mipmap->Lookup(st, dstdx, dstdy);
		Treturn ret;
		convertOut(mem, &ret);
		return ret;
	}
private:
	// ImageTexture Private Methods
	static MIPMap<Tmemory> *GetTexture(const std::string &filename,
		bool doTrilinear, Float maxAniso,
		ImageWrap wrap, Float scale, bool gamma)
	{
		// Return _MIPMap_ from texture cache if present
		TexInfo texInfo(filename, doTrilinear, maxAniso, wrap, scale, gamma);
		if (textures.find(texInfo) != textures.end())
			return textures[texInfo].get();

		// Create _MIPMap_ for _filename_
	//	ProfilePhase _(Prof::TextureLoading);
		Point2i resolution;
		std::unique_ptr<RGBSpectrum[]> texels = ReadImage(filename, &resolution);
		if (!texels) {
			printf("Creating a constant grey texture to replace \"%s\".",
				filename.c_str());
			resolution.x = resolution.y = 1;
			RGBSpectrum *rgb = new RGBSpectrum[1];
			*rgb = RGBSpectrum(0.5f);
			texels.reset(rgb);
		}

		// Flip image in y; texture coordinate space has (0,0) at the lower
		// left corner.
		for (int y = 0; y < resolution.y / 2; ++y)
			for (int x = 0; x < resolution.x; ++x) {
				int o1 = y * resolution.x + x;
				int o2 = (resolution.y - 1 - y) * resolution.x + x;
				std::swap(texels[o1], texels[o2]);
			}

		MIPMap<Tmemory> *mipmap = nullptr;
		if (texels) {
			// Convert texels to type _Tmemory_ and create _MIPMap_
			std::unique_ptr<Tmemory[]> convertedTexels(
				new Tmemory[resolution.x * resolution.y]);
			for (int i = 0; i < resolution.x * resolution.y; ++i)
				convertIn(texels[i], &convertedTexels[i], scale, gamma);
			mipmap = new MIPMap<Tmemory>(resolution, convertedTexels.get(),
				doTrilinear, maxAniso, wrap);
		}
		else {
			// Create one-valued _MIPMap_
			Tmemory oneVal = scale;
			mipmap = new MIPMap<Tmemory>(Point2i(1, 1), &oneVal);
		}
		textures[texInfo].reset(mipmap);
		return mipmap;
	}
	static void convertIn(const RGBSpectrum &from, RGBSpectrum *to, Float scale,
		bool gamma) {
		for (int i = 0; i < RGBSpectrum::nSamples; ++i)
			(*to)[i] = scale * (gamma ? InverseGammaCorrect(from[i]) : from[i]);
	}
	static void convertIn(const RGBSpectrum &from, Float *to, Float scale,
		bool gamma) {
		*to = scale * (gamma ? InverseGammaCorrect(from.y()) : from.y());
	}
	static void convertOut(const RGBSpectrum &from, Spectrum *to) {
		Float rgb[3];
		from.ToRGB(rgb);
		*to = Spectrum::FromRGB(rgb);
	}
	static void convertOut(Float from, Float *to) { *to = from; }

	// ImageTexture Private Data
	unique_ptr<TextureMapping2D> mapping;
	MIPMap<Tmemory> *mipmap;

	static std::map<TexInfo, std::unique_ptr<MIPMap<Tmemory> > > textures;
};

template <typename Tmemory, typename Treturn>
std::map<TexInfo, std::unique_ptr<MIPMap<Tmemory> > > ImageTexture<Tmemory, Treturn>::textures;


extern template class ImageTexture<Float, Float>;
extern template class ImageTexture<RGBSpectrum, Spectrum>;

//template <typename Tmemory, typename Treturn>
//MIPMap<Tmemory> *ImageTexture<Tmemory, Treturn>::GetTexture(
//	const std::string &filename, bool doTrilinear, Float maxAniso,
//	ImageWrap wrap, Float scale, bool gamma) {
//	// Return _MIPMap_ from texture cache if present
//	TexInfo texInfo(filename, doTrilinear, maxAniso, wrap, scale, gamma);
//	if (textures.find(texInfo) != textures.end())
//		return textures[texInfo].get();
//
//	// Create _MIPMap_ for _filename_
////	ProfilePhase _(Prof::TextureLoading);
//	Point2i resolution;
//	std::unique_ptr<RGBSpectrum[]> texels = ReadImage(filename, &resolution);
//	if (!texels) {
//		printf("Creating a constant grey texture to replace \"%s\".",
//			filename.c_str());
//		resolution.x = resolution.y = 1;
//		RGBSpectrum *rgb = new RGBSpectrum[1];
//		*rgb = RGBSpectrum(0.5f);
//		texels.reset(rgb);
//	}
//
//	// Flip image in y; texture coordinate space has (0,0) at the lower
//	// left corner.
//	for (int y = 0; y < resolution.y / 2; ++y)
//		for (int x = 0; x < resolution.x; ++x) {
//			int o1 = y * resolution.x + x;
//			int o2 = (resolution.y - 1 - y) * resolution.x + x;
//			std::swap(texels[o1], texels[o2]);
//		}
//
//	MIPMap<Tmemory> *mipmap = nullptr;
//	if (texels) {
//		// Convert texels to type _Tmemory_ and create _MIPMap_
//		std::unique_ptr<Tmemory[]> convertedTexels(
//			new Tmemory[resolution.x * resolution.y]);
//		for (int i = 0; i < resolution.x * resolution.y; ++i)
//			convertIn(texels[i], &convertedTexels[i], scale, gamma);
//		mipmap = new MIPMap<Tmemory>(resolution, convertedTexels.get(),
//			doTrilinear, maxAniso, wrap);
//	}
//	else {
//		// Create one-valued _MIPMap_
//		Tmemory oneVal = scale;
//		mipmap = new MIPMap<Tmemory>(Point2i(1, 1), &oneVal);
//	}
//	textures[texInfo].reset(mipmap);
//	return mipmap;
//}

//template <typename Tmemory, typename Treturn>
//ImageTexture<Tmemory, Treturn>::ImageTexture(
//	std::unique_ptr<TextureMapping2D> mapping, const std::string &filename,
//	bool doTrilinear, Float maxAniso, ImageWrap wrapMode, Float scale,
//	bool gamma)
//	: mapping(std::move(mapping)) {
//	mipmap = GetTexture(filename, doTrilinear, maxAniso, wrapMode, scale, gamma);
//}



#endif