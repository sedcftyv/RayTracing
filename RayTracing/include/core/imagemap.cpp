#include<core/imagemap.h>

template <typename Tmemory, typename Treturn>
std::map<TexInfo, std::unique_ptr<MIPMap<Tmemory>>>
ImageTexture<Tmemory, Treturn>::textures;

template <typename Tmemory, typename Treturn>
MIPMap<Tmemory> *ImageTexture<Tmemory, Treturn>::GetTexture(
	const std::string &filename, bool doTrilinear, Float maxAniso,
	ImageWrap wrap, Float scale, bool gamma) {
	TexInfo texInfo(filename, doTrilinear, maxAniso, wrap, scale, gamma);
	if (textures.find(texInfo) != textures.end())
		return textures[texInfo].get();
	Point2i resolution;
	std::unique_ptr<RGBSpectrum[]> texels = ReadImage(filename, &resolution);
	if (!texels) {
		printf("Creating a constant grey texture to replace \"%s\".",filename.c_str());
		resolution.x = resolution.y = 1;
		RGBSpectrum *rgb = new RGBSpectrum[1];
		*rgb = RGBSpectrum(0.5f);
		texels.reset(rgb);
	}
						
	MIPMap<Tmemory> *mipmap = nullptr;
	if (texels) {
		std::unique_ptr<Tmemory[]> convertedTexels(new Tmemory[resolution.x * resolution.y]);
		for (int i = 0; i < resolution.x * resolution.y; ++i)
			convertIn(texels[i], &convertedTexels[i], scale, gamma);
		mipmap = new MIPMap<Tmemory>(resolution, convertedTexels.get(),doTrilinear, maxAniso, wrap);
	}
	else {
		Tmemory oneVal = scale;
		mipmap = new MIPMap<Tmemory>(Point2i(1, 1), &oneVal);
	}
	textures[texInfo].reset(mipmap);
	return mipmap;
}


template <typename Tmemory, typename Treturn>
ImageTexture<Tmemory, Treturn>::ImageTexture(
	std::unique_ptr<TextureMapping2D> mapping, const std::string &filename,
	bool doTrilinear, Float maxAniso, ImageWrap wrapMode, Float scale,bool gamma)
	: mapping(std::move(mapping)) {
	mipmap = GetTexture(filename, doTrilinear, maxAniso, wrapMode, scale, gamma);
}


template class ImageTexture<Float, Float>;
template class ImageTexture<RGBSpectrum, Spectrum>;