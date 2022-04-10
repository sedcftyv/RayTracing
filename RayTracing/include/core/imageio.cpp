#include<core/imageio.h>
#include<core/spectrum.h>
#include<ext/lodepng.h>

static RGBSpectrum *ReadImagePNG(const std::string &name, int *w, int *h);

std::unique_ptr<RGBSpectrum[]> ReadImage(const std::string &name,
	Point2i *resolution) {
	if (HasExtension(name, ".png"))
		return std::unique_ptr<RGBSpectrum[]>(
			ReadImagePNG(name, &resolution->x, &resolution->y));
	printf("Unable to load image stored in format \"%s\" for filename \"%s\".",
		strrchr(name.c_str(), '.') ? (strrchr(name.c_str(), '.') + 1)
		: "(unknown)",
		name.c_str());
	return nullptr;
}

static RGBSpectrum *ReadImagePNG(const std::string &name, int *width,
	int *height) {
	unsigned char *rgb;
	unsigned w, h;
	unsigned int error = lodepng_decode24_file(&rgb, &w, &h, name.c_str());
	if (error != 0) {
		printf("Error reading PNG\n");
		return nullptr;
	}
	*width = w;
	*height = h;

	RGBSpectrum *ret = new RGBSpectrum[*width * *height];
	unsigned char *src = rgb;
	for (unsigned int y = 0; y < h; ++y) {
		for (unsigned int x = 0; x < w; ++x, src += 3) {
			Float c[3];
			c[0] = src[0] / 255.f;
			c[1] = src[1] / 255.f;
			c[2] = src[2] / 255.f;
			ret[y * *width + x] = RGBSpectrum::FromRGB(c);
		}
	}

	free(rgb);
	//LOG(INFO) << StringPrintf("Read PNG image %s (%d x %d)",
	//	name.c_str(), *width, *height);
	return ret;
}