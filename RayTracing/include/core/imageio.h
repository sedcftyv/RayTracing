#pragma once
#ifndef PBRT_CORE_IMAGEIO_H
#define PBRT_CORE_IMAGEIO_H

#include "pbrt.h"
#include "geometry.h"

std::unique_ptr<RGBSpectrum[]> ReadImage(const std::string &name,
	Point2i *resolution);
//RGBSpectrum *ReadImageEXR(const std::string &name, int *width,
//	int *height, Bounds2i *dataWindow = nullptr,
//	Bounds2i *displayWindow = nullptr);
//
//void WriteImage(const std::string &name, const Float *rgb,
//	const Bounds2i &outputBounds, const Point2i &totalResolution);






#endif