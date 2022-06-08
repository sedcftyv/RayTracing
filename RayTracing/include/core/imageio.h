#pragma once
#ifndef PBRT_CORE_IMAGEIO_H
#define PBRT_CORE_IMAGEIO_H

#include "pbrt.h"
#include "geometry.h"

std::unique_ptr<RGBSpectrum[]> ReadImage(const std::string &name,
	Point2i *resolution);






#endif