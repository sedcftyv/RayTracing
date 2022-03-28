#pragma once
#ifndef PDF_H
#define PDF_H

//#include "Vec3.h"
#include "core/geometry.h"
#include "rtweekend.h"
#include "onb.h"

class pdf {
public:
	virtual ~pdf() {}

	virtual double value(const Vector3f& direction) const = 0;
	virtual Vector3f generate() const = 0;
};

class cosine_pdf : public pdf {
public:
	cosine_pdf(const Vector3f& w) { uvw.build_from_w(w); }

	virtual double value(const Vector3f& direction) const override {
		auto cosine = Dot(Normalize(direction), uvw.w());
		return (cosine <= 0) ? 0 : cosine / pi;
	}

	virtual Vector3f generate() const override {
		return uvw.local(random_cosine_direction());
	}

public:
	onb uvw;
};

class hittable_pdf : public pdf {
public:
	hittable_pdf(shared_ptr<hittable> p, const Vector3f& origin) : ptr(p), o(origin) {}

	virtual double value(const Vector3f& direction) const override {
		return ptr->pdf_value(o, direction);
	}

	virtual Vector3f generate() const override {
		return ptr->random(o);
	}

public:
	Vector3f o;
	shared_ptr<hittable> ptr;
};

class mixture_pdf : public pdf {
public:
	mixture_pdf(shared_ptr<pdf> p0, shared_ptr<pdf> p1) {
		p[0] = p0;
		p[1] = p1;
	}

	virtual double value(const Vector3f& direction) const override {
		return 0.5 * p[0]->value(direction) + 0.5 *p[1]->value(direction);
	}

	virtual Vector3f generate() const override {
		if (random_double() < 0.5)
			return p[0]->generate();
		else
			return p[1]->generate();
	}

public:
	shared_ptr<pdf> p[2];
};


#endif