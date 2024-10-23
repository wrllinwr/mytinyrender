#ifndef LIGHTING_H
#define LIGHTING_H

#include <stdint.h>

class Lighting{
	public:
		uint32_t apply_light_intensity(uint32_t color, float t);
};

#endif // LIGHTING_H
