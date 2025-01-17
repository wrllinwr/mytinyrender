#include "lighting.h"
#include <algorithm>
// #include <as-ops.h>


uint32_t Lighting::apply_light_intensity(const uint32_t color, const float t) {
  // const float clamped_t = as_clamp_float(t, 0.0f, 1.0f);
  const float clamped_t =std::clamp(t, 0.0f, 1.0f);
  const uint32_t a = color & 0xff000000;
  const uint32_t r = (color & 0x00ff0000) * clamped_t;
  const uint32_t g = (color & 0x0000ff00) * clamped_t;
  const uint32_t b = (color & 0x000000ff) * clamped_t;
  return a | (r & 0x00ff0000) | (g & 0x0000ff00) | (b & 0x000000ff);
}
