#ifndef DISPLAY_H
#define DISPLAY_H

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "geometry.h"
#include "tgaimage.h"
#include "model.h"
#include "lighting.h"

class Display{
	public:
		int32_t fps(void);
		float seconds_per_frame(void);
		double seconds_elapsed(const uint64_t previous_counter, const uint64_t current_counter);
		bool initialize_window(void);
		void draw_pixel(int x, int y, uint32_t color);
		void draw_line(int x0, int y0, int x1, int y1, uint32_t color);
		void draw_vec2_line(Vec2i p0, Vec2i p1, uint32_t color);
		void draw_vec3_line(Vec3i p0, Vec3i p1, uint32_t color);
 		void triangle(Vec2i t0, Vec2i t1, Vec2i t2, uint32_t color);
		void clear_color_buffer(const uint32_t color);
		void clear_depth_buffer(void);
		void render_color_buffer(void);
		void deinitialize_window(void);
		void create_color_buffer(void);
		void destroy_color_buffer(void);
		void create_depth_buffer(void);
		void destroy_depth_buffer(void);
		void renderer_present(void);
		int window_width(void);
		int window_height(void);
		// Vec3f barycentric(Vec2i *pts, Vec2i P);
		Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P);
		// void triangle_barycentric(Vec2i *pts, uint32_t color);
		void triangle_barycentric(Vec3f *pts, float *zbuffer, uint32_t color);
		void triangle_texture_tga(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, float intensity, int *zbuffer, Model* model, Lighting* light);
		void triangle_camera(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, int *zbuffer, Lighting* light);
};



#endif
