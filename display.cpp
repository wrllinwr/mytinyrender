#include "display.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_render.h>

static struct SDL_Window* s_window = NULL;
static struct SDL_Renderer* s_renderer = NULL;
static uint32_t* s_color_buffer = NULL;
static float* s_depth_buffer = NULL;
static struct SDL_Texture* s_color_buffer_texture = NULL;
// default/fallback window width/height
static int s_window_width = 800;
static int s_window_height = 800;

int32_t Display::fps(void) {
	return 60;
}

float Display::seconds_per_frame(void) {
	return 1.0f / (float)fps();
}

double Display::seconds_elapsed(
		const uint64_t previous_counter, const uint64_t current_counter) {
	return (double)(current_counter - previous_counter)
		/ (double)SDL_GetPerformanceFrequency();
}

bool Display::initialize_window(void) {
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
		fprintf(stderr, "Error initializing SDL.\n");
		return false;
	}

	SDL_DisplayMode display_mode;
	SDL_GetCurrentDisplayMode(0, &display_mode);
	const int window_width = display_mode.w;
	const int window_height = display_mode.h;

	const int scale_factor = 1; // increase to create pixelated effect
	s_window_width = window_width / scale_factor;
	s_window_height = window_height / scale_factor;
	s_window_width = 800;
	s_window_height = 800;

	s_window = SDL_CreateWindow(
			NULL,
			SDL_WINDOWPOS_CENTERED,
			SDL_WINDOWPOS_CENTERED,
			s_window_width,
			s_window_height,
			// 1280,
			// 720,
			SDL_WINDOW_BORDERLESS);

	if (!s_window) {
		fprintf(stderr, "Error creating SDL window.\n");
		return false;
	}

	s_renderer = SDL_CreateRenderer(s_window, -1, 0);

	if (!s_renderer) {
		fprintf(stderr, "Error creating SDL renderer.\n");
		return false;
	}

	// Apply a transformation to flip the y-axis
	SDL_RenderSetScale(s_renderer, 1.0f, -1.0f);
	return true;
}

void Display::draw_pixel(int x, int y, uint32_t color)
{
	if (x < 0 || x >= s_window_width || y <= 0
			|| y >= s_window_height) {
		return;
	}
	// printf("draw_pixel\n");
	// s_color_buffer[y * s_window_width + x] = color;
	s_color_buffer[(s_window_height - y - 1) * s_window_width + x] = color; // flip y-axis
}

void Display::draw_line(int x0, int y0, int x1, int y1, uint32_t color)
{
	// t1 1 website
	// for (float t=0.; t<1.; t+=.01) { 
	//     int x = x0 + (x1-x0)*t; 
	//     int y = y0 + (y1-y0)*t; 
	//     draw_pixel(x, y, color); 
	// }

	// t1 source code
	// for (float t=0.; t<1.; t+=.1) {
	//   int x = x0*(1.-t) + x1*t;
	//   int y = y0*(1.-t) + y1*t;
	//   draw_pixel(x, y, color);
	// }

	// t1 2 website
	// for (int x=x0; x<=x1; x++) { 
	//     float t = (x-x0)/(float)(x1-x0); 
	//     int y = y0*(1.-t) + y1*t; 
	//     draw_pixel(x, y, color); 
	// }

	// t1 3 source code
	bool steep = false;
	if (std::abs(x0-x1)<std::abs(y0-y1)) {
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0>x1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}

	for (int x=x0; x<=x1; x++) {
		float t = (x-x0)/(float)(x1-x0);
		int y = y0*(1.-t) + y1*t;
		if (steep) {
			draw_pixel(y, x, color);
		} else {
			draw_pixel(x, y, color);
		}
	}
}

void Display::draw_vec2_line(Vec2i p0, Vec2i p1, uint32_t color){
	bool steep = false;
	if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x>p1.x) {
		std::swap(p0, p1);
	}

	for (int x=p0.x; x<=p1.x; x++) {
		float t = (x-p0.x)/(float)(p1.x-p0.x);
		int y = p0.y*(1.-t) + p1.y*t;
		if (steep) {
			draw_pixel(y, x, color);
		} else {
			draw_pixel(x, y, color);
		}
	}
}

void Display::draw_vec3_line(Vec3i p0, Vec3i p1, uint32_t color) {
	bool steep = false;
	if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x>p1.x) {
		std::swap(p0, p1);
	}

	for (int x=p0.x; x<=p1.x; x++) {
		float t = (x-p0.x)/(float)(p1.x-p0.x);
		int y = p0.y*(1.-t) + p1.y*t+.5;
		if (steep) {
			draw_pixel(y, x, color);
		} else {
			draw_pixel(x, y, color);
		}
	}
}

void Display::triangle(Vec2i t0, Vec2i t1, Vec2i t2, uint32_t color){
	// t21 source code
	// draw_vect_line(t0, t1, color);
	// draw_vect_line(t1, t2, color);
	// draw_vect_line(t2, t0, color);

	// t22 website
	// A good method of drawing a triangle: It should be symmetrical, the picture should not depend on the order of vertices passed to the drawing function.
	// if (t0.y>t1.y) std::swap(t0, t1); 
	// if (t0.y>t2.y) std::swap(t0, t2); 
	// if (t1.y>t2.y) std::swap(t1, t2); 
	// draw_vect_line(t1, t2, 0xff00ff00);
	// draw_vect_line(t0, t1, 0xff00ff00);
	// draw_vect_line(t2, t0, 0xff0000ff);

	// t23 website Rasterization tirangle
	// sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
	// if (t0.y>t1.y) std::swap(t0, t1);
	// if (t0.y>t2.y) std::swap(t0, t2);
	// if (t1.y>t2.y) std::swap(t1, t2);
	// int total_height = t2.y-t0.y;
	// for (int y=t0.y; y<=t1.y; y++) {
	//     int segment_height = t1.y-t0.y+1;
	//     float alpha = (float)(y-t0.y)/total_height;
	//     float beta  = (float)(y-t0.y)/segment_height; // be careful with divisions by zero
	//     Vec2i A = t0 + (t2-t0)*alpha;
	//     Vec2i B = t0 + (t1-t0)*beta;
	//     // draw_pixel(A.x, y, 0xff0000ff); // Draw left line
	//     // draw_pixel(B.x, y, 0xff00ff00); // Draw right line
	//     draw_vect_line(A, B, 0xffffffff);
	// }

	// t24 website
	// sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
	// if (t0.y>t1.y) std::swap(t0, t1);
	// if (t0.y>t2.y) std::swap(t0, t2);
	// if (t1.y>t2.y) std::swap(t1, t2);
	// int total_height = t2.y-t0.y;
	// for (int y=t0.y; y<=t1.y; y++) {
	//     int segment_height = t1.y-t0.y+1;
	//     float alpha = (float)(y-t0.y)/total_height;
	//     float beta  = (float)(y-t0.y)/segment_height; // be careful with divisions by zero
	//     Vec2i A = t0 + (t2-t0)*alpha;
	//     Vec2i B = t0 + (t1-t0)*beta;
	//     if (A.x>B.x) std::swap(A, B);
	//     for (int j=A.x; j<=B.x; j++) {
	//   	  draw_pixel(j, y, color); // attention, due to int casts t0.y+i != A.y
	//     }
	// }
	// for (int y=t1.y; y<=t2.y; y++) {
	//     int segment_height =  t2.y-t1.y+1;
	//     float alpha = (float)(y-t0.y)/total_height;
	//     float beta  = (float)(y-t1.y)/segment_height; // be careful with divisions by zero
	//     Vec2i A = t0 + (t2-t0)*alpha;
	//     Vec2i B = t1 + (t2-t1)*beta;
	//     if (A.x>B.x) std::swap(A, B);
	//     for (int j=A.x; j<=B.x; j++) {
	//   	  draw_pixel(j, y, color); // attention, due to int casts t0.y+i != A.y
	//     }
	// }

	// t25 website
	if (t0.y==t1.y && t0.y==t2.y) return; // I dont care about degenerate triangles
										  // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
	if (t0.y>t1.y) std::swap(t0, t1);
	if (t0.y>t2.y) std::swap(t0, t2);
	if (t1.y>t2.y) std::swap(t1, t2);
	int total_height = t2.y-t0.y;
	for (int i=0; i<total_height; i++) {
		bool second_half = i>t1.y-t0.y || t1.y==t0.y;
		int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
		float alpha = (float)i/total_height;
		float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
		Vec2i A =               t0 + (t2-t0)*alpha;
		Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta;
		if (A.x>B.x) std::swap(A, B);
		for (int j=A.x; j<=B.x; j++) {
			draw_pixel(j, t0.y+i, color); // attention, due to int casts t0.y+i != A.y
		}
	}

	// t2x source code fill triangle
	// if (t0.y==t1.y && t0.y==t2.y) return; // i dont care about degenerate triangles
	// if (t0.y>t1.y) std::swap(t0, t1);
	// if (t0.y>t2.y) std::swap(t0, t2);
	// if (t1.y>t2.y) std::swap(t1, t2);
	// int total_height = t2.y-t0.y;
	// for (int i=0; i<total_height; i++) {
	//     bool second_half = i>t1.y-t0.y || t1.y==t0.y;
	//     int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
	//     float alpha = (float)i/total_height;
	//     float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
	//     Vec2i A =               t0 + (t2-t0)*alpha;
	//     Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta;
	//     if (A.x>B.x) std::swap(A, B);
	//     for (int j=A.x; j<=B.x; j++) {
	//         draw_pixel(j, t0.y+i, color); // attention, due to int casts t0.y+i != A.y
	//     }
	// }
}

void Display::clear_color_buffer(const uint32_t color) {
	for (int col = 0; col < s_window_width; ++col) {
		for (int row = 0; row < s_window_height; ++row) {
			s_color_buffer[row * s_window_width + col] = color;
		}
	}
}

void Display::clear_depth_buffer(void) {
	for (int col = 0; col < s_window_width; ++col) {
		for (int row = 0; row < s_window_height; ++row) {
			s_depth_buffer[row * s_window_width + col] = 1.0f;
		}
	}
}

void Display::render_color_buffer(void) {
	SDL_UpdateTexture(
			s_color_buffer_texture,
			NULL,
			s_color_buffer,
			s_window_width * sizeof(uint32_t));
	SDL_RenderCopy(s_renderer, s_color_buffer_texture, NULL, NULL);
}

void Display::deinitialize_window(void) {
	SDL_DestroyRenderer(s_renderer);
	SDL_DestroyWindow(s_window);
	SDL_Quit();
}

void Display::create_color_buffer(void) {
	s_color_buffer = (uint32_t*)malloc(sizeof(uint32_t) * s_window_width * s_window_height);
	s_color_buffer_texture = SDL_CreateTexture(
			s_renderer,
			SDL_PIXELFORMAT_RGBA32,
			SDL_TEXTUREACCESS_STREAMING,
			s_window_width,
			s_window_height);
}

void Display::destroy_color_buffer(void) {
	SDL_DestroyTexture(s_color_buffer_texture);
	free(s_color_buffer);
}

void Display::create_depth_buffer(void) {
	s_depth_buffer = (float*)malloc(sizeof(float) * s_window_width * s_window_height);
}

void Display::destroy_depth_buffer(void) {
	free(s_depth_buffer);
}

void Display::renderer_present(void) {
	SDL_RenderPresent(s_renderer);
}

int Display::window_width(void) {
	return s_window_width;
}

int Display::window_height(void) {
	return s_window_height;
}

// Vec3f Display::barycentric(Vec2i *pts, Vec2i P) {
Vec3f Display::barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
	Vec3f s[2];
	for (int i=2; i--; ) {
		s[i][0] = C[i]-A[i]; 
		s[i][1] = B[i]-A[i]; 
		s[i][2] = A[i]-P[i]; 
	}
	// printf("%f %f %f %f %f %f\n", s[0][0], s[0][1], s[0][2],s[1][0], s[1][1], s[1][2]);
	// Vec3f u = cross(s[0], s[1]);

	// Vec3f sy(C.y-A.y, B.y-A.y, A.y-P.y);
	// Vec3f sx(C.x-A.x, B.x-A.x, A.x-P.x);
	// Vec3f u = cross(sy, sx);

	// s[0].x = C.x-A.x;
	// s[0].y = B.x-A.x;
	// s[0].z = A.x-P.x;
	// s[1].x = C.y-A.y;
	// s[1].y = B.y-A.y;
	// s[1].z = A.y-P.y;
	// Vec3f u = cross(s[0], s[1]);
	Vec3f u = s[0] ^ s[1];
	if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
		return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
	return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void Display::triangle_barycentric(Vec3f *pts, float *zbuffer, uint32_t color){
	// t2
	// void Display::triangle_barycentric(Vec2i *pts, uint32_t color) {
	// Vec2i bboxmin(200,  200);
	// Vec2i bboxmax(0, 0);
	// Vec2i clamp(200, 200);
	// for (int i=0; i<3; i++) {
	//     bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
	// bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));

	// bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
	// bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
	// }
	// Vec2i P;
	// for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
	//     for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
	//         // Vec3f bc_screen  = barycentric_zh(pts[0], pts[1], pts[2], P);
	//         Vec3f bc_screen  = barycentric(pts, P);
	//         if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
	//         draw_pixel(P.x, P.y, color);
	//     }
	// }

	// t3
	Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max()); // Vec2f bboxmin(∞,∞);
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()); // Vec2f bboxmax(-∞,-∞);
	Vec2f clamp(800, 800);
	for (int i=0; i<3; i++) { // why 3; create a vec2?
		for (int j=0; j<2; j++) { // why 2
			bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}
	Vec3f P;
	for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
		for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
			if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
			P.z = 0;
			for (int i=0; i<3; i++) P.z += pts[i][2]*bc_screen[i];
			if (zbuffer[int(P.x+P.y*window_width())]<P.z) {
				zbuffer[int(P.x+P.y*window_width())] = P.z;
				draw_pixel(P.x, P.y, color);
			}
		}
	}
}

void Display::triangle_texture_tga(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, float intensity, int *zbuffer, Model* model, Lighting* light){
	if (t0.y==t1.y && t0.y==t2.y) return; // i dont care about degenerate triangles
	if (t0.y>t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
	if (t0.y>t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
	if (t1.y>t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }

	int total_height = t2.y-t0.y;
	for (int i=0; i<total_height; i++) {
		bool second_half = i>t1.y-t0.y || t1.y==t0.y;
		int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
		float alpha = (float)i/total_height;
		float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
		Vec3i A   =               t0  + Vec3f(t2-t0  )*alpha;
		Vec3i B   = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
		Vec2i uvA =               uv0 +      (uv2-uv0)*alpha;
		Vec2i uvB = second_half ? uv1 +      (uv2-uv1)*beta : uv0 +      (uv1-uv0)*beta;
		if (A.x>B.x) { std::swap(A, B); std::swap(uvA, uvB); }
		for (int j=A.x; j<=B.x; j++) {
			float phi = B.x==A.x ? 1. : (float)(j-A.x)/(float)(B.x-A.x);
			Vec3i   P = Vec3f(A) + Vec3f(B-A)*phi;
			Vec2i uvP =     uvA +   (uvB-uvA)*phi;
			int idx = P.x+P.y*window_width();
			if (zbuffer[idx]<P.z) {
				zbuffer[idx] = P.z;
				TGAColor color = model->diffuse(uvP);
				uint32_t c = (
						(uint32_t)color.a << 24) | // Shift alpha to the top 8 bits
					((uint32_t)color.r) | // Shift red to the next 8 bits
					((uint32_t)color.g << 8)  | // Shift green to the next 8 bits
					((uint32_t)color.b << 16);       // Leave blue in the bottom 8 bits
				draw_pixel(P.x, P.y, light->apply_light_intensity(c, intensity));
			}
		}
	}
}

void Display::triangle_camera(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, int *zbuffer, Lighting* light) {
    if (t0.y==t1.y && t0.y==t2.y) return; // i dont care about degenerate triangles
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(ity0, ity1); }
    if (t0.y>t2.y) { std::swap(t0, t2); std::swap(ity0, ity2); }
    if (t1.y>t2.y) { std::swap(t1, t2); std::swap(ity1, ity2); }

    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
        Vec3i A    =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B    = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
        float ityA =               ity0 +   (ity2-ity0)*alpha;
        float ityB = second_half ? ity1 +   (ity2-ity1)*beta : ity0 +   (ity1-ity0)*beta;
        if (A.x>B.x) { std::swap(A, B); std::swap(ityA, ityB); }
        for (int j=A.x; j<=B.x; j++) {
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(B.x-A.x);
            Vec3i    P = Vec3f(A) +  Vec3f(B-A)*phi;
            float ityP =    ityA  + (ityB-ityA)*phi;
            int idx = P.x+P.y*window_width();
            if (P.x>=window_width()||P.y>=window_height()||P.x<0||P.y<0) continue;
            if (zbuffer[idx]<P.z) {
                zbuffer[idx] = P.z;
                draw_pixel(P.x, P.y, light->apply_light_intensity(0xffffffff, ityP));
                // draw_pixel(P.x, P.y, 0xffffffff);
            }
        }
    }
}
