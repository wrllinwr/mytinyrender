#include "display.h"
#include "fps.h"
#include <cmath>
#include <limits>
#include <vector>

#define MODEL_SIZE 1

Display* display = nullptr;
uint64_t g_previous_frame_time = 0;
Fps g_fps = {.head_ = 0, .tail_ = FpsMaxSamples - 1};
// Model* model = nullptr;
Model* model[MODEL_SIZE];
Lighting* light = nullptr;
Matrix* matrix;
Vec3f camera(0,0,3);
// int *zbuffer = NULL;
int*  zbuffer = new int[display->window_width()*display->window_height()];
// Vec3f light_dir(0,0,-1);
Vec3f light_dir = Vec3f(1,-1,1).normalize();
Vec3f eye(1,1,3);
Vec3f center(0,0,0);

void setup(void){
	display->create_color_buffer();
	display->create_depth_buffer();

	model[0] = new Model("obj/african_head/african_head.obj");
	// model[1] = new Model("obj/african_head/african_head_eye_inner.obj");
	// model[2] = new Model("obj/african_head/african_head_eye_outer.obj");
	// model = new Model("obj/boggie/body.obj");
	// model = new Model("obj/boggie/head.obj");
	// model = new Model("obj/boggie/eyes.obj");
	// model = new Model("obj/diablo3_pose/diablo3_pose.obj");
	// model[0] = new Model("obj/girl/girl.obj");
	// model[0] = new Model("obj/cube.obj");
	// model[1] = new Model("obj/nude/female_nude.obj"); // f=4
}

bool process_input(void) {
	for (SDL_Event event; SDL_PollEvent(&event) != 0;) {
		switch (event.type) {
			case SDL_QUIT:
				return false;
			case SDL_MOUSEMOTION: 
				break;
			case SDL_MOUSEBUTTONDOWN: 
				break;
			case SDL_MOUSEBUTTONUP: 
				break;
			case SDL_KEYDOWN: 
				if (event.key.keysym.sym == SDLK_ESCAPE) {
					return false;
				}
				break;
			case SDL_KEYUP: 
				break;
			default:
				break;
		}
	}
	return true;
}

void wait_to_update(void) {
	// reference: https://davidgow.net/handmadepenguin/ch18.html
	// seconds elapsed since last update
	const double seconds =
		display->seconds_elapsed(g_previous_frame_time, SDL_GetPerformanceCounter());
	if (seconds < display->seconds_per_frame()) {
		const double remainder_s = (double)display->seconds_per_frame() - seconds;
		// wait 4ms less than actual remainder due to resolution of SDL_Delay
		// (we don't want to delay/sleep too long and get behind)
		const double remainder_pad_s = remainder_s - 0.004;
		const double remainder_pad_ms = remainder_pad_s * 1000.0;
		const double remainder_pad_ms_clamped = fmax(remainder_pad_ms, 0.0);
		const uint32_t delay = (uint32_t)remainder_pad_ms_clamped;
		SDL_Delay(delay);
		// busy wait for the remaining time
		while (display->seconds_elapsed(g_previous_frame_time, SDL_GetPerformanceCounter())
				< display->seconds_per_frame()) {
			;
		}
	}
}

void calculate_framerate(void) {
	// const int64_t time_window =
	//   calculate_window(&g_fps, SDL_GetPerformanceCounter());
	// if (time_window != -1) {
	//   const double framerate =
	//     (double)(FpsMaxSamples - 1)
	//     / ((double)(time_window) / SDL_GetPerformanceFrequency());
	//   // fprintf(stderr, "fps: %f\n", framerate);
	// }
}

void update(void) {
	wait_to_update();

	const int64_t current_counter = SDL_GetPerformanceCounter();
	// const double delta_time =
	//   display->seconds_elapsed(g_previous_frame_time, current_counter);
	g_previous_frame_time = current_counter;

	// calculate_framerate();

	// update_movement(delta_time);
	// printf("update\n");
}

Vec3f world2screen(Vec3f v) {
	return Vec3f(int((v.x+1.)*display->window_width()/2.+.5), int((v.y+1.)*display->window_height()/2.+.5), v.z); // head.model
	// return Vec3f(int((v.x+1.)*display->window_width()/2.+.5), int((v.y+1.)*display->window_height()/2.+.5 -(display->window_height()/2.) + 20), v.z); // girl.model
}

void render(void) {
	display->clear_color_buffer(0xff000000); // alpa blue green red
	display->clear_depth_buffer();
	// t5
	for (int i=0; i<display->window_width()*display->window_height(); i++) {
		zbuffer[i] = std::numeric_limits<int>::min();
	}

	{ // draw the model
		Matrix ModelView  = matrix->lookat(eye, center, Vec3f(0,1,0));
		Matrix Projection = Matrix::identity(4);
		Matrix ViewPort   = matrix->viewport(display->window_width()/8, display->window_height()/8, display->window_width()*3/4, display->window_height()*3/4);
		Projection[3][2] = -1.f/(eye - center).norm();

		std::cerr << ModelView << std::endl;
		std::cerr << Projection << std::endl;
		std::cerr << ViewPort << std::endl;
		Matrix z = (ViewPort*Projection*ModelView);
		std::cerr << z << std::endl;

		for (int m = 0; m < MODEL_SIZE; m++){
			for (int i=0; i<model[m]->nfaces(); i++) {
				std::vector<int> face = model[m]->face(i);
				Vec3i screen_coords[3];
				Vec3f world_coords[3];
				float intensity[3];
				for (int j=0; j<3; j++) {
					Vec3f v = model[m]->vert(face[j]);
					screen_coords[j] =  Vec3f(matrix->m2v(ViewPort*Projection*ModelView*matrix->v2m(v)));
					world_coords[j]  = v;
					intensity[j] = model[m]->norm(i, j)*light_dir;
				}
				display->triangle_camera(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0], intensity[1], intensity[2], zbuffer, light);
			}
		}
	}
	
	{ // dump z-buffer (debugging purposes only)
		TGAImage zbimage(800, 800, TGAImage::GRAYSCALE);
		for (int i=0; i<800; i++) {
			for (int j=0; j<800; j++) {
				zbimage.set(i, j, TGAColor(zbuffer[i+j*800], 1));
			}
		}
		zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		zbimage.write_tga_file("zbuffer.tga");
	}

	// Step 3
	display->render_color_buffer();
	display->renderer_present();
}

void teardown(void) {
	for(int i = 0; i < sizeof(model)/sizeof(model[0]); i++){
		printf("teardwon delete model[%d]\n", i);
		delete model[i];
	}
	delete display;
	delete light;
	delete matrix;
	delete [] zbuffer;
    display->destroy_depth_buffer();
	display->destroy_color_buffer();
	display->deinitialize_window();
}

int main(int argc, char** argv) {
	bool is_running = display->initialize_window();

	setup();

	g_previous_frame_time = SDL_GetPerformanceCounter();
	while (is_running) {
		is_running = process_input();
		update();
		render();
	}

	teardown();

	return 0;
}
