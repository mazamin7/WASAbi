/* WASAbi 2.5D
 * 
 * This is the entrance of the program. 
 * SDL is used as the interface to show the wave propagation.
 */

#include <iostream>
#include <SDL.h>
#include <SDL_ttf.h>
#include <omp.h>
#include <Windows.h>
#undef main		// https://stackoverflow.com/questions/6847360

#include "simulation.h"
#include "partition.h"
#include "boundary.h"
#include "sound_source.h"
#include "gaussian_source.h"
#include "recorder.h"

using namespace std;

bool is_record_response = false;
bool is_record_field = false;

/* Set constant parameters. */

double Partition::boundary_absorption_ = 0.5;	// Absorption coefficients of the boundaries.
double Simulation::air_absorption_alpha1_ = 0.0; // Coefficient of constant part of air absorption.
double Simulation::air_absorption_alpha2_ = 1e-6; // Coefficient of frequency dependent part of air absorption.
double Simulation::duration_ = 5e-1;		// Duration of the whole simulation (seconds).

//double Simulation::dh_ = 0.05;			// Space sampling rate.
//double Simulation::dt_ = 0.625e-4;		// Time sampling rate.

//double Simulation::dh_ = 0.1;
//double Simulation::dt_ = 1.25e-4;

//double Simulation::dh_ = 0.2;
//double Simulation::dt_ = 2e-4;

double Simulation::dh_ = 0.5;
double Simulation::dt_ = 6.25e-4;

double Simulation::c0_ = 3.435e2;		// Speed of sound
int Simulation::n_pml_layers_ = 5;		// Number of pml layers.

int main()
{
	double time1 = omp_get_wtime();		// Record the begining time. Used for showing the consuming time.

	// Check if nested parallelism is enabled
	if (omp_get_nested()) {
		std::cout << "Nested parallelism is already enabled." << std::endl;
	}
	else {
		std::cout << "Nested parallelism is NOT enabled." << std::endl;
		std::cout << "Enabling nested parallelism..." << std::endl;
		omp_set_nested(1);  // Enable nested parallelism
	}

	// Verify nested parallelism
	if (omp_get_nested()) {
		std::cout << "Nested parallelism has been successfully enabled." << std::endl;
	}
	else {
		std::cout << "Failed to enable nested parallelism." << std::endl;
	}

	std::string dir_name = "./output/" + std::to_string(Simulation::dh_) + "_" + std::to_string(Partition::boundary_absorption_) + "_" + std::to_string(Simulation::air_absorption_alpha1_) + "_" + std::to_string(Simulation::air_absorption_alpha2_);
	CreateDirectory(dir_name.c_str(), NULL);	// Prepare for the output folder.
												// ! Without this and the corresponding folder does not exist, the program will not write the output data.

	std::vector<std::shared_ptr<Partition>> partitions;
	std::vector<std::shared_ptr<SoundSource>> sources;
	std::vector<std::shared_ptr<Recorder>> recorders;

	partitions = Partition::ImportPartitions("./assets/hall_parallel_test.txt");			// Read partition properties from file.
	sources = SoundSource::ImportSources("./assets/hall_parallel_test-sources.txt");		// Read source properties from file.
	recorders = Recorder::ImportRecorders("./assets/hall_parallel_test-recorders.txt");	// Read recorder properties from file. Recorder is not mandatory. 

	for (auto record : recorders)
	{
		record->FindPartition(partitions);		// Assign recorders to the corresponding partition.
	}

	auto simulation = std::make_shared<Simulation>(partitions, sources);	// Initialize the simulation.
	simulation->Info();														// Show basic info of the simulation
	simulation->look_from_ = 0;											// FOR DEBUG: show field from another view direction.

	/* Initialize SDL window
	 * simulation_rect: show field.
	 * message_rect: show instant progress during the simulation
	 */
	SDL_Event event;
	SDL_Init(SDL_INIT_VIDEO);
	SDL_PixelFormat* fmt = SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888);
	int resolution_x = 800;
	int resolution_y = resolution_x / simulation->size_x()*simulation->size_y();
	SDL_Window* window = SDL_CreateWindow("WASAbi 2.5D Simulator",
		SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, resolution_x, resolution_y + 20, 0);
	SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
	SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB888, SDL_TEXTUREACCESS_STREAMING,
		simulation->size_x(), simulation->size_y());

	SDL_Rect simulation_rect;
	simulation_rect.x = 0;
	simulation_rect.y = 0;
	simulation_rect.w = resolution_x;
	simulation_rect.h = resolution_y;

	TTF_Init();
	TTF_Font* Sans = TTF_OpenFont("font/SourceSansPro-Regular.ttf", 64); //this opens a font style and sets a size
	SDL_Color White = { 255, 255, 255 };  // this is the color in rgb format, maxing out all would give you the color white, and it will be your text's color
	SDL_Rect Message_rect;
	Message_rect.x = 0;
	Message_rect.y = resolution_y;
	Message_rect.w = 100;
	Message_rect.h = 20;
	SDL_Rect Message_rect2;
	Message_rect2.x = resolution_x - 100;
	Message_rect2.y = resolution_y;
	Message_rect2.w = 100;
	Message_rect2.h = 20;

	bool quit = false;
	int time_step = 0;
	int total_time_steps = Simulation::duration_ / Simulation::dt_;
	std::string message;

	double time2 = omp_get_wtime();
	/*std::cout << omp_get_num_procs() << std::endl;
	omp_set_num_threads(omp_get_num_procs())*/;
	std::cout << "Initialization finished. (" << time2 - time1 << " s)" << std::endl;
	std::cout << "############################################################" << std::endl;

	while (!quit && time_step < total_time_steps)
	{
		while (SDL_PollEvent(&event)) {
			switch (event.type)
			{
			case SDL_QUIT:
				quit = true;
				break;
			}
		}

		time_step = simulation->Update();		// ! Updating sound field.

		if (is_record_response)
		{
			for (auto record : recorders)
			{
				record->RecordResponse(time_step);	// Record room's response.
			}
		}

		if (is_record_field)
		{
			for (auto record : recorders)
			{
				record->RecordField(time_step);	// Record sound field.
			}
		}

		message = std::to_string(time_step) + '/' + std::to_string(total_time_steps);
		SDL_Surface* surfaceMessage = TTF_RenderText_Solid(Sans, message.c_str(), White);
		SDL_Texture* Message = SDL_CreateTextureFromSurface(renderer, surfaceMessage);

		if (simulation->ready())
		{
			SDL_UpdateTexture(texture, nullptr,
				simulation->pixels().data(), simulation->size_x() * sizeof(Uint32));
		}
		SDL_RenderClear(renderer);
		SDL_RenderCopy(renderer, texture, nullptr, &simulation_rect);
		SDL_RenderCopy(renderer, Message, NULL, &Message_rect);

		message = std::to_string(static_cast<int>(floor((omp_get_wtime() - time1) / 60))) + " min, " + std::to_string(static_cast<int>(floor((omp_get_wtime() - time1))) % 60) + " sec";
		surfaceMessage = TTF_RenderText_Solid(Sans, message.c_str(), White);
		Message = SDL_CreateTextureFromSurface(renderer, surfaceMessage);
		SDL_RenderCopy(renderer, Message, NULL, &Message_rect2);

		SDL_FreeSurface(surfaceMessage);
		SDL_DestroyTexture(Message);
		SDL_RenderPresent(renderer);
	}

	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	double time3 = omp_get_wtime();
	std::cout << std::endl << "Simulation finished. (" << time3 - time1 << " s)" << std::endl;
	std::cout << "############################################################" << std::endl;

	return 0;
}

