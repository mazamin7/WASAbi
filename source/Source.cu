/* WASAbi 2.5D
 * 
 * This is the entrance of the program. 
 * SDL is used as the interface to show the wave propagation.
 */

#include <iostream>
#include <SDL.h>
#include <SDL_ttf.h>
#include <omp.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <iomanip>

#undef main
#include "ini.h"
#include "simulation.h"
#include "partition.h"
#include "boundary.h"
#include "sound_source.h"
#include "gaussian_source.h"
#include "recorder.h"

using namespace std;

bool is_record_response = false;
bool is_record_field = false;

double Partition::boundary_absorption_ = 0.5;
double Simulation::air_absorption_alpha1_ = 0.0;
double Simulation::air_absorption_alpha2_ = 1e-6;
double Simulation::duration_ = 2e-2;
double Simulation::c0_ = 343.5;
double Simulation::dh_ = 0.2;
double Simulation::dt_ = 2e-4;
int Simulation::n_pml_layers_ = 5;
int Simulation::viz_skip_ = 10;
float Simulation::max_viz_gain_ = 100.0f;

struct Config {
    string asset_name;
    double boundary_absorption;
    double air_absorption_alpha1;
    double air_absorption_alpha2;
    double duration;
    double c0;
    int n_pml_layers;
    string precision;
    int viz_skip;
    int fixed_panel_size;
    float max_viz_gain;
};

int parse_ini_handler(void* user, const char* section, const char* name, const char* value) {
    Config* config = (Config*)user;
    if (strcmp(section, "simulation") == 0) {
        if (strcmp(name, "asset_name") == 0) config->asset_name = value;
        else if (strcmp(name, "boundary_absorption") == 0) config->boundary_absorption = atof(value);
        else if (strcmp(name, "air_absorption_alpha1") == 0) config->air_absorption_alpha1 = atof(value);
        else if (strcmp(name, "air_absorption_alpha2") == 0) config->air_absorption_alpha2 = atof(value);
        else if (strcmp(name, "duration") == 0) config->duration = atof(value);
        else if (strcmp(name, "c0") == 0) config->c0 = atof(value);
        else if (strcmp(name, "n_pml_layers") == 0) config->n_pml_layers = atoi(value);
        else if (strcmp(name, "precision") == 0) config->precision = value;
        else if (strcmp(name, "viz_skip") == 0) config->viz_skip = atoi(value);
        else if (strcmp(name, "fixed_panel_size") == 0) config->fixed_panel_size = atoi(value);
        else if (strcmp(name, "max_viz_gain") == 0) config->max_viz_gain = atof(value);
    }
    return 1;
}

Config load_config(const string& filename) {
    Config config;
    config.viz_skip = 10;
    config.fixed_panel_size = 400;
    config.max_viz_gain = 100.0f;
    if (ini_parse(filename.c_str(), parse_ini_handler, &config) < 0) {
        cerr << "Can't load " << filename << endl;
    }
    return config;
}

void set_precision_params(const string& precision, double& dh, double& dt) {
    if (precision == "coarse") { dh = 0.5; dt = 6.25e-4; }
    else if (precision == "fine") { dh = 0.2; dt = 2e-4; }
    else if (precision == "finer") { dh = 0.1; dt = 1.25e-4; }
    else if (precision == "finest") { dh = 0.05; dt = 0.625e-4; }
    else { dh = 0.2; dt = 2e-4; }
}

void ensureConfigExists(const std::string& config_path, const std::string& default_path) {
    if (!std::filesystem::exists(config_path)) {
        std::filesystem::copy(default_path, config_path);
    }
}

enum class RunMode { UNKNOWN, SIM_RECORD_FIELD, SIM_RECORD_RESPONSE, SIM_VIZ, VIZ_RECORD };

struct CliArgs {
    RunMode mode = RunMode::UNKNOWN;
    string config_path = "./config/config.ini";
    string playback_file = "";
    int playback_delay = 33;
};

CliArgs parse_cli_args(int argc, char* argv[]) {
    CliArgs args;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) {
            string m = argv[++i];
            if (m == "sim-record-field") args.mode = RunMode::SIM_RECORD_FIELD;
            else if (m == "sim-record-response") args.mode = RunMode::SIM_RECORD_RESPONSE;
            else if (m == "sim-viz") args.mode = RunMode::SIM_VIZ;
            else if (m == "viz-record") args.mode = RunMode::VIZ_RECORD;
        }
        else if (arg == "--config" && i + 1 < argc) args.config_path = argv[++i];
        else if (arg == "--playback-file" && i + 1 < argc) args.playback_file = argv[++i];
        else if (arg == "--playback-delay" && i + 1 < argc) args.playback_delay = atoi(argv[++i]);
    }
    return args;
}

Uint32 calculate_color_playback(double p, float v_coef) {
    double norm = 0.5 * fmax(-1.0, fmin(1.0, p * v_coef)) + 0.5;
    int r, g, b;
    if (norm >= 0.5) {
        double pos = (norm - 0.5) * 2.0;
        r = 255; g = (int)(255 * (1.0 - pos)); b = (int)(255 * (1.0 - pos));
    } else {
        double neg = norm * 2.0;
        r = (int)(255 * neg); g = (int)(255 * neg); b = 255;
    }
    // Color normalization and shift to match SDL_PIXELFORMAT_RGB888 (original)
    // RGB888 expects 24 bits: R:16-23, G:8-15, B:0-7
    return (r << 16) | (g << 8) | b;
}

int main(int argc, char* argv[]) {
    // 1. Shift current directory to root if running from build
    try {
        std::filesystem::path cp = std::filesystem::current_path();
        if (cp.filename() == "build") {
            std::filesystem::current_path("..");
            // cout << "Running from build. Shifted working directory to project root: " << std::filesystem::current_path() << endl;
        }
    } catch (const std::exception& e) {
        cerr << "Directory error: " << e.what() << endl;
    }

    // 2. Parse arguments
    CliArgs cli_args = parse_cli_args(argc, argv);
    if (cli_args.mode == RunMode::UNKNOWN) {
        cout << "Usage: WASAbiApp --mode [sim-record-field|sim-record-response|sim-viz|viz-record] ...\n";
        return 0;
    }

    // 3. Ensure config and load it
    ensureConfigExists(cli_args.config_path, "./config/default.ini");
    Config config = load_config(cli_args.config_path);

    Partition::boundary_absorption_ = config.boundary_absorption;
    Simulation::air_absorption_alpha1_ = config.air_absorption_alpha1;
    Simulation::air_absorption_alpha2_ = config.air_absorption_alpha2;
    Simulation::duration_ = config.duration;
    Simulation::c0_ = config.c0;
    Simulation::n_pml_layers_ = config.n_pml_layers;
    set_precision_params(config.precision, Simulation::dh_, Simulation::dt_);
    Simulation::viz_skip_ = config.viz_skip;
    Simulation::max_viz_gain_ = config.max_viz_gain;

    if (cli_args.mode == RunMode::SIM_RECORD_FIELD) { is_record_field = true; is_record_response = false; }
    else if (cli_args.mode == RunMode::SIM_RECORD_RESPONSE) { is_record_field = false; is_record_response = true; }
    else { is_record_field = false; is_record_response = false; }

    double time1 = omp_get_wtime();
    cout << "Current Working Directory: " << std::filesystem::current_path() << endl;

    string dir_name = "./output/" + to_string(Simulation::dh_) + "_" + to_string(Partition::boundary_absorption_) + "_" + to_string(Simulation::air_absorption_alpha1_) + "_" + to_string(Simulation::air_absorption_alpha2_);
    std::filesystem::create_directories(dir_name);

    vector<shared_ptr<Partition>> partitions;
    vector<shared_ptr<SoundSource>> sources;
    vector<shared_ptr<Recorder>> recorders;

    string asset_path = "./assets/" + config.asset_name + ".txt";
    cout << "Loading assets from: " << asset_path << endl;
    partitions = Partition::ImportPartitions(asset_path);
    if (partitions.empty()) {
        cerr << "FATAL ERROR: No partitions loaded from " << asset_path << "! Check file path and content." << endl;
        return 1;
    }
    sources = SoundSource::ImportSources(asset_path);
    cout << "Loaded " << partitions.size() << " partitions and " << sources.size() << " sources." << endl;

    
    if (cli_args.mode == RunMode::SIM_RECORD_FIELD || cli_args.mode == RunMode::SIM_RECORD_RESPONSE) {
        recorders = Recorder::ImportRecorders(asset_path);
        for (auto r : recorders) r->FindPartition(partitions);
    }

    shared_ptr<Simulation> simulation = nullptr;
    if (cli_args.mode != RunMode::VIZ_RECORD) {
        simulation = make_shared<Simulation>(partitions, sources);
        simulation->Info();
    }

    // Determine dimensions
    int rxs = 1e9, rys = 1e9, rzs = 1e9, rxe = -1e9, rye = -1e9, rze = -1e9;
    for (auto p : partitions) {
        rxs = min(rxs, p->x_start_); rys = min(rys, p->y_start_); rzs = min(rzs, p->z_start_);
        rxe = max(rxe, p->x_end_);   rye = max(rye, p->y_end_);   rze = max(rze, p->z_end_);
    }
    int dim_x = rxe - rxs, dim_y = rye - rys, dim_z = rze - rzs;
    int pml = Simulation::n_pml_layers_;
    int panel_w_sim = max({ dim_x, dim_y, dim_z }) + 2 * pml;
    int panel_h_sim = panel_w_sim;
    int resolution_x = panel_w_sim * 3;
    int resolution_y = panel_h_sim;
    int panel_sz = config.fixed_panel_size;

    SDL_Window* window = nullptr; SDL_Renderer* renderer = nullptr; SDL_Texture* texture = nullptr;
    TTF_Font* Sans = nullptr; SDL_Rect simulation_rect, Message_rect, Message_rect2, label_rects[3], bar_max_label, bar_min_label, bar_zero_label;
    SDL_Color White = { 255, 255, 255 }, Gray = { 180, 180, 180 }, Red = { 255, 0, 0 }, Blue = { 0, 0, 255 };
    string panel_titles[3] = { "XY", "XZ", "YZ" };

    if (cli_args.mode == RunMode::SIM_VIZ || cli_args.mode == RunMode::VIZ_RECORD) {
        const int GAP = 12; // Spacing between panels
        SDL_Init(SDL_INIT_VIDEO); TTF_Init();
        int win_w = (panel_sz + GAP) * 2 + panel_sz + 80; int win_h = panel_sz + 48;
        window = SDL_CreateWindow("WASAbi 2.5D", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, win_w, win_h, 0);
        renderer = SDL_CreateRenderer(window, -1, 0);
        // CRITICAL: MUST BE RGB888 TO MATCH ORIGINAL COLORBAR AND REPO BEHAVIOR
        texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB888, SDL_TEXTUREACCESS_STREAMING, resolution_x, resolution_y);
        Sans = TTF_OpenFont("font/SourceSansPro-Regular.ttf", 48);
        Message_rect = { 4, win_h - 20, 120, 18 }; Message_rect2 = { win_w - 124, win_h - 20, 118, 18 };
        for (int i = 0; i < 3; i++) label_rects[i] = { i * (panel_sz + GAP) + 4, 2, 60, 20 };
        bar_max_label = { (panel_sz + GAP) * 2 + panel_sz + 35, 24, 40, 18 };
        bar_min_label = { (panel_sz + GAP) * 2 + panel_sz + 35, panel_sz + 24 - 18, 40, 18 };
        bar_zero_label = { (panel_sz + GAP) * 2 + panel_sz + 35, panel_sz / 2 + 24 - 9, 40, 18 };
    }

    bool quit = false; int time_step = 0; int total_steps = Simulation::duration_ / Simulation::dt_;
    SDL_Event event;

    if (cli_args.mode != RunMode::VIZ_RECORD) {
        while (!quit && time_step < total_steps) {
            while (SDL_PollEvent(&event)) if (event.type == SDL_QUIT) quit = true;
            time_step = simulation->Update();
            for (auto r : recorders) {
                if (is_record_response) r->RecordResponse(time_step);
                if (is_record_field) r->RecordField(time_step);
            }
            if (cli_args.mode == RunMode::SIM_VIZ && time_step % Simulation::viz_skip_ == 0) {
                const int GAP = 12;
                SDL_RenderClear(renderer);
                SDL_SetRenderDrawColor(renderer, 20, 20, 20, 255);
                SDL_Rect bg = { 0, 0, (panel_sz + GAP) * 2 + panel_sz + 80, panel_sz + 48 }; SDL_RenderFillRect(renderer, &bg);
                
                SDL_UpdateTexture(texture, nullptr, simulation->pixels().data(), simulation->render_w() * sizeof(Uint32));
                for (int i = 0; i < 3; i++) {
                    SDL_Rect src = { i * panel_w_sim, 0, panel_w_sim, panel_h_sim };
                    SDL_Rect dst = { i * (panel_sz + GAP), 24, panel_sz, panel_sz };
                    SDL_RenderCopy(renderer, texture, &src, &dst);
                }
                for (int i = 0; i < 3; i++) {
                    SDL_Surface* s = TTF_RenderText_Solid(Sans, panel_titles[i].c_str(), Gray);
                    SDL_Texture* t = SDL_CreateTextureFromSurface(renderer, s);
                    SDL_RenderCopy(renderer, t, nullptr, &label_rects[i]);
                    SDL_FreeSurface(s); SDL_DestroyTexture(t);
                }
                string m = to_string(time_step) + "/" + to_string(total_steps);
                SDL_Surface* sm = TTF_RenderText_Solid(Sans, m.c_str(), White);
                SDL_Texture* tm = SDL_CreateTextureFromSurface(renderer, sm);
                SDL_RenderCopy(renderer, tm, nullptr, &Message_rect);
                SDL_FreeSurface(sm); SDL_DestroyTexture(tm);
                for (int y = 0; y < panel_sz; y++) {
                    float n = 1.0f - (float)y / panel_sz; int r,g,b;
                    if (n >= 0.5f) { float p = (n-0.5f)*2.0f; r=255; g=255*(1-p); b=255*(1-p); }
                    else { float n2 = n*2.0f; r=255*n2; g=255*n2; b=255; }
                    SDL_SetRenderDrawColor(renderer, r, g, b, 255);
                    SDL_RenderDrawLine(renderer, (panel_sz + GAP) * 2 + panel_sz + 10, 24 + y, (panel_sz + GAP) * 2 + panel_sz + 30, 24 + y);
                }
                double cp = simulation->GetLastMaxPressure();
                stringstream ss; ss << fixed << setprecision(2) << cp;
                SDL_Surface *smax = TTF_RenderText_Solid(Sans, ("+"+ss.str()).c_str(), Red);
                SDL_Texture *tmax = SDL_CreateTextureFromSurface(renderer, smax);
                SDL_RenderCopy(renderer, tmax, nullptr, &bar_max_label);
                SDL_FreeSurface(smax); SDL_DestroyTexture(tmax);
                SDL_Surface *smin = TTF_RenderText_Solid(Sans, ("-"+ss.str()).c_str(), Blue);
                SDL_Texture *tmin = SDL_CreateTextureFromSurface(renderer, smin);
                SDL_RenderCopy(renderer, tmin, nullptr, &bar_min_label);
                SDL_FreeSurface(smin); SDL_DestroyTexture(tmin);
                SDL_Surface *szero = TTF_RenderText_Solid(Sans, "0.0", White);
                SDL_Texture *tzero = SDL_CreateTextureFromSurface(renderer, szero);
                SDL_RenderCopy(renderer, tzero, nullptr, &bar_zero_label);
                SDL_FreeSurface(szero); SDL_DestroyTexture(tzero);
                SDL_RenderPresent(renderer);
            } else if (cli_args.mode != RunMode::SIM_VIZ && time_step % (Simulation::viz_skip_*5) == 0) {
                cout << "Progress: " << time_step << "/" << total_steps << "\r"; cout.flush();
            }
        }
    } else {
        // VIZ_RECORD (simplified match to original look, black background)
        ifstream ifs(cli_args.playback_file); string line; 
        int pml_lay = Simulation::n_pml_layers_;
        int sx_p = sources.empty() ? dim_x/2 : (sources[0]->x() - rxs);
        int sy_p = sources.empty() ? dim_y/2 : (sources[0]->y() - rys);
        int sz_p = sources.empty() ? dim_z/2 : (sources[0]->z() - rzs);
        vector<double> fd(dim_x*dim_y*dim_z); vector<Uint32> pixels(resolution_x*resolution_y, 0x000000);
        float smooth_v = 0.0f; int fi = 0;
        while (getline(ifs, line) && !quit) {
            while (SDL_PollEvent(&event)) if (event.type == SDL_QUIT) quit = true;
            stringstream ss(line); double max_p = 0;
            for (int i=0; i<dim_x*dim_y*dim_z; i++) { ss >> fd[i]; max_p = max(max_p, abs(fd[i])); }
            float tv = 0.7f / (max_p + 1e-9); if (tv > Simulation::max_viz_gain_) tv = Simulation::max_viz_gain_;
            if (fi==0) smooth_v = tv; else smooth_v = smooth_v*0.9f + tv*0.1f;
            float v_c = max(0.001f, min(1000.0f, smooth_v));
            fill(pixels.begin(), pixels.end(), 0x000000);
            for (auto p : partitions) {
                int px1 = p->x_start_-rxs+pml_lay, px2 = p->x_end_-rxs+pml_lay, py1 = p->y_start_-rys+pml_lay, py2 = p->y_end_-rys+pml_lay, pz1 = p->z_start_-rzs+pml_lay, pz2 = p->z_end_-rzs+pml_lay;
                Uint32 c = p->should_render_ ? 0xFFFFFF : 0x808080;
                if (sz_p+pml_lay >= pz1 && sz_p+pml_lay < pz2) for (int j=py1; j<py2; j++) for (int i=px1; i<px2; i++) pixels[j*resolution_x+i] = c;
                if (sy_p+pml_lay >= py1 && sy_p+pml_lay < py2) for (int j=pz1; j<pz2; j++) for (int i=px1; i<px2; i++) pixels[j*resolution_x+i+panel_w_sim] = c;
                if (sx_p+pml_lay >= px1 && sx_p+pml_lay < px2) for (int j=pz1; j<pz2; j++) for (int i=py1; i<py2; i++) pixels[j*resolution_x+i+2*panel_w_sim] = c;
            }
            for (int j=0; j<dim_y; j++) for (int i=0; i<dim_x; i++) {
                double v = fd[sz_p*dim_y*dim_x + j*dim_x + i]; if (abs(v)>1e-12) pixels[(j+pml_lay)*resolution_x+(i+pml_lay)] = calculate_color_playback(v, v_c);
            }
            for (int k=0; k<dim_z; k++) for (int i=0; i<dim_x; i++) {
                double v = fd[k*dim_y*dim_x + sy_p*dim_x + i]; if (abs(v)>1e-12) pixels[(k+pml_lay)*resolution_x+(i+pml_lay+panel_w_sim)] = calculate_color_playback(v, v_c);
            }
            for (int k=0; k<dim_z; k++) for (int j=0; j<dim_y; j++) {
                double v = fd[k*dim_y*dim_x + j*dim_x + sx_p]; if (abs(v)>1e-12) pixels[(k+pml_lay)*resolution_x+(j+pml_lay+2*panel_w_sim)] = calculate_color_playback(v, v_c);
            }
            const int GAP = 12;
            SDL_RenderClear(renderer);
            SDL_SetRenderDrawColor(renderer, 20, 20, 20, 255);
            SDL_Rect bg = { 0, 0, (panel_sz + GAP) * 2 + panel_sz + 80, panel_sz + 48 }; SDL_RenderFillRect(renderer, &bg);

            SDL_UpdateTexture(texture, nullptr, pixels.data(), resolution_x*4);
            for (int i = 0; i < 3; i++) {
                SDL_Rect src = { i * panel_w_sim, 0, panel_w_sim, panel_h_sim };
                SDL_Rect dst = { i * (panel_sz + GAP), 24, panel_sz, panel_sz };
                SDL_RenderCopy(renderer, texture, &src, &dst);
            }
            for (int i = 0; i < 3; i++) {
                SDL_Surface* s = TTF_RenderText_Solid(Sans, panel_titles[i].c_str(), Gray);
                SDL_Texture* t = SDL_CreateTextureFromSurface(renderer, s); SDL_RenderCopy(renderer, t, nullptr, &label_rects[i]);
                SDL_FreeSurface(s); SDL_DestroyTexture(t);
            }
            for (int y = 0; y < panel_sz; y++) {
                float n = 1.0f - (float)y / panel_sz; int r,g,b;
                if (n >= 0.5f) { float p = (n-0.5f)*2.0f; r=255; g=255*(1-p); b=255*(1-p); }
                else { float n2 = n*2.0f; r=255*n2; g=255*n2; b=255; }
                SDL_SetRenderDrawColor(renderer, r, g, b, 255);
                SDL_RenderDrawLine(renderer, (panel_sz + GAP) * 2 + panel_sz + 10, 24 + y, (panel_sz + GAP) * 2 + panel_sz + 30, 24 + y);
            }
            stringstream ss_m; ss_m << fixed << setprecision(2) << max_p;
            SDL_Surface *smax = TTF_RenderText_Solid(Sans, ("+"+ss_m.str()).c_str(), Red);
            SDL_Texture *tmax = SDL_CreateTextureFromSurface(renderer, smax); SDL_RenderCopy(renderer, tmax, nullptr, &bar_max_label);
            SDL_FreeSurface(smax); SDL_DestroyTexture(tmax);
            SDL_Surface *smin = TTF_RenderText_Solid(Sans, ("-"+ss_m.str()).c_str(), Blue);
            SDL_Texture *tmin = SDL_CreateTextureFromSurface(renderer, smin); SDL_RenderCopy(renderer, tmin, nullptr, &bar_min_label);
            SDL_FreeSurface(smin); SDL_DestroyTexture(tmin);
            SDL_Surface *szero = TTF_RenderText_Solid(Sans, "0.0", White);
            SDL_Texture *tzero = SDL_CreateTextureFromSurface(renderer, szero); SDL_RenderCopy(renderer, tzero, nullptr, &bar_zero_label);
            SDL_FreeSurface(szero); SDL_DestroyTexture(tzero);
            SDL_RenderPresent(renderer); fi++; SDL_Delay(cli_args.playback_delay);
        }
    }

    if (window) { SDL_DestroyTexture(texture); SDL_DestroyRenderer(renderer); SDL_DestroyWindow(window); SDL_Quit(); TTF_Quit(); }
    cout << "\nSimulation finished. (" << omp_get_wtime() - time1 << " s)" << endl;
    return 0;
}
