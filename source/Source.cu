/* WASAbi 2.5D
 * 
 * This is the entrance of the program. 
 * SDL is used as the interface to show the wave propagation.
 */

#include <iostream>
#include <omp.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include "json.hpp"

#undef main
#include "simulation.h"
#include "partition.h"
#include "boundary.h"
#include "sound_source.h"
#include "gaussian_source.h"
#include "recorder.h"

// New modular architecture
#include "config_loader.h"
#include "cli_parser.h"
#include "visualizer.h"

using namespace std;

bool is_record_response = false;
bool is_record_field = false;

// Global settings defaults
double Simulation::air_absorption_alpha1_ = 0.0;
double Simulation::air_absorption_alpha2_ = 1e-6;
double Simulation::duration_ = 2e-2;
double Simulation::c0_ = 343.5;
double Simulation::dh_ = 0.2;
double Simulation::dt_ = 2e-4;
int Simulation::n_pml_layers_ = 5;
int Simulation::viz_skip_ = 10;
float Simulation::max_viz_gain_ = 100.0f;

int main(int argc, char* argv[]) {
    // 1. Shift current directory to root if running from build
    try {
        std::filesystem::path cp = std::filesystem::current_path();
        if (cp.filename() == "build") {
            std::filesystem::current_path("..");
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

    // 3. Resolve active config
    std::string active_config = cli_args.config_path;
    
    // Default to 'hall' experiment if nothing is specified
    if (cli_args.experiment_name.empty() && active_config.empty()) {
        cli_args.experiment_name = "hall";
    }

    Config config = load_config(active_config, cli_args.experiment_name);

    Simulation::duration_ = config.duration;
    Simulation::n_pml_layers_ = config.n_pml_layers;
    Simulation::dh_ = config.dh;
    Simulation::dt_ = config.dt;
    Simulation::viz_skip_ = config.viz_skip;
    Simulation::max_viz_gain_ = config.max_viz_gain;

    if (cli_args.mode == RunMode::SIM_RECORD_FIELD) { is_record_field = true; is_record_response = false; }
    else if (cli_args.mode == RunMode::SIM_RECORD_RESPONSE) { is_record_field = false; is_record_response = true; }
    else { is_record_field = false; is_record_response = false; }

    double time1 = omp_get_wtime();
    
    // Determine the experiment root directory by searching parent directories
    std::filesystem::path exp_root;
    bool found_root = false;
    std::filesystem::path current_search = std::filesystem::current_path();
    
    // Check up to 3 levels up for the experiments directory
    for (int i = 0; i < 4; ++i) {
        cout << "Searching for experiments in: " << (current_search / "experiments").string() << endl;
        if (std::filesystem::exists(current_search / "experiments")) {
            exp_root = current_search / "experiments";
            found_root = true;
            cout << "Found experiments root at: " << exp_root.string() << endl;
            break;
        }
        if (current_search.has_parent_path()) {
            current_search = current_search.parent_path();
        } else {
            break;
        }
    }

    if (!found_root) {
        cerr << "WARNING: Could not find 'experiments' directory in current or parent folders." << endl;
        exp_root = "./experiments"; // Fallback to current dir
    }

    string dir_name;
    if (!config.experiment_name.empty()) {
        dir_name = (exp_root / config.experiment_name / "output").string();
    } else {
        dir_name = "./output";
    }
    std::filesystem::create_directories(dir_name);

    vector<shared_ptr<Partition>> partitions;
    vector<shared_ptr<SoundSource>> sources;
    vector<shared_ptr<Recorder>> recorders;

    string asset_path;
    if (!config.experiment_name.empty()) {
        asset_path = (exp_root / config.experiment_name / "asset.json").string();
    } else {
        asset_path = "./assets/" + config.asset_name + ".txt";
    }
    cout << "Loading assets from: " << asset_path << endl;
    partitions = Partition::ImportPartitions(asset_path);
    if (partitions.empty()) {
        cerr << "FATAL ERROR: No partitions loaded from " << asset_path << "! Check file path and content." << endl;
        return 1;
    }
    sources = SoundSource::ImportSources(asset_path, dir_name);
    cout << "Loaded " << partitions.size() << " partitions and " << sources.size() << " sources." << endl;

    // Unconditionally load recorders to know their positions for the visualizer
    recorders = Recorder::ImportRecorders(config, asset_path, Simulation::duration_ / Simulation::dt_, dir_name);
    for (auto r : recorders) r->FindPartition(partitions);

    shared_ptr<Simulation> simulation = nullptr;
    if (cli_args.mode != RunMode::VIZ_RECORD) {
        simulation = make_shared<Simulation>(partitions, sources);
        simulation->Info();
    }

    // Determine dimensions for Visualization
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

    Visualizer visualizer(cli_args.mode, panel_sz, panel_w_sim, panel_h_sim, resolution_x, resolution_y);
    
    // Pass Marker locations to Visualizer
    vector<Visualizer::Marker> src_markers;
    for (auto s : sources) src_markers.push_back({ s->x(), s->y(), s->z() });
    
    vector<Visualizer::Marker> rec_markers;
    for (auto r : recorders) rec_markers.push_back({ r->x(), r->y(), r->z() });
    
    visualizer.SetMarkers(src_markers, rec_markers, rxs, rys, rzs, pml);

    bool quit = false; 
    int time_step = 0; 
    int total_steps = Simulation::duration_ / Simulation::dt_;

    // --- MAIN SIMULATION LOOP ---
    if (cli_args.mode != RunMode::VIZ_RECORD) {
        while (!quit && time_step < total_steps) {
            visualizer.HandleEvents(quit);
            time_step = simulation->Update();
            
            for (auto r : recorders) {
                if (is_record_response) r->RecordResponse(time_step);
                if (is_record_field) r->RecordField(time_step);
            }
            
            if (cli_args.mode == RunMode::SIM_VIZ && Simulation::viz_skip_ > 0 && time_step % Simulation::viz_skip_ == 0) {
                visualizer.RenderSimulationFrame(time_step, total_steps, simulation);
            } else if (cli_args.mode != RunMode::SIM_VIZ && Simulation::viz_skip_ > 0 && time_step % (Simulation::viz_skip_ * 5) == 0) {
                cout << "Progress: " << time_step << "/" << total_steps << "\r"; cout.flush();
            } else if (Simulation::viz_skip_ == 0 && time_step % 50 == 0) {
                cout << "Progress: " << time_step << "/" << total_steps << "\r"; cout.flush();
            }
        }
        // Ensure all buffered response data is written to disk
        for (auto r : recorders) {
            if (is_record_response) r->FlushResponse();
        }
    } 
    // --- PLAYBACK LOOP ---
    else {
        // Read directly from binary stream
        int pml_lay = Simulation::n_pml_layers_;
        int sx_p = sources.empty() ? dim_x / 2 : (sources[0]->x() - rxs);
        int sy_p = sources.empty() ? dim_y / 2 : (sources[0]->y() - rys);
        int sz_p = sources.empty() ? dim_z / 2 : (sources[0]->z() - rzs);

        vector<double> fd(dim_x * dim_y * dim_z); 
        size_t frame_bytes = fd.size() * sizeof(double);
        vector<Uint32> pixels(resolution_x * resolution_y, 0x000000);
        float smooth_v = 0.0f; 
        int fi = 0;

        ifstream ifs(cli_args.playback_file, std::ios::binary); 
        if (!ifs.is_open()) {
            cerr << "FATAL ERROR: Could not open playback file " << cli_args.playback_file << endl;
            return 1;
        }
        // Read directly from binary stream
        while (ifs.read(reinterpret_cast<char*>(fd.data()), frame_bytes) && !quit) {
            visualizer.HandleEvents(quit);
            
            double max_p = 0;
            for (size_t i = 0; i < fd.size(); i++) { 
                max_p = max(max_p, abs(fd[i])); 
            }
            
            float tv = 0.7f / (max_p + 1e-9); 
            if (tv > Simulation::max_viz_gain_) tv = Simulation::max_viz_gain_;
            if (fi == 0) smooth_v = tv; else smooth_v = smooth_v * 0.9f + tv * 0.1f;
            float v_c = max(0.001f, min(1000.0f, smooth_v));
            
            fill(pixels.begin(), pixels.end(), 0x000000);
            
            for (auto p : partitions) {
                int px1 = p->x_start_ - rxs + pml_lay, px2 = p->x_end_ - rxs + pml_lay;
                int py1 = p->y_start_ - rys + pml_lay, py2 = p->y_end_ - rys + pml_lay;
                int pz1 = p->z_start_ - rzs + pml_lay, pz2 = p->z_end_ - rzs + pml_lay;
                Uint32 c = p->should_render_ ? 0xFFFFFF : 0x808080;
                
                if (sz_p + pml_lay >= pz1 && sz_p + pml_lay < pz2) 
                    for (int j = py1; j < py2; j++) for (int i = px1; i < px2; i++) pixels[j * resolution_x + i] = c;
                if (sy_p + pml_lay >= py1 && sy_p + pml_lay < py2) 
                    for (int j = pz1; j < pz2; j++) for (int i = px1; i < px2; i++) pixels[j * resolution_x + i + panel_w_sim] = c;
                if (sx_p + pml_lay >= px1 && sx_p + pml_lay < px2) 
                    for (int j = pz1; j < pz2; j++) for (int i = py1; i < py2; i++) pixels[j * resolution_x + i + 2 * panel_w_sim] = c;
            }
            
            for (int j = 0; j < dim_y; j++) for (int i = 0; i < dim_x; i++) {
                double v = fd[sz_p * dim_y * dim_x + j * dim_x + i]; 
                if (abs(v)>1e-12) pixels[(j + pml_lay) * resolution_x + (i + pml_lay)] = Visualizer::CalculateColorPlayback(v, v_c);
            }
            for (int k = 0; k < dim_z; k++) for (int i = 0; i < dim_x; i++) {
                double v = fd[k * dim_y * dim_x + sy_p * dim_x + i]; 
                if (abs(v)>1e-12) pixels[(k + pml_lay) * resolution_x + (i + pml_lay + panel_w_sim)] = Visualizer::CalculateColorPlayback(v, v_c);
            }
            for (int k = 0; k < dim_z; k++) for (int j = 0; j < dim_y; j++) {
                double v = fd[k * dim_y * dim_x + j * dim_x + sx_p]; 
                if (abs(v)>1e-12) pixels[(k + pml_lay) * resolution_x + (j + pml_lay + 2 * panel_w_sim)] = Visualizer::CalculateColorPlayback(v, v_c);
            }
            
            visualizer.RenderPlaybackFrame(fi, total_steps, max_p, pixels);
            if (fi % 50 == 0) { cout << "Playback Progress: " << fi << "/" << total_steps << "\r"; cout.flush(); }
            fi++; 
            SDL_Delay(cli_args.playback_delay);
        }
    }

    cout << "\nSimulation finished. (" << omp_get_wtime() - time1 << " s)" << endl;
    return 0;
}
