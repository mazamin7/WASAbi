#include "config_loader.h"
#include <iostream>
#include <filesystem>
#include "ini.h"
#include <cstring>

using namespace std;

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

Config load_config(const std::string& filename) {
    Config config;
    config.viz_skip = 10;
    config.fixed_panel_size = 400;
    config.max_viz_gain = 100.0f;
    if (ini_parse(filename.c_str(), parse_ini_handler, &config) < 0) {
        cerr << "Can't load " << filename << endl;
    }
    return config;
}

void set_precision_params(const std::string& precision, double& dh, double& dt) {
    if (precision == "coarse") { dh = 0.5; dt = 6.25e-4; }
    else if (precision == "fine") { dh = 0.2; dt = 2e-4; }
    else if (precision == "finer") { dh = 0.1; dt = 1.25e-4; }
    else if (precision == "finest") { dh = 0.05; dt = 0.625e-4; }
    else { dh = 0.2; dt = 2e-4; }
}

void ensureConfigExists(const std::string& config_path, const std::string& default_path) {
    if (!std::filesystem::exists(config_path)) {
        cout << "Config not found at " << config_path << ". Attempting to copy " << default_path << "..." << endl;
        if (std::filesystem::exists(default_path)) {
            try {
                std::filesystem::copy(default_path, config_path);
                cout << "Successfully created " << config_path << " from default." << endl;
            } catch (const std::exception& e) {
                cerr << "Failed to copy default config: " << e.what() << endl;
            }
        } else {
            cerr << "CRITICAL: Default config not found at " << default_path << endl;
        }
    }
}
