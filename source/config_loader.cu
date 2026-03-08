#include "config_loader.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

Config load_config(const std::string& filename, const std::string& experiment_name) {
    Config config;
    config.experiment_name = experiment_name;
    
    // Explicit Default Fallbacks
    config.asset_name = "hall";
    // ... (rest of defaults)
    config.boundary_absorption = 1.0;
    config.air_absorption_alpha1 = 0.0;
    config.air_absorption_alpha2 = 0.0;
    config.duration = 0.2;
    config.c0 = 343.5;
    config.n_pml_layers = 5;
    config.precision = "coarse";
    config.is_record_response = false;
    config.is_record_field = false;
    
    config.viz_skip = 10;
    config.fixed_panel_size = 400;
    config.max_viz_gain = 100.0f;

    std::string actual_path = filename;
    if (!experiment_name.empty()) {
        actual_path = "./experiments/" + experiment_name + "/config.json";
    }

    std::ifstream file(actual_path);
    if (!file.is_open()) {
        cerr << "WARNING: Can't load " << actual_path << ". Using default parameters." << endl;
        return config;
    }

    try {
        json j;
        file >> j;
        
        if (j.contains("simulation")) {
            auto& s = j["simulation"];
            if (s.contains("asset_name")) config.asset_name = s["asset_name"];
            // ... (keep the existing mapping logic)
            if (s.contains("boundary_absorption")) config.boundary_absorption = s["boundary_absorption"];
            if (s.contains("air_absorption_alpha1")) config.air_absorption_alpha1 = s["air_absorption_alpha1"];
            if (s.contains("air_absorption_alpha2")) config.air_absorption_alpha2 = s["air_absorption_alpha2"];
            if (s.contains("duration")) config.duration = s["duration"];
            if (s.contains("c0")) config.c0 = s["c0"];
            if (s.contains("n_pml_layers")) config.n_pml_layers = s["n_pml_layers"];
            if (s.contains("precision")) config.precision = s["precision"];
            if (s.contains("is_record_response")) config.is_record_response = s["is_record_response"];
            if (s.contains("is_record_field")) config.is_record_field = s["is_record_field"];
        }
        
        if (j.contains("visualization")) {
            auto& v = j["visualization"];
            if (v.contains("viz_skip")) config.viz_skip = v["viz_skip"];
            if (v.contains("fixed_panel_size")) config.fixed_panel_size = v["fixed_panel_size"];
            if (v.contains("max_viz_gain")) config.max_viz_gain = v["max_viz_gain"];
        }
    } catch (const std::exception& e) {
        cerr << "FATAL ERROR parsing " << actual_path << ": " << e.what() << endl;
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


