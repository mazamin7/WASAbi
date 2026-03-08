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
    config.n_pml_layers = 5;
    config.dh = 0.5;
    config.dt = 0.000625;
    
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
            if (s.contains("duration")) config.duration = s["duration"];
            if (s.contains("n_pml_layers")) config.n_pml_layers = s["n_pml_layers"];
            
            // Override with explicit dh/dt if provided
            if (s.contains("dh")) config.dh = s["dh"];
            if (s.contains("dt")) config.dt = s["dt"];
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


