#ifndef CONFIG_LOADER_H
#define CONFIG_LOADER_H

#include <string>

struct Config {
    std::string experiment_name;
    std::string asset_name;
    double duration;
    double c0;
    int n_pml_layers;
    double dh;

    double dt;
    int viz_skip;
    int fixed_panel_size;
    float max_viz_gain;
};

Config load_config(const std::string& filename, const std::string& experiment_name = "");

#endif // CONFIG_LOADER_H
