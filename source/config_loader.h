#ifndef CONFIG_LOADER_H
#define CONFIG_LOADER_H

#include <string>

struct Config {
    std::string asset_name;
    double boundary_absorption;
    double air_absorption_alpha1;
    double air_absorption_alpha2;
    double duration;
    double c0;
    int n_pml_layers;
    std::string precision;
    int viz_skip;
    int fixed_panel_size;
    float max_viz_gain;
};

Config load_config(const std::string& filename);
void set_precision_params(const std::string& precision, double& dh, double& dt);
void ensureConfigExists(const std::string& config_path, const std::string& default_path);

#endif // CONFIG_LOADER_H
