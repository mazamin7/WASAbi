#ifndef VISUALIZER_H
#define VISUALIZER_H

#include <SDL.h>
#undef main
#include <SDL_ttf.h>
#include <string>
#include <vector>
#include "cli_parser.h"
#include "simulation.h"

class Visualizer {
public:
    Visualizer(RunMode mode, int panel_sz, int panel_w_sim, int panel_h_sim, int resolution_x, int resolution_y);
    ~Visualizer();

    bool IsValid() const { return window_ != nullptr && renderer_ != nullptr && texture_ != nullptr && font_ != nullptr; }
    void HandleEvents(bool& quit);
    void RenderSimulationFrame(int time_step, int total_steps, std::shared_ptr<Simulation> simulation);
    void RenderPlaybackFrame(int time_step, int total_steps, double max_p, const std::vector<Uint32>& pixels);
    
    static Uint32 CalculateColorPlayback(double p, float v_coef);

private:
    void RenderUI(int time_step, int total_steps, double max_p);

    SDL_Window* window_;
    SDL_Renderer* renderer_;
    SDL_Texture* texture_;
    TTF_Font* font_;

    SDL_Rect message_rect_;
    SDL_Rect message_rect2_;
    SDL_Rect label_rects_[3];
    SDL_Rect bar_max_label_;
    SDL_Rect bar_min_label_;
    SDL_Rect bar_zero_label_;

    std::string panel_titles_[3] = { "XY", "XZ", "YZ" };
    SDL_Color color_white_ = { 255, 255, 255 };
    SDL_Color color_gray_ = { 180, 180, 180 };
    SDL_Color color_red_ = { 255, 0, 0 };
    SDL_Color color_blue_ = { 0, 0, 255 };

    int panel_sz_;
    int panel_w_sim_;
    int panel_h_sim_;
    int resolution_x_;
    int resolution_y_;
    const int GAP = 12;
};

#endif // VISUALIZER_H
