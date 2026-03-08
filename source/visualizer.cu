#include "visualizer.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

Visualizer::Visualizer(RunMode mode, int panel_sz, int panel_w_sim, int panel_h_sim, int resolution_x, int resolution_y)
    : window_(nullptr), renderer_(nullptr), texture_(nullptr), font_(nullptr),
      panel_sz_(panel_sz), panel_w_sim_(panel_w_sim), panel_h_sim_(panel_h_sim),
      resolution_x_(resolution_x), resolution_y_(resolution_y) 
{
    if (mode == RunMode::SIM_VIZ || mode == RunMode::VIZ_RECORD) {
        SDL_Init(SDL_INIT_VIDEO); 
        TTF_Init();
        int win_w = (panel_sz_ + GAP) * 2 + panel_sz_ + 80; 
        int win_h = panel_sz_ + 48;
        
        window_ = SDL_CreateWindow("WASAbi 2.5D", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, win_w, win_h, 0);
        renderer_ = SDL_CreateRenderer(window_, -1, 0);
        texture_ = SDL_CreateTexture(renderer_, SDL_PIXELFORMAT_RGB888, SDL_TEXTUREACCESS_STREAMING, resolution_x_, resolution_y_);
        font_ = TTF_OpenFont("font/SourceSansPro-Regular.ttf", 48);
        
        message_rect_ = { 4, win_h - 20, 120, 18 }; 
        message_rect2_ = { win_w - 124, win_h - 20, 118, 18 };
        for (int i = 0; i < 3; i++) label_rects_[i] = { i * (panel_sz_ + GAP) + 4, 2, 60, 20 };
        bar_max_label_ = { (panel_sz_ + GAP) * 2 + panel_sz_ + 35, 24, 40, 18 };
        bar_min_label_ = { (panel_sz_ + GAP) * 2 + panel_sz_ + 35, panel_sz_ + 24 - 18, 40, 18 };
        bar_zero_label_ = { (panel_sz_ + GAP) * 2 + panel_sz_ + 35, panel_sz_ / 2 + 24 - 9, 40, 18 };
    }
}

Visualizer::~Visualizer() {
    if (font_) TTF_CloseFont(font_);
    if (texture_) SDL_DestroyTexture(texture_);
    if (renderer_) SDL_DestroyRenderer(renderer_);
    if (window_) SDL_DestroyWindow(window_);
    TTF_Quit();
    SDL_Quit();
}

void Visualizer::HandleEvents(bool& quit) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_QUIT) {
            quit = true;
        }
    }
}

void Visualizer::RenderUI(int time_step, int total_steps, double max_p) {
    // 1. Draw 3 Sub-panels from the large simulation texture
    for (int i = 0; i < 3; i++) {
        SDL_Rect src = { i * panel_w_sim_, 0, panel_w_sim_, panel_h_sim_ };
        SDL_Rect dst = { i * (panel_sz_ + GAP), 24, panel_sz_, panel_sz_ };
        SDL_RenderCopy(renderer_, texture_, &src, &dst);
    }
    
    // 2. Draw panel labels (XY, XZ, YZ)
    for (int i = 0; i < 3; i++) {
        SDL_Surface* s = TTF_RenderText_Solid(font_, panel_titles_[i].c_str(), color_gray_);
        SDL_Texture* t = SDL_CreateTextureFromSurface(renderer_, s);
        SDL_RenderCopy(renderer_, t, nullptr, &label_rects_[i]);
        SDL_FreeSurface(s); 
        SDL_DestroyTexture(t);
    }
    
    // 3. Draw progress text
    string m = to_string(time_step) + "/" + to_string(total_steps);
    SDL_Surface* sm = TTF_RenderText_Solid(font_, m.c_str(), color_white_);
    SDL_Texture* tm = SDL_CreateTextureFromSurface(renderer_, sm);
    SDL_RenderCopy(renderer_, tm, nullptr, &message_rect_);
    SDL_FreeSurface(sm); 
    SDL_DestroyTexture(tm);
    
    // 4. Draw colorbar gradient manually
    for (int y = 0; y < panel_sz_; y++) {
        float n = 1.0f - (float)y / panel_sz_; 
        int r, g, b;
        if (n >= 0.5f) { float p = (n - 0.5f) * 2.0f; r = 255; g = 255 * (1 - p); b = 255 * (1 - p); }
        else { float n2 = n * 2.0f; r = 255 * n2; g = 255 * n2; b = 255; }
        SDL_SetRenderDrawColor(renderer_, r, g, b, 255);
        SDL_RenderDrawLine(renderer_, (panel_sz_ + GAP) * 2 + panel_sz_ + 10, 24 + y, (panel_sz_ + GAP) * 2 + panel_sz_ + 30, 24 + y);
    }
    
    // 5. Draw colorbar text labels
    stringstream ss; 
    ss << fixed << setprecision(2) << max_p;
    
    SDL_Surface* smax = TTF_RenderText_Solid(font_, ("+" + ss.str()).c_str(), color_red_);
    SDL_Texture* tmax = SDL_CreateTextureFromSurface(renderer_, smax);
    SDL_RenderCopy(renderer_, tmax, nullptr, &bar_max_label_);
    SDL_FreeSurface(smax); SDL_DestroyTexture(tmax);
    
    SDL_Surface* smin = TTF_RenderText_Solid(font_, ("-" + ss.str()).c_str(), color_blue_);
    SDL_Texture* tmin = SDL_CreateTextureFromSurface(renderer_, smin);
    SDL_RenderCopy(renderer_, tmin, nullptr, &bar_min_label_);
    SDL_FreeSurface(smin); SDL_DestroyTexture(tmin);
    
    SDL_Surface* szero = TTF_RenderText_Solid(font_, "0.0", color_white_);
    SDL_Texture* tzero = SDL_CreateTextureFromSurface(renderer_, szero);
    SDL_RenderCopy(renderer_, tzero, nullptr, &bar_zero_label_);
    SDL_FreeSurface(szero); SDL_DestroyTexture(tzero);
}

void Visualizer::RenderSimulationFrame(int time_step, int total_steps, std::shared_ptr<Simulation> simulation) {
    if (!IsValid()) return;
    
    SDL_RenderClear(renderer_);
    SDL_SetRenderDrawColor(renderer_, 20, 20, 20, 255); // Background
    SDL_Rect bg = { 0, 0, (panel_sz_ + GAP) * 2 + panel_sz_ + 80, panel_sz_ + 48 }; 
    SDL_RenderFillRect(renderer_, &bg);
    
    // Update central texture with CUDA calculation output
    SDL_UpdateTexture(texture_, nullptr, simulation->pixels().data(), simulation->render_w() * sizeof(Uint32));
    
    RenderUI(time_step, total_steps, simulation->GetLastMaxPressure());
    
    SDL_RenderPresent(renderer_);
}

void Visualizer::RenderPlaybackFrame(int time_step, int total_steps, double max_p, const std::vector<Uint32>& pixels) {
    if (!IsValid()) return;
    
    SDL_RenderClear(renderer_);
    SDL_SetRenderDrawColor(renderer_, 20, 20, 20, 255); // Background
    SDL_Rect bg = { 0, 0, (panel_sz_ + GAP) * 2 + panel_sz_ + 80, panel_sz_ + 48 }; 
    SDL_RenderFillRect(renderer_, &bg);
    
    // Update central texture with pre-calculated playback frame
    SDL_UpdateTexture(texture_, nullptr, pixels.data(), resolution_x_ * 4);
    
    RenderUI(time_step, total_steps, max_p);
    
    SDL_RenderPresent(renderer_);
}

Uint32 Visualizer::CalculateColorPlayback(double p, float v_coef) {
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
