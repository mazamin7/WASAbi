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
        texture_ = SDL_CreateTexture(renderer_, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, resolution_x_, resolution_y_);
        font_ = TTF_OpenFont("source/font/SourceSansPro-Regular.ttf", 48);
        
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
    DrawMarkers();
    
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
    DrawMarkers();
    
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
    // Color normalization and shift to match SDL_PIXELFORMAT_RGBA8888 (A:24-31, B:16-23, G:8-15, R:0-7)
    // To match the CUDA kernel RGBAToUint32: (a << 24) | (b << 16) | (g << 8) | r
    uint32_t a = 255;
    return (a << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
}

void Visualizer::SetMarkers(const std::vector<Marker>& sources, const std::vector<Marker>& receivers, 
                            int rxs, int rys, int rzs, int pml) {
    source_markers_ = sources;
    receiver_markers_ = receivers;
    rxs_ = rxs; rys_ = rys; rzs_ = rzs; pml_ = pml;
}

void Visualizer::DrawCross(int x, int y, int size, SDL_Color color) {
    SDL_SetRenderDrawColor(renderer_, color.r, color.g, color.b, 255);
    SDL_RenderDrawLine(renderer_, x - size, y - size, x + size, y + size);
    SDL_RenderDrawLine(renderer_, x - size, y + size, x + size, y - size);
}

void Visualizer::DrawCircle(int x, int y, int radius, SDL_Color color) {
    SDL_SetRenderDrawColor(renderer_, color.r, color.g, color.b, 255);
    int offsetx, offsety, d;
    offsetx = 0;
    offsety = radius;
    d = radius - 1;

    while (offsety >= offsetx) {
        SDL_RenderDrawPoint(renderer_, x + offsetx, y + offsety);
        SDL_RenderDrawPoint(renderer_, x + offsety, y + offsetx);
        SDL_RenderDrawPoint(renderer_, x - offsetx, y + offsety);
        SDL_RenderDrawPoint(renderer_, x - offsety, y + offsetx);
        SDL_RenderDrawPoint(renderer_, x + offsetx, y - offsety);
        SDL_RenderDrawPoint(renderer_, x + offsety, y - offsetx);
        SDL_RenderDrawPoint(renderer_, x - offsetx, y - offsety);
        SDL_RenderDrawPoint(renderer_, x - offsety, y - offsetx);

        if (d >= 2 * offsetx) {
            d -= 2 * offsetx + 1;
            offsetx += 1;
        } else if (d < 2 * (radius - offsety)) {
            d += 2 * offsety - 1;
            offsety -= 1;
        } else {
            d += 2 * (offsety - offsetx - 1);
            offsety -= 1;
            offsetx += 1;
        }
    }
}

void Visualizer::DrawMarkers() {
    if (!IsValid()) return;
    
    // Convert simulation pixel coordinate space to actual screen space
    float scale = (float)panel_sz_ / (float)panel_w_sim_;
    
    // Common colors: Sources = Green Crosses, Receivers = Cyan Circles
    SDL_Color src_color = { 50, 255, 50 };
    SDL_Color rec_color = { 0, 255, 255 };

    int marker_size = max(2, (int)(3.0f * scale)); // Scale marker size based on panel resolution
    
    auto draw_markers_for_list = [&](const std::vector<Marker>& markers, SDL_Color color, bool is_source) {
        for (const auto& marker : markers) {
            // Map global coords to simulation field coords
            int sim_x = marker.x - rxs_ + pml_;
            int sim_y = marker.y - rys_ + pml_;
            int sim_z = marker.z - rzs_ + pml_;

            // XY Panel (sub-panel 0): Panel logic renders x-axis along panel_width, y-axis along panel_height
            int screen_xy_x = (0 * (panel_sz_ + GAP)) + (sim_x * scale);
            int screen_xy_y = 24 + (sim_y * scale);
            
            // XZ Panel (sub-panel 1): Panel logic renders x-axis along panel_width, z-axis along panel_height
            int screen_xz_x = (1 * (panel_sz_ + GAP)) + (sim_x * scale);
            int screen_xz_y = 24 + (sim_z * scale);
            
            // YZ Panel (sub-panel 2): Panel logic renders y-axis along panel_width, z-axis along panel_height
            int screen_yz_x = (2 * (panel_sz_ + GAP)) + (sim_y * scale);
            int screen_yz_y = 24 + (sim_z * scale);
            
            // Draw marker shapes based on type
            if (is_source) {
                DrawCross(screen_xy_x, screen_xy_y, marker_size, color);
                DrawCross(screen_xz_x, screen_xz_y, marker_size, color);
                DrawCross(screen_yz_x, screen_yz_y, marker_size, color);
            } else {
                DrawCircle(screen_xy_x, screen_xy_y, marker_size, color);
                DrawCircle(screen_xz_x, screen_xz_y, marker_size, color);
                DrawCircle(screen_yz_x, screen_yz_y, marker_size, color);
            }
        }
    };
    
    // Draw all collected markers
    draw_markers_for_list(source_markers_, src_color, true);
    draw_markers_for_list(receiver_markers_, rec_color, false);
}
