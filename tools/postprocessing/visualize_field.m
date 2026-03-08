clear all, close all, clc;

% Default experiment path
if ~exist('experiment_path', 'var')
    experiment_path = '../../experiments/hall';
end

data_filename = [experiment_path, '/output/record_data_0.bin'];
% Load Experiment Config
config_str = fileread([experiment_path, '/config.json']);
config = jsondecode(config_str);

% Load Assets
asset_str = fileread([experiment_path, '/asset.json']);
asset = jsondecode(asset_str);

% Simulation Params
dh = config.simulation.dh;
n_pml = config.simulation.n_pml_layers;

% Calculate Grid Size from Partitions (Bounding Box)
partitions = asset.partitions;
x_min = inf; x_max = -inf;
y_min = inf; y_max = -inf;
z_min = inf; z_max = -inf;

for i = 1:length(partitions)
    p = partitions(i);
    % Partitions in JSON are in WORLD meters. Convert to grid indices.
    x_min = min(x_min, floor(p.x / dh));
    x_max = max(x_max, ceil((p.x + p.w) / dh));
    y_min = min(y_min, floor(p.y / dh));
    y_max = max(y_max, ceil((p.y + p.h) / dh));
    z_min = min(z_min, floor(p.z / dh));
    z_max = max(z_max, ceil((p.z + p.d) / dh));
end

% Total dimensions (excluding PML, as they are not recorded)
Lx = x_max - x_min;
Ly = y_max - y_min;
Lz = z_max - z_min;

% Stream Binary Chunk
fileID = fopen(data_filename, 'r');
A = fread(fileID, inf, 'double');
fclose(fileID);

% Reshape the data into a 4D array (Time, X, Y, Z)
% We calculate the number of time steps (N) based on the total elements
points_per_frame = Lx * Ly * Lz;
N = round(length(A) / points_per_frame);

pressure_values = zeros(N, Lx, Ly, Lz);

for iT = 1:N
    start_idx = (iT - 1) * points_per_frame + 1;
    end_idx = iT * points_per_frame;
    pressure_values(iT,:,:,:) = reshape(A(start_idx:end_idx), Lx, Ly, Lz);
end

%%
close all;

for iT = 50:1:54
    figure()
    p = squeeze(pressure_values(iT,:,:,:));
    
    slice(p, Ly*5/6, Lx/2, Lz*3/10);  % Create a slice plot
    
    shading interp
    % title(['Pressure Distribution at Time Instant ', num2str(time_instants(iT))]);
    daspect([1, 1, 1]); % Equal aspect ratio for x, y, and z dimensions
    xlabel('y')
    ylabel('x')
    zlabel('z')
    clim([-1 1])
    colormap jet
    % Change font size
    set(gca, 'FontSize', 24)
end
drawnow;