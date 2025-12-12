function [Nx, Ny, Nz] = calculate_grid_size(filename, dh)
% CALCULATE_GRID_SIZE Determines the total grid size (Nx x Ny x Nz) 
% from a file containing partition dimensions and a given spatial resolution (dh).
%
% Inputs:
%   filename - String, path to the text file containing partition data.
%   dh       - Double, the spatial resolution (delta h).
%
% Output:
%   Nx, Ny, Nz - Integers, the number of nodes in each dimension.
%
% File Format (Example Line):
%   x_start y_start z_start width height depth

    % --- 1. Load Data ---
    try
        % Read all data into a matrix
        data = load(filename);
    catch
        error('Error: Could not load data from file "%s". Check file path and format.', filename);
    end

    % Check if dh is valid
    if dh <= 0
        error('Error: Spatial resolution (dh) must be positive and non-zero.');
    end

    % --- 2. Determine Overall Physical Extent ---
    
    % The file columns are: [x_start, y_start, z_start, width, height, depth]
    x_starts = data(:, 1);
    y_starts = data(:, 2);
    z_starts = data(:, 3);
    widths   = data(:, 4);
    heights  = data(:, 5);
    depths   = data(:, 6);

    % Calculate the end coordinates for all partitions
    x_ends = x_starts + widths;
    y_ends = y_starts + heights;
    z_ends = z_starts + depths;

    % Find the global minimum and maximum physical coordinates
    X_start_global = min(x_starts);
    Y_start_global = min(y_starts);
    Z_start_global = min(z_starts);

    X_end_global = max(x_ends);
    Y_end_global = max(y_ends);
    Z_end_global = max(z_ends);

    % Calculate the total physical lengths of the simulation domain
    Lx = X_end_global - X_start_global;
    Ly = Y_end_global - Y_start_global;
    Lz = Z_end_global - Z_start_global;

    % --- 3. Calculate Number of Nodes ---
    
    % The number of nodes is (Length / dh) + 1, rounded to the nearest integer 
    % to handle floating point precision errors, assuming L must be a multiple of dh.
    
    Nx = round(Lx / dh);
    Ny = round(Ly / dh);
    Nz = round(Lz / dh);
    
    % Convert to integer type, as grid size must be discrete
    Nx = int32(Nx);
    Ny = int32(Ny);
    Nz = int32(Nz);

    % --- 4. Display Results ---
    disp('--- Global Simulation Extent ---');
    fprintf('Total physical length (Lx x Ly x Lz): %.2f x %.2f x %.2f\n', Lx, Ly, Lz);
    fprintf('Spatial Resolution (dh): %.4f\n', dh);
    disp('--------------------------------');
    fprintf('Calculated Grid Size (Nx x Ny x Nz): %d x %d x %d\n', Nx, Ny, Nz);

end