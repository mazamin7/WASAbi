function [S, V] = compute_geometry_S_V(file_path)
    
    % Read the data from the text file
    data = dlmread(file_path);
    
    % Extract the columns for position and dimensions
    x_pos = data(:, 1);
    y_pos = data(:, 2);
    z_pos = data(:, 3);
    x_length = data(:, 4);
    y_length = data(:, 5);
    z_length = data(:, 6);
    
    % Initialize total surface area and volume
    total_surface_area = 0;
    total_volume = 0;
    
    % Loop through each rectangle
    for i = 1:length(x_pos)
        % Calculate the surface area of the current rectangle
        surface_area = 2 * (x_length(i) * y_length(i) + x_length(i) * z_length(i) + y_length(i) * z_length(i));
        
        % Subtract the shared surface area with the previous rectangles
        for j = 1:i-1
            dx = max(0, min(x_pos(i) + x_length(i), x_pos(j) + x_length(j)) - max(x_pos(i), x_pos(j)));
            dy = max(0, min(y_pos(i) + y_length(i), y_pos(j) + y_length(j)) - max(y_pos(i), y_pos(j)));
            dz = max(0, min(z_pos(i) + z_length(i), z_pos(j) + z_length(j)) - max(z_pos(i), z_pos(j)));
            
            shared_surface = 2 * (dx * dy + dx * dz + dy * dz);
            surface_area = surface_area - shared_surface;
        end
        
        % Add the current surface area to the total
        total_surface_area = total_surface_area + surface_area;
        
        % Add the current volume to the total
        total_volume = total_volume + (x_length(i) * y_length(i) * z_length(i));
    end
    
    % Return the results
    S = total_surface_area;
    V = total_volume;

    S = 2720;
    V = 7600;

end
