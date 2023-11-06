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
    total_perimeter = 0;
    total_surface_area = 0;

    % Create a cell array to store the segments
    segments = {};
    
    % Create an array to store the results of check_segment_overlap
    segmentOverlaps = {};

    % Loop through each rectangle
    for i = 1:length(x_pos)
        % Calculate the surface area of the current rectangle
        surface_area = x_length(i) * y_length(i);
        
        % Add the current surface area to the total
        total_surface_area = total_surface_area + surface_area;

        % Calculate the perimeter of the current rectangle
        perimeter = 2 * (x_length(i) + y_length(i));

        % Add the current perimeter to the total
        total_perimeter = total_perimeter + perimeter;


        % Define and add the coordinates of the individual segments one by one
        top_segment = [x_pos(i), y_pos(i); x_pos(i) + x_length(i), y_pos(i)]; % Top
        segments{end + 1} = top_segment;
        
        
        right_segment = [x_pos(i) + x_length(i), y_pos(i); x_pos(i) + x_length(i), y_pos(i) + y_length(i)]; % Right
        segments{end + 1} = right_segment;
        

        bottom_segment = [x_pos(i) + x_length(i), y_pos(i) + y_length(i); x_pos(i), y_pos(i) + y_length(i)]; % Bottom
        segments{end + 1} = bottom_segment;
        

        left_segment = [x_pos(i), y_pos(i) + y_length(i); x_pos(i), y_pos(i)]; % Left
        segments{end + 1} = left_segment;
        
    end

    % Loop through the existing segments
    for i = 1:numel(segments)
        segmentOverlaps{i} = compute_segment_total_overlap(segments, segments{i});
        total_perimeter = total_perimeter - segmentOverlaps{i};
    end
    
    % Return the results
    S = total_perimeter * z_length(1) + 2 * total_surface_area;
    V = total_surface_area * z_length(1);
end