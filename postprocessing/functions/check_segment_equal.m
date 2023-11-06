function equalExists = check_segment_equal(segments, oneSegment)
    % Initialize a flag to check if an equal segment exists
    count = 0;

    % Loop through the existing segments
    for i = 1:numel(segments)
        if are_segments_equal(segments{i}, oneSegment)
            count = count + 1;
        end
    end

    equalExists = count > 1;
end