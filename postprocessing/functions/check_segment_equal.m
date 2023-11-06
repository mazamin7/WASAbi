function equalExists = check_segment_equal(segments, newSegment)
    % Initialize a flag to check if an equal segment exists
    equalExists = false;

    % Loop through the existing segments
    for i = 1:numel(segments)
        if are_segments_equal(segments{i}, newSegment)
            equalExists = true;
            break;  % An equal segment was found, no need to check further
        end
    end
end