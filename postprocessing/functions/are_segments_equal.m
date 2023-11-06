function equal = are_segments_equal(segment1, segment2)
    % Check if the two segments are equal, considering opposite directions
    equal = isequal(segment1, segment2) || isequal(segment1, flipud(segment2));
end
