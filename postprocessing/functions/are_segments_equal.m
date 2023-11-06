function equal = are_segments_equal(segment1, segment2)
    % Check if the two segments are equal, considering opposite directions
    equal = isequal(sort(segment1), sort(segment2));
end
