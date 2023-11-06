function overlapLength = compute_segment_overlap(segment1, segment2)
    % Check if there is an overlap between two segments
    % Calculate the overlap length if an overlap exists

    % Sort the segments to ensure consistent ordering
    segment1 = sort(segment1);
    segment2 = sort(segment2);

    % Check for overlap based on the orientation and coordinate match of the segments
    if is_horizontal(segment1) && is_horizontal(segment2) && segment1(1, 2) == segment2(1, 2)
        % Both segments are horizontal and have the same y-coordinate
        overlapLength = compute_horizontal_overlap(segment1, segment2);
    elseif is_vertical(segment1) && is_vertical(segment2) && segment1(1, 1) == segment2(1, 1)
        % Both segments are vertical and have the same x-coordinate
        overlapLength = compute_vertical_overlap(segment1, segment2);
    else
        % Segments have different orientations or coordinates, so no overlap
        overlapLength = 0;
    end
end

function isHorizontal = is_horizontal(segment)
    % Check if a segment is horizontal (parallel to the x-axis)
    isHorizontal = segment(1, 2) == segment(2, 2);
end

function isVertical = is_vertical(segment)
    % Check if a segment is vertical (parallel to the y-axis)
    isVertical = segment(1, 1) == segment(2, 1);
end

function overlapLength = compute_horizontal_overlap(segment1, segment2)
    % Calculate the overlap length between two horizontal segments
    overlapLength = max(0, min(segment1(2, 1), segment2(2, 1)) - max(segment1(1, 1), segment2(1, 1)));
end

function overlapLength = compute_vertical_overlap(segment1, segment2)
    % Calculate the overlap length between two vertical segments
    overlapLength = max(0, min(segment1(2, 2), segment2(2, 2)) - max(segment1(1, 2), segment2(1, 2)));
end
