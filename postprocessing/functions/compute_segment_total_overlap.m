function totalOverlap = compute_segment_total_overlap(segments, oneSegment)
    % Initialize the total overlap length
    totalOverlap = 0;

    % Loop through the existing segments
    for i = 1:numel(segments)
        overlapLength = compute_segment_overlap(segments{i}, oneSegment);
        totalOverlap = totalOverlap + overlapLength;
    end

    % Subtract the own segment length as it's part of the total overlap
    ownSegmentLength = compute_segment_overlap(oneSegment, oneSegment);
    totalOverlap = totalOverlap - ownSegmentLength;
end
