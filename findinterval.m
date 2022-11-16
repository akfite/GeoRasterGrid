function ibin = findinterval(lower, upper, query)
%FINDINTERVAL Identify which limits contain a value.
%
%   Usage:
%
%       INDEX = FINDINTERVAL(LOWER, UPPER, QUERY)
%
%   Inputs:
%
%       LOWER <double vector of length N>
%           - the lower bounding points (such that all(LOWER <= UPPER) is true)
%           - must be sorted in increasing order
%
%       UPPER <double vector of length N>
%           - the upper bound for each pair
%           - must be sorted in increasing order
%
%       QUERY <double matrix>
%           - points to check within the lower/upper bounds
%
%   Outputs:
%
%       INDEX <double matrix (same size as QUERY)>
%           - the index into vector N for the first bounding pair found
%           - zero if no match was found
%
%   Description:
%
%       For each query point Q in QUERY, this function finds the first INDEX that
%       satisfies Q >= LOWER(INDEX) && Q < UPPER(INDEX).  For it to work, the bounding
%       vectors must be sorted in increasing order.  There may be gaps between the
%       bounding pairs, but the ranges of the bounding pairs cannot overlap.
%
%   Example:
%
%       % 0.5 second gap between intervals
%       t1 = 1:10;
%       t2 = 1.5:10.5;
%       query = [-1 1.25 1.75 4.25 15];
%       idx = findinterval(t1, t2, query)

%   Author:     Austin Fite
%   Contact:    akfite@gmail.com

    % flatten bins to a 1-d vector of contiguous edges
    edges = [lower(:)'; upper(:)'-eps(upper(end))];
    edges = edges(:);

    assert(issorted(edges),...
        'MATLAB:findinterval:notsorted',...
        'The bounding points are not sorted or overlap exists.');

    [~,~,ibin] = histcounts(query, edges);

    % map bins back to their original index
    iz = mod(ibin,2) == 0;
    ibin = (ibin + 1) / 2;
    ibin(iz) = 0;
    
end
