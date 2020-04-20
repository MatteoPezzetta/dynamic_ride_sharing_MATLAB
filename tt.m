% This function compute the travel time between two points (called 'origin' and 'destination').
% The velocity parameter is assumed arbitrarily: for real implementation,
% we should compute how much time the vehicle spend to travel 1 meter
% distance.

% distance d in this case is computed as the length of the difference
% vector between the two locations. We should use the Manhattan distance
% for a more realistic implementation.

function travel_time = tt(orgn,dest)
    vel_k = 0.08;
    d = norm(orgn-dest);
    travel_time = vel_k*d;
end