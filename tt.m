% This function compute the travel time between two points (called 'origin' and 'destination').
% The velocity parameter is assumed arbitrarily: for real implementation,
% we should compute how much time the vehicle spend to travel 1 meter
% distance.

% distance d in this case is computed as the Manhattan distance ( following
% taxicab geometry) between two points

function travel_time = tt(orgn,dest)
    vel_k = 0.08;
    d = norm((orgn-dest),1); % manhattan distance between the two points
    travel_time = vel_k*d;
end