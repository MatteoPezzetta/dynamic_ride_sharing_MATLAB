% This function compute feasibility of a trip of size 1.
% It returns '1' if the trip is feasible, '0' otherwise

% Feasibility is evaluated on the travel time between the position of the
% vehicle and the pick-up time of the task, depending on the accepted
% delay Delta

function [value] = travel(V,R)
    global Delta;
    v_pos = [V.x;V.y];
    R_orgn = [R.xo;R.yo];
    if tt(v_pos,R_orgn)<=(R.st+Delta);
        value = 1;
    else
        value = 0;
    end
end

% the distance required for the task, from its starting point until the
% dropoff, must be coherent with the arrival time with respect to the
% starting time fo the task itself