% This function compute feasibility of a trip of size 1.
% It returns '1' if the trip is feasible, '0' otherwise

% Feasibility is evaluated on the travel time between the position of the
% vehicle and the pick-up time of the task, depending on the accepted
% delay Delta

function [value,cost_v_pick,priority] = travel(V,R,priority_ass)

    global Delta;
    global DeltaInf;
    v_pos = [V.x;V.y];
    R_orgn = [R.xo;R.yo];
    R_dest = [R.xd;R.yd];
    t1 = tt(v_pos,R_orgn);
    t2 = tt(R_orgn,R_dest);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the case when the couple v-r has been assigned in the pre-assignment
    % for high priority tasks
    if length(priority_ass)>0
        if ismember(V.ID,[priority_ass(:,1)])
            v_indice = find([priority_ass(:,1)] == V.ID, 1);
            if priority_ass(v_indice,2) == R.ID
                value = 1;
                cost_v_pick = t1 + t2;
                priority = 2;
                return % check if this return is right or not
            end
        end
    end
    
    if R.priority == 1
        Delta = DeltaInf;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This is the standard case when the v-r couple is not on the priority task list
    if (t1<=(R.st+Delta) && R.pass <= V.c)
        value = 1;
        cost_v_pick = t1 + t1;
        priority = 1;
    else
        value = 0;
        cost_v_pick = 0;
        priority = 0;
    end
end

% the distance required for the task, from its starting point until the
% dropoff, must be coherent with the arrival time with respect to the
% starting time fo the task itself