% This function compute feasibility of a trip of size 2.
% It returns '1' if the trip is feasible, '0' otherwise

% Feasibility is evaluated on the travel time between different points.
% Tasks can be combined with order of type 1 (o_i->o_j->d_i->d_j)
% or order of type 2 (o_i->o_j->d_j-d_i). The two cases are treated with a
% 'switch'. The link type is entered as argument (c) of the function.
% Other arguments are: vehicle position(V), task 1 pick-up and delivery
% points (R1) and task 2 pick-up and delivery points (R2).
% The function returns the feasibility of the trip (value = '0'/'1'), the
% cost of the trip (in time), the delay on arrival time of task 1 (c_t_1)
% and on task 2 (c_t_2).

% Feasibility is ensured (for case 1) when:
% -> The arrival time of R1 is satisfied, given maximum delay Delta;
% -> The arrival time of R2 is satisfied, given maximum delay Delta;
% -> The starting time of R2 is satisfied, given maximum delay Delta.

% The satisfaction of arrival time of R1 is not taken into account. This
% becuase trips of size 2 require, BY DEFINITION, that the starting point of both the tasks
% that compose the trip can be reached by the specific vehicle taken into
% account.
% Explanation: if, absurdly, task 2 cannot be directly reached by the
% vehicle taken into account, for sure, it would not be reachable in the
% situation where the vehicle has to pick-up task 1 first. This beacuse the
% travelled distance beween vehicle V and starting point of R2 is always shorter than 
% (the distance between vehicle V and starting point of R2) + (the distance
% between starting point of R1 and starting point of R2).

% Feasibility for case 2 is similar:
% -> The arrival time of R2 is satisfied, given maximum delay Delta;
% -> The arrival time of R1 is satisfied, given maximum delay Delta;
% -> The starting time of R2 is satisfied, given maximum delay Delta.

% Notice that case 1 and case 2 have the order of drop-off switched: This
% is the only difference between the two cases.


function [value,cost,c_t_1,c_t_2,pass] = travel2(V,R1,R2,KoL)

    global Delta;
    global DeltaInf;
    
    pass = R1.pass+R2.pass;
    
    if pass <= V.c
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % In this way I am saying that the second task, if has a high priority, can accept
        % to be satisfied also with a big delay. The important thing is that I increase the
        % chance for it to be served
        if R1.priority == 0 && R2.priority == 1
            Delta1 = Delta;
            Delta2 = DeltaInf; % a great constant
        elseif R1.priority == 1 && R2.priority == 1
            Delta1 = DeltaInf;
            Delta2 = DeltaInf;
        elseif R1.priority == 1 && R2.priority == 0 % this is the case when first task has priority
            Delta1 = DeltaInf;
            Delta2 = Delta;
        else
            Delta1 = Delta;
            Delta2 = Delta;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % maybe the conditions are not totally right because in RV graph I stated that
        % the second task cannot be the one that is already assigned/picked up
%         if (R1.status == 1 || R1.status == 2) && (R2.status == 2 || R2.status == 1) % so maybe this condition is not needed at all
%             
%             value = 0;
%             cost = 0;
%             c_t_1 = 0;
%             c_t_2 = 0;
            
        if (R1.status == 0 || R1.status == 1)

            v_pos = [V.x;V.y];
            R1_orgn = [R1.xo;R1.yo];
            R1_dest = [R1.xd;R1.yd];
            R2_orgn = [R2.xo;R2.yo];
            R2_dest = [R2.xd;R2.yd];

            %%%% HERE I SHOULD ENTER A CONDITION TO VERIFY IF THE FIRST TASK HAS ALREADY BEEN PICKED UP.
            %%%% IF YES, I SHOULD CHANGE THE CRITERIA WITH WHICH I COMPUTE THE TRAVEL DISTANCES.
            %%%% WHAT COST TO TAKE INTO ACCOUNT? ALSO THE ALREADY TRAVELES DISTANCE BY THE VEHICLE OR ONYL THE DISTANCE TRAVELED
            %%%% FORM THE EXACT MOMENT IN TIME???

            % maybe I already have the information on the trip_length
            switch KoL
                case 1
                    if (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R1_dest)<=R1.at+Delta1)...
                            && (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R1_dest)+tt(R1_dest,R2_dest)<=R2.at+Delta2)...
                            && (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)<=R2.st+Delta2)
                        value = 1;
                        cost1 = tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R1_dest); % cost on reaching destination R1
                        cost2 = cost1+tt(R1_dest,R2_dest); % cost of reaching destionation of R2
                        cost = cost2;
                        c_t_1 = cost1-R1.at; % delay on reaching destination 1
                        c_t_2 = cost2-R2.at; % delay on reaching destionation 2
                    else
                        value = 0;
                        cost = 0;
                        c_t_1 = 0;
                        c_t_2 = 0;
                    end
                case 2
                    if (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R2_dest)<=R2.at+Delta2)...
                            && (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R2_dest)+tt(R2_dest,R1_dest)<=R1.at+Delta1)...
                            && (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)<=R2.st+Delta2)
                        value = 1;
                        cost2 = tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R2_dest); % time taken to reach destination of R2
                        cost1 = cost2+tt(R2_dest,R1_dest); % time taken to reach destination of 1
                        cost = cost1;
                        c_t_2 = cost2-R2.at; % delay on destination 2
                        c_t_1 = cost1-R1.at; % delaty on destination 1
                    else
                        value = 0;
                        cost = 0;
                        c_t_2 = 0;
                        c_t_1 = 0;
                    end     
            end
%%%%%%    
% when one of the two tasks is already picked up and the other is not assigned yet
        elseif R1.status == 2 % the if regardng the status (when: (R1.status == 2 && R2.status == 0) || (R1.status == 0 %% R2.status == 2)
            
            v_pos = [V.x;V.y];
            R1_orgn = [R1.xo;R1.yo];
            R1_dest = [R1.xd;R1.yd];
            R2_orgn = [R2.xo;R2.yo];
            R2_dest = [R2.xd;R2.yd];
            
            switch KoL
                case 1
                    if (tt(v_pos,R2_orgn)+tt(R2_orgn,R1_dest)<=R1.at+Delta1)...
                            && (tt(v_pos,R2_orgn)+tt(R2_orgn,R1_dest)+tt(R1_dest,R2_dest)<=R2.at+Delta2)...
                            && (tt(v_pos,R2_orgn)<=R2.st+Delta2)
                        value = 1;
                        cost1 = tt(v_pos,R2_orgn)+tt(R2_orgn,R1_dest); % cost on reaching destination R1 %%%% LOOK HOW I WANT TO ACTUALLY COMPUTE THAT COST
                        cost2 = cost1+tt(R1_dest,R2_dest); % cost of reaching destionation of R2
                        cost = cost2;
                        c_t_1 = cost1-R1.at; % delay on reaching destination 1
                        c_t_2 = cost2-R2.at; % delay on reaching destionation 2
                    else
                        value = 0;
                        cost = 0;
                        c_t_1 = 0;
                        c_t_2 = 0;
                    end
                case 2
                    if (tt(v_pos,R2_orgn)+tt(R2_orgn,R2_dest)<=R2.at+Delta2)...
                            && (tt(v_pos,R2_orgn)+tt(R2_orgn,R2_dest)+tt(R2_dest,R1_dest)<=R1.at+Delta1)...
                            && (tt(v_pos,R2_orgn)<=R2.st+Delta2)
                        value = 1;
                        cost2 = tt(v_pos,R2_orgn)+tt(R2_orgn,R2_dest); % time taken to reach destination of R2
                        cost1 = cost2+tt(R2_dest,R1_dest); % time taken to reach destination of 1
                        cost = cost1;
                        c_t_2 = cost2-R2.at; % delay on destination 2
                        c_t_1 = cost1-R1.at; % delaty on destination 1
                    else
                        value = 0;
                        cost = 0;
                        c_t_2 = 0;
                        c_t_1 = 0;
                    end     
            end
        end %(end of the if, else if and else that takes into account already assigned/picked up tasks
%%%%%%
    else
        value = 0;
        cost = 0;
        c_t_2 = 0;
        c_t_1 = 0;
    end