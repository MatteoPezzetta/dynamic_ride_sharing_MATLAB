function [value, trip_cost, trip, pass] = travel3(task_id,V,pu_order,do_order,KoL,R)
% body
    global Delta;
%     cost_min = 5e4;
%     index_min = 0;
    
    id_r1 = task_id(1);
    id_r2 = task_id(2);
    id_r3 = task_id(3);
    
    num_pass = R(id_r1).pass+R(id_r2).pass+R(id_r3).pass;
    if num_pass <= V.c
        
        number_of_trips = 0;
        for i_t3=1:length(pu_order)
            for j_t3=1:length(do_order)
                
                number_of_trips = number_of_trips + 1;
                
                if pu_order(i_t3)~=0 && do_order(j_t3)~=0 % when both are equal to 1, I have a feasible combination
                     % It could be that more size3 trips are generated from that function
                    % order of tasks for pick-up
                    
                    switch i_t3
                        case 1
                            task_1 = id_r1;
                            task_2 = id_r2;
                            task_3 = id_r3;
                        case 2
                            task_1 = id_r1;
                            task_2 = id_r3;
                            task_3 = id_r2;
                        case 3
                            task_1 = id_r3;
                            task_2 = id_r1;
                            task_3 = id_r2;
                    end % end of switch i 

                    if KoL == 1
                        % order of tasks for drop-off
                        switch j_t3
                            case 1
                                task_4 = id_r1;
                                task_5 = id_r2;
                                task_6 = id_r3;
                            case 2
                                task_4 = id_r1;
                                task_5 = id_r3;
                                task_6 = id_r2;
                            case 3
                                task_4 = id_r3;
                                task_5 = id_r1;
                                task_6 = id_r2;
                        end
                    elseif KoL == 2 % In this case 1-2 are connected in a different way and so also 1-3 and 2-3
                        switch j_t3
                            case 1
                                task_4 = id_r2;
                                task_5 = id_r1;
                                task_6 = id_r3;
                            case 2
                                task_4 = id_r2;
                                task_5 = id_r3;
                                task_6 = id_r1;
                            case 3
                                task_4 = id_r3;
                                task_5 = id_r2;
                                task_6 = id_r1;
                        end
                    end
                   
                    % The capacity is a problem to deal with with size3 trips

                    v_pos = [V.x,V.y];
                    R1_orgn = [R(task_1).xo; R(task_1).yo];
                    R4_dest = [R(task_4).xd; R(task_4).yd];
                    R2_orgn = [R(task_2).xo; R(task_2).yo];
                    R5_dest = [R(task_5).xd; R(task_5).yd];
                    R3_orgn = [R(task_3).xo; R(task_3).yo];
                    R6_dest = [R(task_6).xd; R(task_6).yd];
                    
                    % Verify if I can shorten this condition in some way
                    if (tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)<=R(task_2).st+Delta...
                            && tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R3_orgn)<=R(task_3).st+Delta...
                            && tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R3_orgn)+tt(R3_orgn,R4_dest)<=R(task_4).at+Delta...
                            && tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R3_orgn)+tt(R3_orgn,R4_dest)+tt(R4_dest,R5_dest)<=R(task_5).at+Delta...
                            && tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R3_orgn)+tt(R3_orgn,R4_dest)+tt(R4_dest,R5_dest)+tt(R5_dest,R6_dest)<=R(task_6).at+Delta)
%                         index_min = number_of_trips;
                        value = 1;
                        trip_cost = tt(v_pos,R1_orgn)+tt(R1_orgn,R2_orgn)+tt(R2_orgn,R3_orgn)+tt(R3_orgn,R4_dest)+tt(R4_dest,R5_dest)+tt(R5_dest,R6_dest); % Here I should use the assignment inside the condition
                        trip = [R(task_1).ID, R(task_2).ID, R(task_3).ID, R(task_4).ID, R(task_5).ID, R(task_6).ID]; % trips should be made of the IDs
                        pass = num_pass;
                        return
                    else
                        value = 0; % vectorial datum or not?
                        trip_cost = 0;
                        trip = [0, 0, 0, 0, 0, 0];
                        pass = 0;
                    end
                else
                    value = 0; % vectorial datum or not?
                    trip_cost = 0;
                    trip = [0, 0, 0, 0, 0, 0];
                    pass = 0;
                end
            end
        end
    else % in the case the number of passengers exceed the maximum capacity of the vehicle
        value = 0;
        trip_cost = 0;
        trip = [0, 0, 0, 0, 0, 0];
        pass = 0;
    end
end % end of function