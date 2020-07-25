% Create vector of tasks of size 1 already assigned that will enter again the problem.
% This vector is needed so that to keep track of what tasks are assigned and what
% vehicles are on their way to the tasks. Once a task is picked up, I should update
% the status of the task, saying that it is not only assigned but also picked up.
% In this way, the function travel2() will consider the fact that the vehicle is already
% half way to the end of one task, and it would deviate to get the other task.
% This will take more time to satisfy the first task, and therefore the possibility of combination
% cannot be very high.

% I should also take into account of creating a set of tasks that are gone so that to have its history
% it will be used to track the tasks during the simulation in dynamic environment.
% R_fixed_assignment tells which are the tasks belonging to trips of size 2 that therefore are no more assignable

j = 0; % to build a vector of assigned tasks in trips of size 1: therefore they can enter again the problem;
f_a = 0; % to build a vector of assigned tasks in a trip of size 2 and size 3;
f_a_2 = 0;
v_a = 0; % to build the vector of assigned vehicles with 2 tasks;

[r_opt_ass,c_opt_ass] = size(opt_ass_final);
for i=1:r_opt_ass % There could also be the case in which I don't have any assignment
    if opt_ass_final(i,5) == 0 && opt_ass_final(i,6) == 0 % it means tasks of size 1
        
        j = j+1;
        size1_assigned(j,:) = [opt_ass_final(i,2) opt_ass_final(i,4) opt_ass_final(i,1) opt_ass_final(i,6)]; % on second column I have the IDs of tasks
        r_indice = find([R(:).ID] == opt_ass_final(i,4),1);
        R(r_indice).status = 1; % update status of the task that is assigned to 1
        R_ass(j) = R(r_indice); % start creating the new set of tasks. The others will be put in queue
        
        v_indice = find([V(:).ID] == opt_ass_final(i,2),1);
        V(v_indice).status = 1;
        V_ass(j) = V(v_indice);
        
        % HERE I HAVE TO UPDATE THE VECTOR OF ALL_V. V IS ONLY THE VECTOR OF VEHICLES CONSIDERED AT THE ALG. CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_v_indice = find([ALL_V(:).ID] == opt_ass_final(i,2),1);
        ALL_V(all_v_indice).status = 1;
        
        size1_exist = 1; % There is at least a size1 assignment from the previous call
                         % ACTUALLY I SHOULD ATTACH IT TO SIZE1 ASSIGNMENTS OF EARLIER CALLS.
                         % WHAT I AM DOING RIGHT NOW IS TO CONSIDER ONLY THE CALL JUST BEFORE THE NEXT ONE
        
    elseif opt_ass_final(i,5)~=0 && opt_ass_final(i,6)==0 % this to keep track of assigned routes that we do not account anymore in the optimization problem beacause they are full
        f_a_2 = f_a_2 + 1;
        R_fixed_assignment(f_a_2,:) = [opt_ass_final(i,2) opt_ass_final(i,4) opt_ass_final(i,5) 0]; % it tells the vehicle and tasks of the route
        r1_indice = find([R(:).ID] == opt_ass_final(i,4),1);
        r2_indice = find([R(:).ID] == opt_ass_final(i,5),1);
        R(r1_indice).status = 1; % update the status of the request to 1 (assigned)
        R(r2_indice).status = 1; % update the status of the request to 1 (assigned)
        f_a = f_a + 1;
        R_fa(f_a) = R(r1_indice);
        f_a = f_a + 1;
        R_fa(f_a) = R(r2_indice);
        
        % update the vector of assigned vehicles with already 2 tasks
        v_a = v_a + 1;
        v_indice = find([V(:).ID] == opt_ass_final(i,2),1);
        V(v_indice).status = 1; % update the status of the vehicle to 1 (assigned)
        V_fa(v_a) = V(v_indice); % This vehicles are the one serving trips of size 2. They need to end their tasks
                                 % before becoming available again
        % HERE I HAVE TO UPDATE THE VECTOR OF ALL_V. V IS ONLY THE VECTOR OF VEHICLES CONSIDERED AT THE ALG. CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_v_indice = find([ALL_V(:).ID] == opt_ass_final(i,2),1);
        ALL_V(all_v_indice).status = 1;
    elseif opt_ass_final(i,6)~=0
        f_a_2 = f_a_2 + 1;
        R_fixed_assignment(f_a_2,:) = [opt_ass_final(i,2) opt_ass_final(i,4) opt_ass_final(i,5) opt_ass_final(i,6)];
        r1_indice = find([R(:).ID] == opt_ass_final(i,4),1);
        r2_indice = find([R(:).ID] == opt_ass_final(i,5),1);
        r3_indice = find([R(:).ID] == opt_ass_final(i,6),1);
        R(r1_indice).status = 1; % update the status of the request to 1 (assigned)
        R(r2_indice).status = 1; % update the status of the request to 1 (assigned)
        R(r3_indice).status = 1; % update the status of the request to 1 (assigned)
        f_a = f_a + 1;
        R_fa(f_a) = R(r1_indice);
        f_a = f_a + 1;
        R_fa(f_a) = R(r2_indice);
        f_a = f_a + 1;
        R_fa(f_a) = R(r3_indice);
        
        % update the vector of assigned vehicles with already 2 tasks
        v_a = v_a + 1;
        v_indice = find([V(:).ID] == opt_ass_final(i,2),1);
        V(v_indice).status = 1; % update the status of the vehicle to 1 (assigned)
        V_fa(v_a) = V(v_indice); % This vehicles are the one serving trips of size 2. They need to end their tasks
                                 % before becoming available again
        % HERE I HAVE TO UPDATE THE VECTOR OF ALL_V. V IS ONLY THE VECTOR OF VEHICLES CONSIDERED AT THE ALG. CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_v_indice = find([ALL_V(:).ID] == opt_ass_final(i,2),1);
        ALL_V(all_v_indice).status = 1;
    end
end

% Generate a vector of unassigned tasks -> they will have priority at the next call
all_assigned_tasks = [opt_ass_final(:,4); opt_ass_final(:,5); opt_ass_final(:,6)]; % there are some zeros in this vector
priority_tasks = [];
j = 0;
for i=1:N
    if ismember(R(i).ID,all_assigned_tasks)==0 % When a task of the set is not assigned at this call, it enters a priority list
        j = j+1;
        priority_tasks = [priority_tasks; R(i).ID]; % Enters the priority list
        R(i).priority = 1; % Task structure is updated to 'priority'
        R_priority(j) = R(i); % these are the tasks used in 'Prioriy_assignment_priorities'
        flag = 1;
    end
end

% if I have at least a task in priority_ass, then I need to notify the next call
% that some tasks are coming from the previous
if flag == 1
%     size1_exist = 1; % This maybe should be updated when actually assigning vehicles and tasks. So in the file 'Priority_assignment'
    priority_exist = 1;
%     R_PRIORITY = [R_PRIORITY R_priority];
%     V_PRIORITY = [V_PRIORITY V_priority];
else % In the case NO prioritized assignment has been done
    R_priority = []; % or at least remains the original one
    V_priority = []; % or at least remains the original one
    priority_exist = 0;
end

% In a real case I must always maintain a R_ass and R_priority that I update at every call

old_N = N;

% clear old tasks

clear R;

% clear old variables except the one needed at next call

clearvars -EXCEPT size1_assigned size1_exist priority_exist R_ass R_priority priority_tasks V_ass ALL_V old_N num_v L_R_OK

% maybe I should not clear vehicles since they will remain the same

% check for available vehicles (done in Generate_V_and_R that comes before Priority assignment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Here I should receive new tasks if they exist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%% IMPORTANT %%%
%%%%%%%%%%%%%%%%%

% WHEN A TASK IS COMPLETED, IT MUST BE DESTROYED. SO THAT IT CANNOT ENTER THE PROBLEM ANYMORE
% --> TASKS BELONGING TO TRIPS OF SIZE1, WILL ALWAYS REMAIN PRESENT UNTIL COMPLETION OR UNTIL THEY ARE
% ASSIGNED TO TRIPS OF SIZE 2
% --> TASKS THAT ARE NOT ASSIGNED WILL ALWAYS ENTER AGAIN THE PROBLEM AND WILL NOT BE LOST UNTIL THEY ARE USSIGNED