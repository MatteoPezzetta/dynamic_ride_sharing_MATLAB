% build vector of tasks that need priority

% deisgn a linprog to assign the best maximum number of vehicles to them

% I can assign a maximum number of tasks equal to the number of available vehicle.
% These tasks for sure will be served, the others will be given a good chance of
% being assigned in the RV-RTV-ILP process

% size1_exist should be put at 1 also in the case a preassignment is called

% in this case update R = [R; R_priority]

% If there are tasks that need priority and if there are new available vehicles
    
if priority_exist == 1 && new_vehicles > 0
    
    A2_priority = [];
    R_pre = 1;
    new_v = 0;
    
    r_pa = length(priority_tasks); % dimensions 'of priority_ass' vector
    priority_ass = [0,0];
    
    for v=1:new_vehicles
        for r=1:r_pa
            % the vetor of new_vehicles should be already available
            t_cost_priority(1,R_pre) = tt([New_V(v).x,New_V(v).y],[R_priority(r).xo,R_priority(r).yo]); % here I should also be able to bypass the internal checks of travel()
            e_pa(R_pre,:) = [New_V(v).ID R_priority(r).ID];
            A_priority(v,r+new_v*r_pa) = 1;
            R_pre = R_pre+1;
        end
        A2_priority = [A2_priority eye(r_pa,r_pa)];
        new_v = new_v+1;
    end
    
    A1 = A_priority;
    b1 = ones(new_vehicles,1);
    A2 = A2_priority;
    b2 = ones(r_pa,1);
    A_priority_ineq = [A1;A2];
    b_priority_ineq = [b1;b2];
    
    A_pa = ones(1,r_pa*new_vehicles);
    b_pa = min(r_pa,new_vehicles);
    
    num_var = (r_pa*new_vehicles);
    
    X_pa = cplexlp(t_cost_priority,A_priority_ineq,b_priority_ineq,A_pa,b_pa,zeros(num_var,1),ones(num_var,1));
    
    Xmat_pa = diag(X_pa);
    
    % This will be the vector that tells which vehicle is linked to which task
    % that was asking for priority
    priority_ass = Xmat_pa*e_pa;
    priority_ass = priority_ass(any(priority_ass,2),:);
    
    priority_const_exist = 1; % This variable is needed to notify 'Data_for_optimal_assignment' that we need an additional set of constraints
    
    % Here I can build the new set of V and R
    
    clear A1 b1 A2 b2
else
    % There are not prioritized assignments. It could be because of no tasks asking for priority
    % or for the lack of new free vehicles.
    priority_const_exist = 0;
    rows_priority_ass = 0;
    priority_ass = [];
end

%%% HERE I NEED TO UPDATE R IN THE CASE WE HAVE SOME TASKS