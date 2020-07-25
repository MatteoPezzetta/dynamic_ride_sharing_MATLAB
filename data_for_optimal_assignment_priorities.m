% PREPARING DATA FOR OPTIMAL ASSIGNMENT

if kk~=0 % distinction between the cases in which there are/aren't trips of size 2
    e = [e1_save; e2_save];
    c_i_j = [e1_save(:,1); e2_save(:,1)];
else
    e = e1_save;
    c_i_j = [e1_save(:,1)];
end

if kkk~=0
    e = [e; e3_save];
    c_i_j = [c_i_j; e3_save(:,1)];
end

%%%%%%% SOME TRIPS OF SIZE 3 DO NOT HAVE THE COST !!!!! %%%%%%%% WHY ??? %%%%%%%%%%

[e_r,e_c] = size(e); % dimensions of matrix e of feasible edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXAMPLE OF CASE IN WHICH I KNOW WHICH IS THE ALREADY ASSIGNED TASK THAT
% IS ENTERING AGAIN THE PROBLEM. I SHOULD MAKE IT MORE GENERAL SO THAT I
% CAN DO THIS FOR THE ENTIRE SET OF TASKS THAT ARE ONLY ON SIZE 1 TRIPS SO
% THAT CAN BE TAKEN INTO ACCOUNT (actually I should check also for the
% capacity of the vehicle)

% I should get a list of size one trips that can be reassigned (also in
% RV-graph I have to take into account them

% In the following, I am building 'vector1' and 'vector2'.
% ->vector1: used to force to assign at least one of the trips involving the already
% assigned request and the vehicle serving it;
% ->vector2: used to force trips involving the assigned task but other vehicles to NOT be assigned
% and trips involving the assigned vehicle but only other tasks to NOT be assigned.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These operations about vector1 should be modified so that I generate vector1 depending on the priority different form zero ALSO

jj = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size1_exist == 1 % I could put together the if and for with a while a pop the vector size1_assigned
    
    [rows_size1,cols_size1] = size(size1_assigned);
    
    s = 1;
    ee = [e zeros(e_r,1)];
    for j=1:rows_size1 % size1_assigned HAS TO BE GENERATED WITH REQUESTS IDs

        %%% Here I should consider also the tasks that has been assigned in the preassignment phase
        %%% in a way that in the worst case I assign those tasks as in the pre-assignment phase
        vehicle = size1_assigned(j,1);
        request = size1_assigned(j,2);

        for i=1:e_r
            % to force to assign again at least one of the trips linked to v-r couple
            if e(i,2) == vehicle && (e(i,4) == request || e(i,5) == request || e(i,6) == request) % the e(i,5) maybe is not needed because it will never be
                vector1(j,i) = 1;
            else
                vector1(j,i) = 0;
            end
            % to force the trips involing r or v but other vehicles or other
            % requests to not be assigned
             % to add a marker to the set so that to understand if I already took off one edge
            if ((ee(i,2) == vehicle && (ee(i,4) ~= request && ee(i,5) ~= request && ee(i,6) ~= request)) || ((ee(i,4) == request || ee(i,5) == request || ee(i,6) == request) && ee(i,2) ~= vehicle)) && (ee(i,7) == 0)
                vector2(s,:) = zeros(1,e_r);
                vector2(s,i) = 1;
                ee(i,7) = 1;
                s = s+1;
            end
        end
    end
    
    r_vector1 = j;
    r_vector2 = s-1;
end

%%% Loop to forcibly assign the edges containing the priority_ass tasks
if priority_const_exist == 1
    % Here I have to assign at least one of the trips I have assigned in priority assignment
    % I have done the Priority assignment to decide which trips to assign.
    % I could instead have computed the total distance form all vehicles of each requests with priority
    % and assign at least the first n tasks of the list starting form the smallest total distance or from the
    % largest total distance.
    % It is like to create a queue of tasks with priority
    % But the priority assignment is good if I NEED to assign tasks with priority. I should therefore make a list of
    % tasks with priority starting from the highest priority possible, and assign them first
    
    [rows_priority_ass,cols_priority_ass] = size(priority_ass);
    
    for j=1:rows_priority_ass
        
        request = priority_ass(j,2);
        
        for i=1:e_r
            
            if (e(i,4) == request || e(i,5) == request || e(i,6) == request)
                vector3(j,i) = 1;
            else
                vector3(j,i) = 0;
            end
        end
    end
    [r_vector3,c_vector3] = size(vector3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONSTRAINTS MATRICES %%%

% This matrices has been built so that to be able to put in matrix form the
% constraints and objective function of the optimization problem. They
% exploit information of the edges, telling if requests are linked with
% each other and if vehicles are linked to requests.

% They are needed so that to assign only one trip to requests and to assign
% only one vehicle to one trip

% I_V matrix

for v=1:num_v
    for z=1:e_r
        if v==e(z,2)
            I_V(v,z)=1;
        else
            I_V(v,z)=0;
        end
    end
end

% I_T matrix
% maybe if I put it wrong, I have no chances for that trips to be reassigned.
% But the priority constraint are forcing them to be assigned. I should adjust it with the trip numbers
% In practice this error do not lead to errors. This because the fact that is assigned or not is told by the
% additional lines of constraints that we put in A_const2 and b_const2 that force them to be assigned and not "not assigned".
% So this problem is not actually a problem.

for t=1:tot
    for z=1:e_r
        if t==e(z,3)
            I_T(t,z) = 1;
        end
    end
end

I_R = RTV_Adj(1:N,N+num_v+1:end); % part of the RTV_Adj matrix (the symmetric) is not build
                                  % beacuse of timing issues encountered for trips of size 3

% BUILD THE MATRICES OF THE CONSTRAINTS FOR THE ILP

% inequality contraint
A_const1 = [I_V zeros(num_v,N)];
b_const1 = ones(num_v,1);

% A_const2 = [I_R*I_T eye(N)];
% b_const2 = [ones(N,1)];

% equality constraint (two otpions depending if there are size 1 trips from previous algorithm calls)
if size1_exist == 1 % in the case I can try to replay with already assigned trips
    A_const2 = [I_R*I_T eye(N); vector1 zeros(r_vector1,N); vector2 zeros(r_vector2,N)];
    b_const2 = [ones(N,1); ones(r_vector1,1); zeros(r_vector2,1)];
else
    A_const2 = [I_R*I_T eye(N)]; % its column size is the size of variables epsilon_i_j + variables X_k
    b_const2 = ones(N,1);
end

if priority_const_exist == 1
    A_const2 = [A_const2; vector3 zeros(r_vector3,N)];
    b_const2 = [b_const2; ones(r_vector3,1)];
end

%%% COST FUNCTION MATRICES %%%

% c_i_j = [e1_save(:,1); e2_save(:,1)]; % cost vector that associate each epsilon (every edge) with a cost
c_k0 = 1e2*ones(N,1); % large constant to penalize unassigned requests: as many columns as the number of requests
C = [c_i_j; c_k0]; % This matrix has both the cost of the edges and the cost of unassigned requests
C = cast(C,'double');
eps_init = [Init_guess; X_k_init]; % Initial guess for the variables: epsilon and Xk (look Word file)
%eps_init = zeros(1,e_r+N);
X_length = length(eps_init);

intcon = 1:X_length; % This tells the ILP which variables we want to be integers
                     % actually I am omitting this information


% When there are more tasks with high priority than
if rows_priority_ass<length(priority_tasks)
    for i=1:e_r
        if ismember(e(i,5),priority_tasks) || ismember(e(i,6),priority_tasks)
            C(i) = C(i); % It seems that this command is not doing very much
        end
    end
end

%%% OPTIMAL ASSIGNMENT COMPUTATION

% The following calls the solver to get X (optimal solution)

% tic
% X = intlinprog(C,intcon,A_const1,b_const1,A_const2,b_const2,zeros(X_length,1),ones(X_length,1),eps_init);
X = cplexbilp(C,A_const1,b_const1,A_const2,b_const2,eps_init);
% A1 = [1 1; 0 1];
% A2 = [1 0; 1 1];
% b1 = [1 1]';
% b2 = [1 1]';
% costo = [2 5];
% X = cplexbilp(costo,A1,b1,A2,b2,[]);
% toc

% The following matrices are built so that to get from the optimal
% assignment X (vector of zeros and ones) the solution in terms of assigned
% edges from the total bunch of edges e

% assigned edges
X_diag = diag(X);
opt_ass = X_diag(1:e_r,1:e_r)*e; % e_r is the number of rows of e
opt_ass = opt_ass(any(opt_ass,2),:);
opt_ass = cast(opt_ass,'single');

% set of assigned requests
R_OK = [opt_ass*[0 0 0 1 0 0 0]'; opt_ass*[0 0 0 0 1 0 0]'; opt_ass*[0 0 0 0 0 1 0]']; % here I just exptrapolate info on the requests' IDs
R_OK = R_OK(any(R_OK,2),:);
% set of assigned vehicles
V_OK = opt_ass*[0 1 0 0 0 0 0]'; % here I extract info on the vehicles of the assignment

% number of assigned requests and number of assigned vehicles
L_R_OK = length(R_OK);
L_V_OK = length(V_OK);