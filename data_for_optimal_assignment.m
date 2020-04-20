% PREPARING DATA FOR OPTIMAL ASSIGNMENT

if kk~=0 % distinction between the cases in which there are/aren't trips of size 2
    e = [e1_save; e2_save];
    c_i_j = [e1_save(:,1); e2_save(:,1)];
else
    e = e1_save;
    c_i_j = [e1_save(:,1)];
end

%%% CONSTRAINTS MATRICES %%%

% This matrices has been built so that to be able to put in matrix form the
% constraints and objective function of the optimization problem. They
% exploit information of the edges, telling if requests are linked with
% each other and if vehicles are linked to requests.

% They are needed so that to assign only one trip to requests and to assign
% only one vehicle to one trip

% I_V matrix
[e_r,e_c] = size(e);
for v=1:num_v
    for z=1:e_r
        if v==e(z,2)
            I_V(v,z)=1;
        end
    end
end

% I_T matrix
for t=1:tot
    for z=1:e_r
        if t==e(z,3)
            I_T(t,z) = 1;
        end
    end
end

I_R = RTV_Adj(1:N,N+num_v+1:end);

A_const1 = [I_V zeros(num_v,N)]; % its column size is the size of the variables epsilon_i_j
b_const1 = ones(num_v,1);
A_const2 = [I_R*I_T eye(N)]; % its column size is the size of variables epsilon_i_j + variables X_k
b_const2 = ones(N,1);

%%% COST FUNCTION MATRICES %%%

% c_i_j = [e1_save(:,1); e2_save(:,1)]; % cost vector that associate each epsilon (every edge) with a cost
c_k0 = 1e2*ones(N,1); % large constant to penalize unassigned requests: as many columns as the number of requests
C = [c_i_j; c_k0]; % This matrix has both the cost of the edges and the cost of unassigned requests
eps_init = [Init_guess; X_k_init]; % Initial guess for the variables: epsilon and Xk (look Word file)
X_length = length(eps_init);

intcon = 1:X_length; % This tells the ILP which variables we want to be integers

%%% OPTIMAL ASSIGNMENT COMPUTATION

% The following calls the solver to get X (optimal solution)
X = intlinprog(C,intcon,A_const1,b_const1,A_const2,b_const2,zeros(X_length,1),ones(X_length,1),eps_init);

% The following matrices are built so that to get from the optimal
% assignment X (vector of zeros and ones) the solution in terms of assigned
% edges from the total bunch of edges e

% assigned edges
X_diag = diag(X);
opt_ass = X_diag(1:length(e),1:length(e))*e;
opt_ass = opt_ass(any(opt_ass,2),:);
opt_ass = cast(opt_ass,'single');
% toc

% set of assigned requests
R_OK = [opt_ass*[0 0 0 1 0 0]'; opt_ass*[0 0 0 0 1 0]'];
R_OK = R_OK(any(R_OK,2),:);
% set of assigned vehicles
V_OK = opt_ass*[0 1 0 0 0 0]';

% number of assigned requests and number of assigned vehicles
L_R_OK = length(R_OK);
L_V_OK = length(V_OK);