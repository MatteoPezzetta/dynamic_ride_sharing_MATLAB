% This script runs only if both the number of assigned vehicles and the
% number of assigned requests are less than the corresponing total number
% of vehicles and total number of requests

if L_R_OK<N && L_V_OK<num_v
     
    R_KO = 0; % vector of unassigned requests initialized to zero
    V_idle = 0; % vector of idle vehicles (unassigned) initialized to zero

    R_OK = cast(R_OK,'single');
    V_OK = cast(V_OK,'single');

    % construction of a vector of unassigned requests, containing the
    % indeces of those requests
    if L_R_OK<N
        R_KO = setdiff([1:N],R_OK);
        L_R_KO = length(R_KO);
    end

    % construction of a vector of idle vehicles, containing the indeces of
    % those vehicles
    if L_V_OK<num_v
        V_idle = setdiff([1:num_v],V_OK);
        L_V_idle = length(V_idle);
    end
    
    % Now we want to build the matrices to solve the linear program for the
    % assignment of unassigned requests to idle vehicles. As for the ILP,
    % one request cannot be assigned to more than one vehicles, and one
    % vehicle cannot be assigned to more than one request. In this case we
    % assume only trips of size one, so we don't build trips but we
    % consider the set of unassigned requests as a set of trips of size 1
    
    A2 = [];
    R_pre = 1;
    new_v = 0;
    for v=1:L_V_idle
        idx_v = V_idle(v);
        for r=1:L_R_KO
            idx_r = R_KO(r);
            t_cost(1,R_pre) = tt([V(idx_v).x,V(idx_v).y],[R(idx_r).xo,R(idx_r).yo]);
            e_reb(R_pre,:) = [t_cost(1,R_pre),idx_v,idx_r]; % construction of a vector of trips, with cost, vehicle index and request index
            Ar(v,r+new_v*L_R_KO) = 1;
            R_pre = R_pre+1;
        end
        A2 = [A2 eye(L_R_KO,L_R_KO)];
        new_v = new_v+1;
    end
    A1 = Ar;
    b1 = ones(L_V_idle,1);
    b2 = ones(L_R_KO,1);
    A_reb_in = [A1;A2];
    b_reb_in = [b1;b2];
    A_reb = ones(1,L_R_KO*L_V_idle); % A_reb,b_reb for the inequality about rebalancing
    b_reb = min(L_R_KO,L_V_idle);

    % LINEAR PROGRAM for the assignment
    
    % intcon = [1:(L_V_idle*L_R_KO)]';
    num_var = (L_V_idle*L_R_KO);
    % X_reb = intlinprog(C,intcon,A1,b1,A_reb,b_reb,zeros(num_var,1),ones(num_var,1),[]);
    X_reb = linprog(t_cost,A_reb_in,b_reb_in,A_reb,b_reb,zeros(num_var,1),ones(num_var,1));

    Xmat_reb = diag(X_reb);
    opt_reb = Xmat_reb*e_reb;
    opt_reb = opt_reb(any(opt_reb,2),:);
    [r_opt_reb,c_opt_reb] = size(opt_reb);
    
    opt_ass_final = [opt_ass; opt_reb(:,1:2) zeros(r_opt_reb,1) opt_reb(:,3) zeros(r_opt_reb,1)];
else
    opt_ass_final = opt_ass;
end