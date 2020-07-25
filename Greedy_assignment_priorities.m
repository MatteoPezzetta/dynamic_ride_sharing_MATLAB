%%%

a = 0;
b = 0;
c = 0;

% creation of edges as vectors of vectors (each different row is a
% different edge, containing the edge information on the columns)
% The edges are build only looking at '1' entries of the RTV_graph

%%%%% I SHOULD ADD TO e1 THE TRIPS OF SIZE 1 THAT I MUST DO BECAUSE ALREADY ASSIGNED.
%%%%% BECAUSE THEY MAY HAVE BEEN SKIPPED WHEN CHOOSING TRIPS RANDOMLY
for t=1:tot
    for v=1:num_v
        % edges for trips of size 2
        if (RTV_Adj(N+v,N+num_v+t)==1) && (T(t).size==3) && (kkk~=0)
            a = a+1;
            e3(a,:) = [cost_edge(v,t) V(v).ID t T(t).rPU(1,1) T(t).rPU(1,2) T(t).rPU(1,3) T(t).pass]; %%%% SHOULD BE V(v).ID
        end
        if (RTV_Adj(N+v,N+num_v+t)==1) && (T(t).size==2) && (kk~=0)
            b = b+1;
            e2(b,:) = [cost_edge(v,t) V(v).ID t T(t).rPU(1,1) T(t).rPU(1,2) 0 T(t).pass]; %%%% SHOULD BE V(v).ID
%             e2(a,:) = [c_t_edge(v,t) v t T(t).r1 T(t).r2]; % when cost is measured as delay on the arrival
        end
        % edges for trips of size 1
        if (RTV_Adj(N+v,N+num_v+t)==1) && (T(t).size==1) && (k~=0)
            c = c+1;
            e1(c,:) = [cost_edge(v,t) V(v).ID t T(t).rPU(1,1) 0 0 T(t).pass]; %%%% SHOULD BE V(V).ID
%             e1(b,:) = [c_t_edge(v,t) v t T(t).r1]; % when cost is measured as delay on the arrival
        end
    end
end

% I need to add ALL the size 1 trips coming from the previou step
% because they MUST BE ASSIGNED.
% This is needed because in the RTV graph formulation, I pick random trips of size 1
% to check for cliques.
% Another good method should be to add in queue to the randoms, (during the creation of T1)
% the v-r that are already assigned. Doing so, maybe the complexity becomes too high.

% Here I am forcibly modifying the set e1 so that to add trips (of size 1) already assigned at the previous call
[r_e1,c_e1] = size(e1);

if size1_exist == 1;
    [rows_R_ass,cols_R_ass] = size(size1_assigned);
    for i=1:rows_R_ass
        flag = 0;
        for j=1:r_e1
            if (size1_assigned(i,1) == e1(j,2)) && (size1_assigned(i,2) == e1(j,4))
                flag = 1;
            end
        end
        if flag == 0
            % Forcibly add the edges to the entire set of edges
            e1 = [size1_assigned(i,3) size1_assigned(i,1) 0 size1_assigned(i,2) 0 0 size1_assigned(i,4); e1];
        end
    end     
end

% maybe here I should do the same thing for size2 trips to reconsider them in the problem

e1 = sortrows(e1); % sort e1 so that to assign edges in increasing cost 
e1_save = e1; % save the edges of size 1 in a new variable because we are
              % goind to destroy the e1 vector in the greedy assignmen
              % so that to exclude already assigned edges

% return
i = 1;
j = 1;
f = 0; % number of vehicles that serve trips of size 1
g = 0; % number of vehicles that serve trips of size 2

R_OK = 0;
V_OK = 0;

i_g = b+c; % to count indeces for the initial guess. The vector of all the edges has
                  % first the trips of size 1 and then the trips of size 2.
                  % The greedy assignment, though, assign first trips of
                  % size 2, so that we need to skip the indeces of trips of
                  % size 1. Then, i_g = 1 just before assign trips of size
                  % 1

if kkk~=0 % we enter here only when trips of size 2 exist
    e3 = sortrows(e3); % sort e2 so that to assign edges in increasing cost
    e3_save = e3; % save the edges of size 1 in a new variable because we are
                  % goind to destroy the e1 vector in the greedy assignmen
                  % so that to exclude already assigned edges
    
    while (length(e3)>0) % as long as we have finished to check trips of size 2

        i_g = i_g+1;
        % pop e2(T,v)
        sel = e3(1,:);
        e3(1,:) = [];

        % assignment of trips if they do not already belong to the
        % assigned routes/assigned vehicles
        if (ismember(sel(1,4),R_OK)==0 && ismember(sel(1,5),R_OK)==0 && ismember(sel(1,6),R_OK)==0 && ismember(sel(1,2),V_OK)==0)
            % We create a vector containing the indeces of the assigned
            % tasks. Trips of size 2 involve 2 requests at a time
            R_OK(i) = sel(1,4); %%%% HERE I AM SAVING THE REQUEST ID, NOT THE INDEX IN THE R VECTOR
            R_OK(i+1) = sel(1,5);
            R_OK(i+2) = sel(1,6);
            % We create a vecotr containing the indeces of the assigned
            % vehicles
            V_OK(j) = sel(1,2);
            i = i+3;
            Greedy_ass(j,:) = [sel(1,2),sel(1,3)]; % Here we create a 2D vector containing the vehicle and the assigned trip
            flag = 1; % this flag says that the edge has been assigned, so we put 1 in Init guess
            j = j+1;
            T(sel(1,3)).status = 1; % the trip is assigned (the status of the trip is updated
            g = g+1; % number of vehicles that serve trips of size 2;
        end

        % In this step we put 1 to assigned edges and 0 to unassigned edges
        if flag==1
            Init_guess(i_g,1) = 1;
            flag = 0;
        else
            Init_guess(i_g,1) = 0;
        end
    end
end

% We start by assigning trips of size 2, if present. We assign trips
% composed of tasks that has not been already assigned and to vehicles that
% has not been already assigned. Once there are no more trips of size 2
% that we can assign, we start assign trips of size 1.
%%% Can this concept hold even in the case of reconsideration of previous assignments???

i_g = c; % to place the elements in the right position in the init_assignment vector

if kk~=0 % we enter here only when trips of size 2 exist
    e2 = sortrows(e2); % sort e2 so that to assign edges in increasing cost
    e2_save = e2; % save the edges of size 1 in a new variable because we are
                  % goind to destroy the e1 vector in the greedy assignmen
                  % so that to exclude already assigned edges
    
    while (length(e2)>0) % as long as we have finished to check trips of size 2

        i_g = i_g+1;
        % pop e2(T,v)
        sel = e2(1,:);
        e2(1,:) = [];

        % assignment of trips if they do not already belong to the
        % assigned routes/assigned vehicles
        if (ismember(sel(1,4),R_OK)==0 && ismember(sel(1,5),R_OK)==0 && ismember(sel(1,2),V_OK)==0)
            % We create a vector containing the indeces of the assigned
            % tasks. Trips of size 2 involve 2 requests at a time
            R_OK(i) = sel(1,4); %%%% HERE I AM SAVING THE REQUEST ID, NOT THE INDEX IN THE R VECTOR
            R_OK(i+1) = sel(1,5);
            % We create a vecotr containing the indeces of the assigned
            % vehicles
            V_OK(j) = sel(1,2);
            i = i+2;
            Greedy_ass(j,:) = [sel(1,2),sel(1,3)]; % Here we create a 2D vector containing the vehicle and the assigned trip
            flag = 1; % this flag says that the edge has been assigned, so we put 1 in Init guess
            j = j+1;
%             T(sel(1,3)).status = 1; % the trip is assigned (the status of the trip is updated
            g = g+1; % number of vehicles that serve trips of size 2;
        end

        % In this step we put 1 to assigned edges and 0 to unassigned edges
        if flag==1
            Init_guess(i_g,1) = 1;
            flag = 0;
        else
            Init_guess(i_g,1) = 0;
        end
    end
end

i_g = 0;

while (length(e1)>0) % as long as we have finished to check trips of size 1
    
    i_g = i_g+1;
    % pop e2(T,v)
    sel = e1(1,:);
    e1(1,:) = [];
    
    % assignment of trips if they do not already belong to the
    % assigned routes/assigned vehicles
    if (ismember(sel(1,4),R_OK)==0 && ismember(sel(1,2),V_OK)==0)
        % We add tasks to the vector of assigned tasks
        R_OK(i) = sel(1,4); %%% HERE I AM SAVING THE R.ID
        % We add vehicles to the vector of assigned vehicles
        V_OK(j) = sel(1,2); %%%% HERE I AM SAVING THE V.ID
        i = i+1;
        % I do not need this datum so I can also cancel it
        Greedy_ass(j,:) = [sel(1,2),sel(1,3)]; % NOT GOOD BECAUSE THE TRIP INDEX IS ZERO(BEACUSE OF THE REALANCING PHASE)
        flag = 1; % this flag says that the edge has been assigned, so we put 1 in Init guess
        j = j+1;
%         T(sel(1,3)).status = 1; % the trip is assigned
        f = f+1; % number of vehcles that serve trips of size 1;
    end
    
    % In this step we put 1 to assigned edges and 0 to unassigned edges
    if flag==1
        Init_guess(i_g,1) = 1;
        flag = 0;
    else
        Init_guess(i_g,1) = 0;
    end
end

% Here we create a vector of size as the total number of requests. It is
% '1' when the corresponding task is not assigned. Assigned tasks can be
% found in vector R_OK
q = 1;
X_k_init = zeros(N,1);
for w=1:N %%%% IT WILL NOT WORK ANYMORE. ismember(R(w).ID,R_OK)
    if ismember(R(w).ID,R_OK)==0 %%%% w IS JUST A NUMBER FROM 1 TO N, NOT THE R.ID, WHILE R_OK IS CONTAINING R.IDs
        X_k_init(q,1) = 1;
    end
    q = q+1;
end