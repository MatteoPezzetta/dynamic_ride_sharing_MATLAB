%%%

a = 0;
b = 0;

% creation of edges as vectors of vectors (each different row is a
% different edge, containing the edge information on the columns)
% The edges are build only looking at '1' entries of the RTV_graph
for t=1:tot
    for v=1:num_v
        % edges for trips of size 2
        if (RTV_Adj(N+v,N+num_v+t)==1) && (T(t).size==2) && (kk~=0)
            a = a+1;
            e2(a,:) = [cost_edge(v,t) v t T(t).r1 T(t).r2];
%             e2(a,:) = [c_t_edge(v,t) v t T(t).r1 T(t).r2]; % when cost is measured as delay on the arrival
        end
        % edges for trips of size 1
        if (RTV_Adj(N+v,N+num_v+t)==1) && (T(t).size==1)
            b = b+1;
            e1(b,:) = [cost_edge(v,t) v t T(t).r1];
%             e1(b,:) = [c_t_edge(v,t) v t T(t).r1]; % when cost is measured as delay on the arrival
        end
    end
end

e1 = sortrows(e1); % sorte e1 so that to assign edges in increasing cost 
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

i_g = length(e1); % to count indeces for the initial guess. The vector of all the edges has
                  % first the trips of size 1 and then the trips of size 2.
                  % The greedy assignment, though, assign first trips of
                  % size 2, so that we need to skip the indeces of trips of
                  % size 1. Then, i_g = 1 just before assign trips of size
                  % 1

% We start by assigning trips of size 2, if present. We assign trips
% composed of tasks that has not been already assigned and to vehicles that
% has not been already assigned. Once there are no more trips of size 2
% that we can assign, we start assign trips of size 1.

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
            R_OK(i) = sel(1,4);
            R_OK(i+1) = sel(1,5);
            % We create a vecotr containing the indeces of the assigned
            % vehicles
            V_OK(j) = sel(1,2);
            i = i+2;
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
        R_OK(i) = sel(1,4);
        % We add vehicles to the vector of assigned vehicles
        V_OK(j) = sel(1,2);
        i = i+1;
        Greedy_ass(j,:) = [sel(1,2),sel(1,3)];
        flag = 1; % this flag says that the edge has been assigned, so we put 1 in Init guess
        j = j+1;
        T(sel(1,3)).status = 1; % the trip is assigned
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
for w=1:N
    if ismember(w,R_OK)==0
        X_k_init(q,1) = 1;
    end
    q = q+1;
end