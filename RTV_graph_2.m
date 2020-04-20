%%% RTV-graph design (construction of feasible trips)

k = 0; % number of trips of size 1 found
kk = 0; % number of trips of size 2 found
out = 0;
tot = 0; % number of total trips found (k+kk)
p = 0; % variable used to count, for each vehicle, the number of trips of size 1

n = 9; % to limit the number of combination and speed up the computation

num_v = num_v; % number of vehicles of the problem

% Initialization of Trip data structure (both size 1 and size 2 together)
T(1).r1 = 0;
T(1).r2 = 0;
T(1).order = 0;
T(1).status = 0;
T(1).size = 0;

for v=1:num_v % for each vehicle do:
    
    % ADD TRIPS OF SIZE ONE
    
    idx = find(Adj(N+v,:)); % The 'N+v-th' row, tells the requests that can be satisfied by vehicle v
    if length(idx)>n % when the size of the problem is big
        idx = datasample(idx,n,'Replace',false); % to pick just 'n' random requests so that to speed up the computation
    end
    
    for i=1:length(idx)
        index = idx(i);
        cost = tt([V(v).x,V(v).y],[R(index).xo,R(index).yo])+tt([R(index).xo,R(index).yo],[R(index).xd,R(index).yd]);
        c_t = cost-R(index).at; % cost in terms of delay on expected arrival
        %%% HERE CHECK FOR PRESENCE OF T
        if ismember(idx(i),[T(:).r1])==1
            indice = find([T(:).r1]==idx(i),1);  %%% find(...,1) beacuse for sure it is in the first position if present
            % update the edge between the existing Trip(indice) and vehicle v
            RTV_Adj(N+v,N+num_v+indice) = 1;
            RTV_Adj(N+num_v+indice,N+v) = 1;
            cost_edge(v,indice) = cost; % we save the travel cost of the edge
            c_t_edge(v,indice) = c_t; % cost in terms of delay on expected arrival
            p = p+1;
            T1_idx(v,p) = index;
        else
            % this is the case when we create a new Trip since it
            % does not exist yet
            k = k+1; % increase number of trips of size 1
            tot = tot+1; % increase total number of trips
            p = p+1;
            index = idx(i);
            % creation of new trip of index tot
            T(tot).r1 = index; % generation of trip from request
            T(tot).r2 = 0;
            T(tot).order = 0;
            T(tot).size = 1;
            cost_edge(v,tot) = cost; % we save the travel cost of the edge
            c_t_edge(v,tot) = c_t; % cost in terms of delay on expected arrival
            T(tot).status = 0;
            T1_idx(v,p) = index;
%             trip1_idx(k) = k; % this is not used maybe -> check if I can throw it
            % Add the edge between the task and the new Trip to wich it belongs;
            % Add the edge between the Trip we have just created and the
            % vehicle v
            RTV_Adj(index,N+num_v+tot) = 1; % e(r,T1)
            RTV_Adj(N+num_v+tot,index) = 1;
            RTV_Adj(N+num_v+tot,N+v) = 1; % e(T1,v)
            RTV_Adj(N+v,N+num_v+tot) = 1;
        end
    end

% ADD TRIPS OF SIZE TWO

    % Trips of size 2 are searched among the trips of size 1 of vehicle v
    % (This because of feasiblity issues (Explained in the Word file)
    for i=1:length(idx) % these are the indeces of trips of size 1 (line 24 of this file)
        for j=1:length(idx)
                  index1 = T1_idx(v,i);
                  index2 = T1_idx(v,j);
                if (index1~=index2) && (Real_Adj(index1,index2)==1) % Explained in the Adj file
                    
                    % check for feasibility of the trip
                    [value,cost,c_t_1,c_t_2] = travel2(V(v),R(index1),R(index2),Kind_of_link(index1,index2));
                    
                    if value==1 % if trip is feasible
                            
                        % check if this Trips of size 2 (the combination of
                        % the two tasks) already exists). If YES, we only
                        % create the edge between the already existing Trip
                        % and the vehicle v
                        if ismember(index1,[T(:).r1])==1
                            indice = find([T(:).r1]==index1);
                            for s=1:length(indice) % to control the second index
                                if T(indice(s)).r2 == index2 % in this case the trip is already present so I only add a 1 to link the trip to the vehicle
                                    RTV_Adj(N+v,N+num_v+indice(s));
                                    RTV_Adj(N+num_v+indice(s),N+v);
                                    cost_edge(v,indice(s)) = cost; % we save the travel cost of the edge
                                    c_t_edge(v,indice(s)) = c_t_1+c_t_2;
                                    out = 1;
                                end
                            end
                        end
                        
                        if out~=1 % if the trip was not already existing
                            kk = kk+1;
                            tot = tot+1;
                            % creation of new trip of size 2
                            T(tot).r1 = index1;
                            T(tot).r2 = index2;
                            T(tot).order = 1;
                            T(tot).status = 0;
                            T(tot).size = 2;
                            cost_edge(v,tot) = cost; % we save the travel cost of the edge
                            c_t_edge(v,tot) = c_t_1+c_t_2;
%                             trip2_idx(kk) = kk;
                            % Add the edges between the single tasks
                            % and the Trip of size 2 we have just created.
                            % Add the edge between this Trip and vehicle v
                            RTV_Adj(index1,N+num_v+tot) = 1;
                            RTV_Adj(N+num_v+tot,index1) = 1;
                            RTV_Adj(index2,N+num_v+tot) = 1;
                            RTV_Adj(N+num_v+tot,index2) = 1;
                            RTV_Adj(N+num_v+tot,N+v) = 1;
                            RTV_Adj(N+v,N+num_v+tot) = 1;
                        end
                        out = 0;
                    end % if of feasibility
                end % end of index1~=index2) && (Adj(index1,index2)==1)
        end % end of 'for j'
    end  % end of 'for i' -- end of 'for' for 'size two trips'
    clear indice;
    p = 0;
end % end of 'for v'

% We can speed it up by giving a certain amount of iteration or maximum
% number of combination so that to not rise the computation time a lot.
% An example could be of limiting the number of feasible trips per vehicle:
% this is done using the 'dataset' function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I should have as a performance parameter, the number of unassigned
% requests with respect to the total. As this parameter becomes smaller,
% the quality would be better
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% As we can decrease the number of combinations with the use of the dataset
% function to pick randomly elements from idx, we should also have less
% trips to use in the optimization problem, so that also the optimization
% problem would benefit from this approach.
% Doing so, we are taking less optimal solution, so we are accepting that
% our solution is not optimal, but given the less computational effort,
% this could be still good. We should find a way to measure the optimality
% of our approach.

% Greedy_assignment