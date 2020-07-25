%%% RTV-graph design (construction of feasible trips)
% sss = 1; % used just as a check
% ffflag = 1; % used just as a check
% exist = 0; % used just as a check
% not_exist = 0; % used just as a check
% exist2 = 0; % used just as a check
% travel3_calls = 0; % used just as a check
number_value = 0;

% task_ID = [];
% all_trips = [];
k = 0; % number of trips of size 1 found
kk = 0; % number of trips of size 2 found
kkk = 0; % number of trips of size 3 found
out = 0;
tot = 0; % number of total trips found (k+kk)
p = 0; % variable used to count, for each vehicle, the number of trips of size 1
pp = 0; % variable used to count, for each vehicle, the number of trips of size 2
ppp = 0;

T1_idx = zeros(num_v,1);
T2_idx = zeros(num_v,1);
T3_idx = zeros(num_v,1);

n = 7; % to limit the number of combination and speed up the computation
n_2 = 6;
num_v = num_v; % number of vehicles of the problem

% Initialization of Trip data structure (both size 1 and size 2 together)
T(1).rPU = [0 0];
T(1).order = 0;
T(1).status = 0;
T(1).size = 0;

T1_index = 0;
T2_index = 0;
T3_index = 0;

for v=1:num_v % for each vehicle do:
    
    % ADD TRIPS OF SIZE ONE
    
    idx_original = find(Adj(N+v,:)); % The 'N+v-th' row, tells the requests that can be satisfied by vehicle v
    
    if length(idx_original)>n % when the size of the problem is big
%          idx_sorted = sort([TripLength(N+v,idx)' idx']);
%          idx = idx_sorted(1:n,2);
         idx = datasample(idx_original,n,'Replace',false); % to pick just 'n' random requests so that to speed up the computation
         % here I should add to idx the tasks that have priorities
         %
         % verify for elements = 2 in Prioriy_Mat
         % if yes, take the index of the R for which it happens
         % then check if that index exists on the idx vector
         % if not, add it
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if ismember(2,priority_mat(N+v,:))
             r_indice = find([priority_mat(N+v,:)] == 2,1);
             if ismember(r_indice,idx)==0
                 idx = [idx,r_indice]; % check if it is a column or row vector
             end
         end    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        idx = idx_original;
    end
    
    for i=1:length(idx)
        index = idx(i); % this is the actual index from the list of requests
        cost = tt([V(v).x,V(v).y],[R(index).xo,R(index).yo])+tt([R(index).xo,R(index).yo],[R(index).xd,R(index).yd]);
        c_t = cost-R(index).at; % cost in terms of delay on expected arrival
        %%% HERE CHECK FOR PRESENCE OF T
        found_T1 = 0;
        for t=1:length(T1_index)
            t_index = T1_index(t);
            if t_index~=0
%                 sono_qui = 1
                if all([R(index).ID 0]==[T(t_index).rPU(1:2)])
                    % update the edge between the existing Trip(indice) and vehicle v
                    RTV_Adj(N+v,N+num_v+t_index) = 1; %%%% INDICE SHOULD REMAIN BECAUSE IT IS THE POSITION IN THE VECTOR T(:)
%                     RTV_Adj(N+num_v+t_index,N+v) = 1;
                    cost_edge(v,t_index) = cost; % we save the travel cost of the edge
                    c_t_edge(v,t_index) = c_t; % cost in terms of delay on expected arrival
                    found_T1 = 1;
                end
            end
        end
        
        if found_T1 == 1
%             sono_qui = 1
            p = p+1;
            T1_idx(v,p) = index; %%%% THIS SHOULD REMAIN BECAUSE DEALS WITH POSITION ON THE VECTOR OF TRIPS/REQUESTS NOT THE ACTUAL REQUEST ID
        end 
        
        if found_T1 == 0
            % this is the case when we create a new Trip since it
            % does not exist yet
            k = k+1; % increase number of trips of size 1
            tot = tot+1; % increase total number of trips
            p = p+1;
            index = idx(i);
            % creation of new trip of index tot
            T(tot).rPU = [R(index).ID 0];
            T(tot).index1 = index;
            T(tot).index2 = 0;
            T(tot).order = 0;
            T(tot).possible_order = 0;
            T(tot).size = 1;
            T(tot).pass = R(index).pass;
            cost_edge(v,tot) = cost; % we save the travel cost of the edge
            c_t_edge(v,tot) = c_t; % cost in terms of delay on expected arrival
            T(tot).status = 0;
            T1_idx(v,p) = index;
            T1_index(k) = tot;
%             trip1_idx(k) = k; % this is not used maybe -> check if I can throw it
            % Add the edge between the task and the new Trip to wich it belongs;
            % Add the edge between the Trip we have just created and the
            % vehicle v
            RTV_Adj(index,N+num_v+tot) = 1; % e(r,T1)
            %RTV_Adj(N+num_v+tot,index) = 1;
            %RTV_Adj(N+num_v+tot,N+v) = 1; % e(T1,v)
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
                    [value,cost,c_t_1,c_t_2,pass] = travel2(V(v),R(index1),R(index2),Kind_of_link(index1,index2));
                    
%                     % WHAT IS THIS????? I THINK ONLY A CHECK FOR DEBUGGING
%                     if value == 0
%                         verification(sss,:) = [R(index1).ID, R(index2).ID, V(v).ID];
%                         sss = sss+1;
%                         ffflag = 0;
%                     end
                    if value==1 % if trip is feasible
                        out = 0;
                        % check if this Trip of size 2 (the combination of
                        % the two tasks) already exists). If YES, we only
                        % create the edge between the already existing Trip
                        % and the vehicle v
%                         T2_idx_indeces = find(T2_idx(v,:));
%                         T2_idx = T2_idx(T2_idx_indeces);
                        for t=1:length(T2_index)
                            t_index = T2_index(t);
                            if t_index ~= 0
                                if all([R(index1).ID R(index2).ID] == [T(t_index).rPU(1:2)])
%                                     exist = exist + 1;
                                    RTV_Adj(N+v,N+num_v+t_index) = 1;
                                    RTV_Adj(N+num_v+t_index,N+v) = 1;
                                    
                                    cost_edge(v,t_index) = cost; % we save the travel cost of the edge
                                    c_t_edge(v,t_index) = c_t_1 + c_t_2;
                                    out = 1;
                                    break;
                                end
                            end
                        end
                        
                        if out == 1
                            pp = pp+1;
                            T2_idx(v,pp) = t_index;
                        end
                        
                        if out~=1 % if the trip was not already existing
%                             not_exist = not_exist + 1;
                            kk = kk+1;
                            tot = tot+1;
                            % creation of new trip of size 2
                            T(tot).rPU = [R(index1).ID R(index2).ID];
                            T(tot).index1 = index1;
                            T(tot).index2 = index2;
                            T(tot).order = Kind_of_link(index1,index2);
                            T(tot).possible_order = Kind_of_link2(index1,index2);
                            T(tot).status = 0;
                            T(tot).size = 2;
                            T(tot).pass = pass;
                            
                            % to save which trips of size 2 the vehicle has
                            pp = pp+1;
                            T2_idx(v,pp) = tot; % because it was not existing
                            T2_index(kk) = tot;
                            
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
    
% [ADD TRIPS OF SIZE 3]
    
    % In the C++ code maybe I do not need this step because the vertices have different
    % dimensions and to not have zeros as ending elements
    T1_indeces = find(T1_idx(v,:));
    T1_indeces = T1_idx(v,T1_indeces);
    
    T_indeces = find(T2_idx(v,:));
    
    if length(T_indeces)>n_2 % when the size of the problem is big
%          idx_sorted = sort([TripLength(N+v,idx)' idx']);
%          idx = idx_sorted(1:n,2);
         T_indeces = datasample(T_indeces,n_2,'Replace',false);
    end
    
    if length(T_indeces)>0
        T_indeces = T2_idx(v,T_indeces); % These are the indeces of trips of size 2 relative to that vehicle v
    
        for i=1:length(T1_indeces)
            for j=1:length(T_indeces)
                id_r3 = T1_indeces(i);
                id_T2 = T_indeces(j);
                id_r1 = T(id_T2).index1; % I take the index of the Graph, not the actual R.ID
                id_r2 = T(id_T2).index2;
                % Maybe I can avoid this if. Worst case I check an unreal case
                if (id_r1~=id_r2 && id_r2~=id_r3 && id_r1~=id_r3)
                    T3_trial = [id_r1,id_r2,id_r3];
                    
                    flag1 = 1;
                    flag2 = 1;
                    
                    % Here we check that all the T/r_i are trips of size2 of the current vehicle
                    for y=1:length(T_indeces)
                        % The following condition is enough, because using the RTV_Adj, doing so I know if a generic size 2 trip
                        % for that vehicle is containing the two tasks, no matter the order.
                        if (RTV_Adj(N+num_v+T_indeces(y),id_r1)==1 && RTV_Adj(N+num_v+T_indeces(y),id_r3)==1)
                            flag1 = flag1+1;
                        end
                        
                        if (RTV_Adj(N+num_v+T_indeces(y),id_r2)==1 && RTV_Adj(N+num_v+T_indeces(y),id_r3)==1)
                            flag2 = flag2+1;
                        end
                    end
                    
                    % Here I should verify the order in which the tasks are combined. Then I will have two cases
                    % depending on 1'-2' or 2'-1';
                    
                    % Entering the following condition means that the 3 tasks both belong to some
                    % trips of size2 of the current vehicle
                    if (flag1 > 1) && (flag2 > 1) % This means that all tasks belong to trips of size 2 of the current vehicle
                        % In the following conditions the oder of the tasks is really important
                        % since Real_Adj conserves that actual order of the tasks
                        
                        pu_order1 = 0;
                        pu_order2 = 0;
                        pu_order3 = 0;
                        
                        if Real_Adj(id_r1,id_r3)==1 && Real_Adj(id_r2,id_r3)==1
                            pu_order1 = 1; % order 1,2,3 is okay
                        end
                        if Real_Adj(id_r1,id_r3)==1 && Real_Adj(id_r3,id_r2)==1
                            pu_order2 = 1; % order 1,3,2 is okay
                        end
                        if Real_Adj(id_r3,id_r1)==1 && Real_Adj(id_r3,id_r2)==1
                            pu_order3 = 1; % order 3,1,2 is okay
                        end
                        % Here I should go on with verifying order1, order2, order3 for dropoffs
                        % The conditions change depending on the order of dropoff of the trip T2(r1,r2)
                        % When the subdivision between order1 and order2 trips will be implemented, I'll need to
                        % look at T.order data instead of Kind_of_Link. I could already do so...
                        
                        do_order1 = 0;
                        do_order2 = 0;
                        do_order3 = 0;
                        
                        if Kind_of_link(id_r1,id_r2)==1 % Here I should look into T().order
                            if (Kind_of_link2(id_r1,id_r3)~=2 || Kind_of_link2(id_r3,id_r1)~=1)...
                                    && (Kind_of_link2(id_r2,id_r3)~=2 || Kind_of_link2(id_r3,id_r2)~=1)
                                do_order1 = 1;
                            end
                            if (Kind_of_link2(id_r1,id_r3)~=2 || Kind_of_link2(id_r3,id_r1)~=1)...
                                    && (Kind_of_link2(id_r3,id_r2)~=2 || Kind_of_link2(id_r2,id_r3)~=1)
                                do_order2 = 1;
                            end
                            if (Kind_of_link2(id_r3,id_r1)~=2 || Kind_of_link2(id_r1,id_r3)~=1)...
                                    && (Kind_of_link2(id_r3,id_r2)~=2 || Kind_of_link2(id_r2,id_r3)~=1)
                                do_order3 = 1;
                            end
                            
                        elseif Kind_of_link(id_r1,id_r2)==2 % Here I should look into T().order
                            if (Kind_of_link2(id_r2,id_r3)~=2 || Kind_of_link2(id_r3,id_r2)~=1)...
                                    && (Kind_of_link2(id_r1,id_r3)~=2 || Kind_of_link2(id_r3,id_r1)~=1)
                                do_order1 = 1;
                            end
                            if (Kind_of_link2(id_r2,id_r3)~=2 || Kind_of_link2(id_r3,id_r2)~=1)...
                                    && (Kind_of_link2(id_r3,id_r1)~=2 || Kind_of_link2(id_r1,id_r3)~=1)
                                do_order2 = 1;
                            end
                            if (Kind_of_link2(id_r3,id_r2)~=2 || Kind_of_link2(id_r2,id_r3)~=1)...
                                    && (Kind_of_link2(id_r3,id_r1)~=2 || Kind_of_link2(id_r1,id_r3)~=1)
                                do_order3 = 1;
                            end
                        end 
                        
                        % Then here call function travel3() (if I am not already including the code above into travel3()
                        
                        task_id = [id_r1,id_r2,id_r3];
                        pu_order = [pu_order1, pu_order2, pu_order3];
                        do_order = [do_order1, do_order2, do_order3];
                        KoL = Kind_of_link(id_r1,id_r2);
                        
%                         travel3_calls = travel3_calls + 1;
                        
                        % Call travel3() to check for feasibility of the possible combinations for the size3 trip
                        % I do not want it to be called many times
                        [value, trip_cost, trip, pass] = travel3(task_id,V(v),pu_order,do_order,KoL,R);
                        
                        % Generation of the trips of size 3
                        if value==1
                            
%                             t = index_min;
                            out_3 = 0;
                            
                            % Now I have to search if the trips already exist, and if not, generate a new one.
                            % When the trip is already present I have to link it the thasks and to the vehicle
                            for x=1:length(T3_index)
                                t_index = T3_index(x);
                                if t_index ~= 0
                                    % This condition tells if the size3 trip alreadu exists or not
                                    if all(trip == [T(t_index).rPU T(t_index).rDO])
%                                         exist2 = exist2+1;
                                        out_3 = 1;
                                        RTV_Adj(N+v,N+num_v+t_index) = 1;
                                        RTV_Adj(N+num_v+t_index,N+v) = 1;

                                        cost_edge(v,t_index) = trip_cost;
                                    end
                                end
                                
                                if out_3 == 1
                                    ppp = ppp+1;
                                    T3_idx(v,ppp) = t_index;
                                end
                                
                                if out_3~=1 % In this case the size3 trips does not exist already and we have to add it
                                    kkk = kkk+1;
                                    tot = tot+1;
                                    T3_index(kkk) = tot; % This set of indeces is helpful in tracking the whole set of size 3 trips
                                    
                                    ppp = ppp+1;
                                    T3_idx(v,ppp) = tot;
                                    
                                    cost_edge(v,tot) = trip_cost;
                                    
                                    % Generation of new trip of size 3
                                    T(tot).rPU = [trip(1) trip(2) trip(3)];
                                    T(tot).rDO = [trip(4) trip(5) trip(6)];
                                    
                                    T(tot).size = 3;
                                    T(tot).pass = pass;
                                    
                                    % Update the RTV_Adj matrix
                                    
                                    RTV_Adj(trip(1),N+num_v+tot) = 1;
%                                     RTV_Adj(N+num_v+tot,trip(t,1)) = 1;
                                    RTV_Adj(trip(2),N+num_v+tot) = 1;
%                                     RTV_Adj(N+num_v+tot,trip(t,2)) = 1;
                                    RTV_Adj(trip(3),N+num_v+tot) = 1;
%                                     RTV_Adj(N+num_v+tot,trip(t,3)) = 1;
                                    
                                    RTV_Adj(N+v,N+num_v+tot) = 1;
%                                     RTV_Adj(N+num_v+tot,N+v) = 1;
                                    
                                end
                            end
                        end % for t=1:length(value)
                    end
                end
            end
        end
    end
    clear indice;
    p = 0;
    pp = 0;
    ppp = 0;
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