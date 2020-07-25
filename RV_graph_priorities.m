%%% RV-graph design

% clear all
% close all
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I SHOULD TAKE INTO ACCOUNT WHEN num_v = 0 or R = 0; So that the algortihm will not run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
N = N; % length(R_ass) + length(R_priority) + length(new_tasks)
num_v = num_v; % num_v = new_vehicles + length(V_ass) % this is the number of actually available vehicles that enter the algorithm call
max_cap = max_cap;

% NOW HERE I NEED TO ADD TO THE SET OF VEHICLES AND THE SET OF REQUESTS THE
% OLD REQUESTS AND VEHICLES I HAD. Maybe I need to add a request and
% vehicle ID

% shareability graph between routes

global Delta; % This Delta is used in the standard case. Both tasks of normal priority
Delta = 5;
global DeltaInf; % This Delta will be used for the second task when it has high priority
DeltaInf = 10;

Adj = zeros(N+num_v,N+num_v);
Real_Adj = zeros(N,N);
TripLength = zeros(N,N);
priority_mat = zeros(N+num_v,N+num_v);
Kind_of_link = zeros(N,N);
Kind_of_link2 = zeros(N,N);
% TripLength2 = zeros(N,N);

n_1 = N;
total_indeces = [1:N];
for i=1:N
    
    %total_idx = datasample(total_indeces,n_1,'Replace',false);
    total_idx = total_indeces;
    for j_total_idx=1:n_1
  
        j = total_idx(j_total_idx);
        
        if (R(i).pass + R(j).pass) <= max_cap
            orgn_i = [R(i).xo,R(i).yo];
            dest_i = [R(i).xd,R(i).yd];
            orgn_j = [R(j).xo,R(j).yo];
            dest_j = [R(j).xd,R(j).yd];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This if works when the second task as higher priority and I try to put it as second task
            if R(i).priority == 0 && R(j).priority == 1 % In this case the second task has high priority so I facilitate it to be assigned as second task
                Delta1 = Delta;
                Delta2 = DeltaInf; % a big constant. It is a magic number. I should modify it as a global variable
            elseif R(i).priority == 1 && R(j).priority == 1
                Delta1 = DeltaInf;
                Delta2 = DeltaInf;
            elseif R(i).priority == 1 && R(j).priority == 0
                Delta1 = DeltaInf;
                Delta2 = Delta;
            else
                Delta1 = Delta;
                Delta2 = Delta;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            d1 = tt(orgn_i,orgn_j);
            d2 = tt(orgn_j,dest_i);
            d3 = tt(dest_i,dest_j);
            
            %%% for oi->oj->di->dj %%% (first kind of link between tasks)
            b1 = [R(i).st+Delta1;
                -R(i).st;
                R(j).st+Delta2-d1;
                -R(j).st;
                R(i).at+Delta1-d1-d2;
                R(j).at+Delta2-d1-d2-d3];
            a1 = [1 -1 1 -1 1 1]';

            ubd = min([b1(1),b1(3),b1(5),b1(6)]);
            lbd = max(b1(2),b1(4));

            % update the Adj matrix
            if (ubd>lbd) && (i~=j) && (ubd>=0) && (R(j).status == 0) %%%% it tells that the second task is not already assigned/picked up
                %create Adj between trip i and j
                Adj(i,j) = 1;
                Adj(j,i) = 1;
                Real_Adj(i,j) = 1;
                Adj1(i,j) = 1;
                clear ubd lbd; % to then run also the other kind of combination
                TripLength(i,j) = d1+d2+d3;
                Kind_of_link(i,j) = 1;
                Kind_of_link2(i,j) = 1;
            end

            d2 = tt(orgn_j,dest_j);
            d3 = tt(dest_j,dest_i);
            
            %%% for oi->oj->dj->di %%% (second kind of link between tasks)
            b2 = [R(i).st+Delta1;
                -R(i).st;
                R(j).st+Delta2-d1;
                -R(j).st;
                R(j).at+Delta2-d1-d2;
                R(i).at+Delta1-d1-d2-d3];
            a2 = [1 -1 1 -1 1 1]';

            ubd = min([b2(1),b2(3),b2(5),b2(6)]);
            lbd = max(b2(2),b2(4));

            % update the Adj matrix
            if (ubd>lbd) && (i~=j) && (ubd>=0) && (R(j).status == 0) %%%% it tells that the second task is not already assigned/picked up
                %create Adj between trip i and j
                Adj(i,j) = 1;
                Adj(j,i) = 1;
                Real_Adj(i,j) = 1;
                Adj2(i,j) = 1;
%                 clear ubd lbd; % to then run also the other combinations
                T_length = d1+d2+d3;
                if Kind_of_link(i,j) == 1
                    Kind_of_link2(i,j) = 3;
                else
                    Kind_of_link2(i,j) = 2;
                end
                % check if link of type 2 is more convenient than type 1 (in
                % travel time). If yes, the two tasks are combined with link of
                % type 2), otherwise the link remains 1.
                if (TripLength(i,j) == 0 || T_length<TripLength(i,j))
                    Kind_of_link(i,j) = 2;
                    TripLength(i,j) = T_length;
                end
            end
        end
    end
end

% ASSIGNMENT OF VEHICLES TO REQUESTS (single requests)

% assignment of vehicles to requests depending on 'return' of travel(V,R)
for i=1:num_v
    for j=1:N
        % If the vehicle and the request can be combined, the edge between
        % the node of the vehicle and the node of the request is generated
        [Adj(i+N,j),TripLength(N+i,j),priority_mat(N+i,j)] = travel(V(i),R(j),priority_ass); % it returns 1 if they can be assigned
        Adj(j,i+N) = Adj(i+N,j);
    end
end
% toc
% % % % % % % % HERE ENTERS A VECTOR THAT TELLS EACH ALREADY ASSIGNED VEHICLE TO WHICH REQUEST IS LINKED

% % % % % % % % If the RV-graph has not considered the couple, we add it
% % % % % % % % the 2D vector is called for example: priority_ass = [ vehicle , task ];
% % % % % % % for i=1:priority_ass_length
% % % % % % %     r_indice = find([R(:).ID] == priority_ass(i,2),1);
% % % % % % %     % if this is zero, then we make it 1. So we are forcing the fact that these vehicles have to serve these requests
% % % % % % %     % We could simply add the edges lated, but doing so, we can then also build size 2 trips.
% % % % % % %     % MAYBE I CAN DO ALL THIS IS FUNCTION TRAVEL 1: IF I SEE THAT THE COUPLE BELONG TO THE priority_ass VECTOR,
% % % % % % %     % THEN I FOCRIBLY ASSIGN IT
% % % % % % %     if Adj(priorit_ass(i,1),r_indice) == 0 % this is the position  relative to the vehicle and to the task
% % % % % % %         