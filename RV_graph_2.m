%%% RV-graph design

clear all
close all
clc

N = 50;
num_v = 30;
max_cap = 3;

% trips initialization as structures
% - we should then add the number of passengers
for i=1:floor(N/4)
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo + 50*rand(1,1);
    R(i).yd = R(i).yo + 50*rand(1,1);
    R(i).st = 1;
    R(i).at = 5;
    R(i).pass = randi(2,1); % number of passengers: a random number form 1 to 2
end
for i=floor(N/4)+1:floor(N/2)
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo - 50*rand(1,1);
    R(i).yd = R(i).yo - 50*rand(1,1);
    R(i).st = 1;
    R(i).at = 5;
    R(i).pass = randi(2,1);
end
for i=floor(N/2)+1:floor(3*N/4)
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo + 50*rand(1,1);
    R(i).yd = R(i).yo - 50*rand(1,1);
    R(i).st = 1;
    R(i).at = 5;
    R(i).pass = randi(2,1);
end
for i=floor(3*N/4)+1:N
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo - 50*rand(1,1);
    R(i).yd = R(i).yo + 50*rand(1,1);
    R(i).st = 1;
    R(i).at = 5;
    R(i).pass = randi(2,1);
end

% shareability graph between routes
global Delta;
Delta = 4;
Adj = zeros(N+num_v,N+num_v);
Real_Adj = zeros(N,N);
TripLength1 = zeros(N,N);
TripLength2 = zeros(N,N);

% tic
for i=1:N
    for j=1:N
  
        if (R(i).pass + R(j).pass) <= max_cap
            orgn_i = [R(i).xo,R(i).yo];
            dest_i = [R(i).xd,R(i).yd];
            orgn_j = [R(j).xo,R(j).yo];
            dest_j = [R(j).xd,R(j).yd];

            %%% for oi->oj->di->dj %%% (first kind of link between tasks)
            b1 = [R(i).st+Delta;
                -R(i).st;
                R(j).st+Delta-tt(orgn_i,orgn_j);
                -R(j).st;
                R(i).at+Delta-tt(orgn_i,orgn_j)-tt(orgn_j,dest_i);
                R(j).at+Delta-tt(orgn_i,orgn_j)-tt(orgn_j,dest_i)-tt(dest_i,dest_j)];
            a1 = [1 -1 1 -1 1 1]';

            ubd = min([b1(1),b1(3),b1(5),b1(6)]);
            lbd = max(b1(2),b1(4));

            % update the Adj matrix
            if (ubd>lbd) && (i~=j) && (ubd>=0) % '(i~=j)' excludes the case of the task combined with itself
                %create Adj between trip i and j
                Adj(i,j) = 1;
                Adj(j,i) = 1;
                Real_Adj(i,j) = 1;
                Adj1(i,j) = 1;
                clear ubd lbd; % to then run also the other kind of combination
                TripLength1(i,j) = tt(orgn_i,orgn_j)+tt(orgn_j,dest_i)+tt(dest_i,dest_j);
                Kind_of_link(i,j) = 1;
            end

            %%% for oi->oj->dj->di %%% (second kind of link between tasks)
            b2 = [R(i).st+Delta;
                -R(i).st;
                R(j).st+Delta-tt(orgn_i,orgn_j);
                -R(j).st;
                R(j).at+Delta-tt(orgn_i,orgn_j)-tt(orgn_j,dest_j);
                R(i).at+Delta-tt(orgn_i,orgn_j)-tt(orgn_j,dest_j)-tt(dest_j,dest_i)];
            a2 = [1 -1 1 -1 1 1]';

            ubd = min([b2(1),b2(3),b2(5),b2(6)]);
            lbd = max(b2(2),b2(4));

            % update the Adj matrix
            if (ubd>lbd) && (i~=j) && (ubd>=0)
                %create Adj between trip i and j
                Adj(i,j) = 1;
                Adj(j,i) = 1;
                Real_Adj(i,j) = 1;
                Adj2(i,j) = 1;
                clear ubd lbd; % to then run also the other combinations
                TripLength2(i,j) = tt(orgn_i,orgn_j)+tt(orgn_j,dest_j)+tt(dest_j,dest_i);
                % check if link of type 2 is more convenient than type 1 (in
                % travel time). If yes, the two tasks are combined with link of
                % type 2), otherwise the link remains 1.
                if (TripLength1(i,j) == 0 || TripLength2(i,j)<TripLength1(i,j))
                    Kind_of_link(i,j) = 2;
                end
            end
        end
    end
end

% ASSIGNMENT OF VEHICLES TO REQUESTS (single requests)

% initialization of vehicles

for i=1:floor(num_v/4);
    V(i).x = 70*rand(1,1);
    V(i).y = 20*rand(1,1);
    V(i).c = max_cap;
end
for i=floor(num_v/4)+1:floor(num_v/2)
    V(i).x = 80+20*rand(1,1);
    V(i).y = 80+20*rand(1,1);
    V(i).c = max_cap;
end
for i=floor(num_v/2)+1:floor(3*num_v/4)
    V(i).x = 20+20*rand(1,1);
    V(i).y = 80+20*rand(1,1);
    V(i).c = max_cap;
end
for i=floor(3*num_v/4)+1:num_v
    V(i).x = 50+20*rand(1,1);
    V(i).y = 50+100*rand(1,1);
    V(i).c = max_cap;
end

% assignment of vehicles to requests depending on 'return' of travel(V,R)
for i=1:num_v
    for j=1:N
        % If the vehicle and the request can be combined, the edge between
        % the node of the vehicle and the node of the request is generated
        Adj(i+N,j) = travel(V(i),R(j)); % it returns 1 if they can be assigned
        Adj(j,i+N) = Adj(i+N,j);
    end
end