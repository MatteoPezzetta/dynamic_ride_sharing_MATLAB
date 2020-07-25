% script to run all the assignment problem.

% This scrit runs in the correct order all the needed scripts to obtain the
% optimal assignment.
% tic-toc does not account for the generation of the plots in evaluationg
% the execution time.

% IF YOU WANT TO SEE THE PLOTS OF ALL THE STEPS (RV_graph, RTV-graph, map
% with requests and vehicles positions) JUST COMMENT THE 'return'

% PERFORMANCE: with 100 tasks and 30 vehicles I got about 1.6 seconds to
% obtain the optimal assignment in my laptop. Implemented in C++ (or outside matlab in general)
% maybe it will take even less time

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This part is of interest only when I need to generate new tasks or generate vehicles

N = 4;
num_v = 4;
max_cap = 5;
old_N = 0;
% trips initialization as structures
% - we should then add the number of passengers
for i=1:floor(N/4)
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo + 50*rand(1,1);
    R(i).yd = R(i).yo + 50*rand(1,1);
    R(i).st = 1;
    R(i).at = R(i).st+tt([R(i).xo,R(i).yo],[R(i).xd,R(i).yd]);
    R(i).pass = randi(2,1); % number of passengers: a random number form 1 to 2
    R(i).ID = old_N + i; %%%% request ID
    R(i).status = 0;
    R(i).priority = 0;
end
for i=floor(N/4)+1:floor(N/2)
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo - 50*rand(1,1);
    R(i).yd = R(i).yo - 50*rand(1,1);
    R(i).st = 1;
    R(i).at = R(i).st+tt([R(i).xo,R(i).yo],[R(i).xd,R(i).yd]);
    R(i).pass = randi(2,1);
    R(i).ID = old_N + i;
    R(i).status = 0;
    R(i).priority = 0;
end
for i=floor(N/2)+1:floor(3*N/4)
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo + 50*rand(1,1);
    R(i).yd = R(i).yo - 50*rand(1,1);
    R(i).st = 1;
    R(i).at = R(i).st+tt([R(i).xo,R(i).yo],[R(i).xd,R(i).yd]);
    R(i).pass = randi(2,1);
    R(i).ID = old_N + i;
    R(i).status = 0;
    R(i).priority = 0;
end
for i=floor(3*N/4)+1:N
    R(i).xo = 120*rand(1,1);
    R(i).yo = 120*rand(1,1);
    R(i).xd = R(i).xo - 50*rand(1,1);
    R(i).yd = R(i).yo + 50*rand(1,1);
    R(i).st = 1;
    R(i).at = R(i).st+tt([R(i).xo,R(i).yo],[R(i).xd,R(i).yd]);
    R(i).pass = randi(2,1);
    R(i).ID = old_N + i;
    R(i).status = 0;
    R(i).priority = 0;
end

fileID1 = fopen('Requests_Data2.txt','a');
for i=1:N
    %fprintf(fileID1, '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\r\n', [R(i).xo, R(i).yo,...
        %R(i).xd, R(i).yd, R(i).st, R(i).at, R(i).pass, (R(i).ID-1), R(i).status R(i).priority]);
    fprintf(fileID1, '%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\r\n', [R(i).xo, R(i).yo,...
        R(i).xd, R(i).yd, R(i).st, R(i).at, R(i).pass, (R(i).ID-1), R(i).status R(i).priority]);
end
fclose(fileID1);

%load('R_k_k_kkk.mat');
new_R = R; % these are the actual new tasks
new_R_exist = 1; % To say that new tasks have arrived
clear R; % I do this only to not change the name to the code ubove
% initialization of vehicles

% (the vehicles should be generated once, and then just make their status 0 when free)

for i=1:floor(num_v/4);
    V(i).x = 70*rand(1,1);
    V(i).y = 20*rand(1,1);
    V(i).c = max_cap;
    V(i).ID = i; %%%% vehicle ID
    V(i).status = 0;
end
for i=floor(num_v/4)+1:floor(num_v/2)
    V(i).x = 80+20*rand(1,1);
    V(i).y = 80+20*rand(1,1);
    V(i).c = max_cap;
    V(i).ID = i;
    V(i).status = 0;
end
for i=floor(num_v/2)+1:floor(3*num_v/4)
    V(i).x = 20+20*rand(1,1);
    V(i).y = 80+20*rand(1,1);
    V(i).c = max_cap;
    V(i).ID = i;
    V(i).status = 0;
end
for i=floor(3*num_v/4)+1:num_v
    V(i).x = 50+20*rand(1,1);
    V(i).y = 50+100*rand(1,1);
    V(i).c = max_cap;
    V(i).ID = i;
    V(i).status = 0;
end

fileID2 = fopen('Vehicles_Data2.txt','a');
for i=1:num_v
    %fprintf(fileID2, '%6.3f %6.3f %6.3f %6.3f %6.3f\r\n', [V(i).x, V(i).y,...
        %V(i).c, (V(i).ID-1), V(i).status]);
    fprintf(fileID2, '%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\r\n', [V(i).x, V(i).y,...
        V(i).c, (V(i).ID-1), V(i).status]);
end
fclose(fileID2);

%load('V_k_k_kkk.mat');

ALL_V = V;


size1_exist = 0; % This tells to 'data_for_optimal_assignment' whether to add constraints or not
priority_exist = 0;
priority_tasks = []; % This when I do not have tasks which asks for priority. I should make a for cycle
                     % to search for this tasks
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following deals with the introduction of already assigned trips of size 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here I check if I have trips of size 1 that I can again use in the optimization problem at the next call of the algorithm
% [rows_R_ass,clmn_R_ass] = size(size1_assigned); % to the how many already assigned tasks can be reconsidered

% % When already assigned trips of size 1 exist
% if rows_R_ass > 0 % if I have already assigned tasks, I attach them to the new tasks
%     R = [R_ass R]; % new set of tasks, with also the already assigned tasks
%     V = [V_ass New_V]; % new set of vehicles, with also the already assigned vehicles. V here is the set of vehicles that are now free, after the priority assignment phase
%     N = N + length(size1_assigned); % update the value of N since we have requests that comes from previous assignments
% %     num_v = num_v + length(size1_assigned) % this should match the number of AVAILABLE vehicles. But the number of total vehicles is fixed
%     num_v = length([V(:).ID]);
%     size1_exist = 1; % it means that we are taking again in consideratin previously assinged trips
% end

% These should be the R and the V entering the problem
% R = [R_ass R_priority new_R];
% V = [new_V V_ass];

% All the following runs ONLY if I have available vehicles (dim(V)>0) and if there are tasks to be assigned (dim(R)>0)

%%% WHEN NEW TASKS DOES NOT EXIST

% new_R_exist = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tic
run Generate_V_and_R
% toc

tic
run Priority_assignment_priorities
% toc

% tic
run RV_graph_priorities
% toc

% tic
run RTV_graph_priorities
% toc

% tic
run Greedy_assignment_priorities
% toc

% tic
run data_for_optimal_assignment_priorities


% tic
run Rebalancing_priorities
toc

heading = {'Trips cost' 'Vehicle ID' 'Trip ID' 'First task','Second task','Third task','Passengers'};
t = uitable('data',opt_ass_final,'columnname',heading,'columnwidth',{80});
return
run prepare_for_next_call_priorities % ofter toc because this can be done also after having assigned vehicles

L_R_OK % just to have information on the number of assigned requests with respect to the total of requests.
       % Maybe this information is more important in percentage: L_R_OK/N_old*100
       
return % comment this line to obtain plots

% PLOTS

% plot of trips with pick-up (stars) and drop-off (cross) points
figure(1)
plot([R.xo],[R.yo],'*b','linewidth',0.9);
hold on
plot([R.xd],[R.yd],'+r','linewidth',0.9);
plot([R.xo;R.xd],[R.yo;R.yd],'--g');
grid on
legend('pick-up','drop-off');
title('map with tasks and vehicles');
hold on
axis equal

for ii = 1:N
    text(R(ii).xo,R(ii).yo,num2str(ii),'Color','r')
end
hold on
for ii = 1:N
    text(R(ii).xd,R(ii).yd,num2str(ii),'Color','r')
end
xlabel('longitude');
ylabel('latitude');

% plot of vehicles in the map
figure(1)
hold on
plot([V.x],[V.y],'^k','linewidth',0.9)
hold on
for ii = 1:num_v
    text(V(ii).x,V(ii).y,num2str(ii),'Color','k')
end

% RV-graph plot

A = [zeros(N,N) Adj(1:N,N+1:end);
    Adj(N+1:end,1:N),zeros(num_v,num_v)];
Adj_plot = Adj+eye(N+num_v,N+num_v);
xy_nodes = [R(:).xo V(:).x; R(:).yo V(:).y]';
figure(2)
% plot of requests and the links between them
gplot(Adj_plot(1:end-num_v,1:end-num_v),xy_nodes(1:end-num_v,:),'r:*')
hold on
%plot of vehicles and matching with single requests
gplot(A,xy_nodes,'b:o')
legend('pick-ups','vehicles')
title({'RV-graph: shareable requests are connected','feasible vehicles-requests are connected'});
grid on

%%% RTV-graph plot

Adj_plot = Adj+eye(N+v,N+v);

e_r_T_Adj = [zeros(N,N) zeros(N,num_v) RTV_Adj(1:N,N+num_v+1:end);
    zeros(num_v,N+num_v+tot);
    RTV_Adj(N+num_v+1:end,1:N) zeros(tot,num_v+tot)];
e_T_v_Adj = [zeros(N,N+num_v+tot);
    zeros(num_v,N+num_v) RTV_Adj(N+1:N+num_v,N+num_v+1:end);
    zeros(tot,N) RTV_Adj(N+num_v+1:end,N+1:N+num_v) zeros(tot,tot)];

Trips_mat = eye(tot,tot);

xy_nodes = [R(:).xo V(:).x 130*ones(1,tot); R(:).yo V(:).y linspace(0,130,tot)]';

% plot of e(r,T)
figure
gplot(e_r_T_Adj,xy_nodes,'m:o') % matrix to plot only the single requests
hold on
% plot of single requests
gplot(eye(N,N),xy_nodes,'b:*')
% plot of trips
gplot(Trips_mat,xy_nodes(N+num_v+1:end,:),'k:*')
% plot of e(T,v)
gplot(e_T_v_Adj,xy_nodes,'g:^')

legend('pick-ups','requests','trips','vehicles')
title({'RTV-graph: representation of feasible trips'});
grid on