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

N = 25;
num_v = 10;
max_cap = 4;

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
    R(i).ID = i; %%%% request ID
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
    R(i).ID = i;
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
    R(i).ID = i;
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
    R(i).ID = i;
    R(i).status = 0;
    R(i).priority = 0;
end

% initialization of vehicles

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

% load('R_d2_testing','R')
% load('V_d2_testing','V')

ALL_V = V;

tic

run RV_graph_dynamics1

run RTV_graph_dynamics1

run Greedy_assignment_dynamics1

run data_for_optimal_assignment_dynamics1

run Rebalancing_dynamics1

toc

figure(3)
heading = {'Trips cost' 'Vehicle ID' 'Trip ID' 'First task','Second task'};
t = uitable('data',opt_ass_final,'columnname',heading,'columnwidth',{80});

L_R_OK
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
figure(2)
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
figure(4)
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