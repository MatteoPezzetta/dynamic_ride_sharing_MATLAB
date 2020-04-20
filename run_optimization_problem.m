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

tic

run RV_graph_2

run RTV_graph_2

run Greedy_assignment_2

run data_for_optimal_assignment

run Rebalancing

toc

heading = {'Trips cost' 'Vehicle ID' 'Trip ID' 'First task','Second task'};
t = uitable('data',opt_ass_final,'columnname',heading,'columnwidth',{80});

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