run run_optimization_problem_dynamics1
run prepare_for_next_call_priorities

% update ALL_V

new_R_exist = 0;

% free some vehicles here
% Need to account for the case in which there is no availabel vehicle: -->
% --> so when priority_ass does not exist.

run Generate_V_and_R
run Priority_assignment_priorities
run RV_graph_priorities
run RTV_graph_priorities
run Greedy_assignment_priorities
run data_for_optimal_assignment_priorities

run Rebalancing_priorities

% I should enforce the assignment of trips of size 2 when I have tasks with priority.

% I can for example decrease the cost of assign trips of size 2 by assigning low cost
% to the trips of size 2 that contain tasks with priorities