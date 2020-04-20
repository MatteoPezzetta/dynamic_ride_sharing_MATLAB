# TaskAssignment
Algorithm for optimal task assignment

That algorithm solves the optimal task assignment of a fleet of robots/vehicles for pick-up&delivery tasks.
This first compute the feasible trips from the set of all the reqeusts trying to combine them.
Constraints on a maximum delay tell if tasks can be combined together.

First: an RV-graph (shareability graph) is build;
Second: an RTV-graph is build. It tells which are the feasible trips that will enter the optimization as optimization variables;
Third: a Greedy assignment is computed. It will be used as a starting point for the optimal assignment;
Fourth: the optimal assignment is computed though the formulation of an ILP;
Fifth: a rebalancing phase assign idle vehicle to unassigned requests, when present. 
