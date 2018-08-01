%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run cage on a two polygon object
display('Clearing workspace')
clc; clear all; close all;

% reads the 2 triange shape
p = PlanarShape('Poly2');
samples = 9;

% sets-up the optimization for 2 fingers
planner = MixedIntegerFullCagePlanningProblem(p,3,samples);

% adds the slice constraints
planner = planner.addSliceConstraints();

% adds the limit orientation constraints
planner = planner.addLimitOrientationConstraints();

% adds the convex region constraints
planner = planner.addNoCollisionConstraint();

% adds the closed graph constraints
planner = planner.addCircleConstraint();

% adds the enclosing constraints
planner = planner.addEnclosingConstraint();

% adds the continuous boundary constriants
planner = planner.addContinuousBoundaryVariationConstraints();

% adds the cost function
% planner = planner.addCostFunction();

% solves the optimization
disp('solving...');
planner = planner.solve();

display('results')
for i = 1:planner.n_pushers
	planner.vars.p.value(:,i) + planner.vars.p_ref.value(:,(samples+1)/2)
end