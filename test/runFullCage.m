%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run cage on a two polygon object
display('Clearing workspace')
clc; clear all; close all;

% reads a planar shape
p = PlanarShape('Poly1');
samples = 9;
pushers = 5;

% sets-up the optimization program
planner = MixedIntegerFullCagePlanningProblem(p,pushers,samples);

% adds the slice constraints
planner = planner.addSliceConstraints();

% adds the limit orientation constraints
planner = planner.addLimitOrientationConstraints();

% adds the convex region constraints
planner = planner.addNoCollisionConstraint();

% adds the closed graph constraints
planner = planner.addLoopConstraint();

% adds the enclosing constraints
planner = planner.addEnclosingConstraint();

% adds the continuous boundary constraints
planner = planner.addContinuousBoundaryVariationConstraints();

% adds the cost function
% planner = planner.addCostFunction();

% solves the optimization
disp('solving...');
planner = planner.solve();

display('results')
r = [];
for i = 1:planner.n_pushers
	planner.vars.p.value(:,i)
	r = [r, planner.vars.p.value(:,i) + planner.vars.p_ref.value(:,(samples+1)/2)];
end