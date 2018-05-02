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

% sets-up the optimization for 2 fingers
planner = MixedIntegerCagePlanningProblem(p,2);

% adds the convex region constraints
planner = planner.addNoCollisionConstraint();

% adds the closed graph constraints
planner = planner.addCircleConstraint();

% adds the enclosing constraints
planner = planner.addEnclosingConstraint();

% solves the optimization
planner = planner.solve();

display('results')
planner.vars.p.value