%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% MIT MCube Lab				   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run cage on a two polygon object
display('Clearing workspace')
clc; clear all; close all;

% reads a planar shape
p = PlanarShape('Poly2');
T = 4;
pushers = 4;

% sets-up the optimization program
planner = MixedIntegerCertifiedGraspingProblem(p,pushers,T);

% adds constraints
planner = planner.addNoCollisionConstraint();
planner = planner.addCircleConstraint();
planner = planner.addEnclosingConstraint();
planner = planner.addIntegrationConstraints();
planner = planner.addCspaceContractionConstraints();
planner = planner.addTerminalConstraints();
planner = planner.addYumiKinConstraints();
planner = planner.addObservabilityConstraints();
% adds the cost function
planner = planner.addCostFunction();

% solves the optimization
disp('solving...');
planner = planner.solve();

display('results')
figure(1)
v = [];
for t = 1:planner.nT
	for i = 1:length(p.polygons)
		R = [cos(planner.theta(1,t)),sin(planner.theta(1,t)); -sin(planner.theta(1,t)),cos(planner.theta(1,t))];
		v = R*p.polygons{i}.v;

		hullIndices = convhull(v(1,:),v(2,:));
		h = patch(v(1,hullIndices),v(2,hullIndices),'b');
		set(h,'facealpha',.1);
		axis square
		hold on

		R = [cos(planner.theta(end,t)),sin(planner.theta(end,t)); -sin(planner.theta(end,t)),cos(planner.theta(end,t))];
		v = R*p.polygons{i}.v;
		h = patch(v(1,hullIndices),v(2,hullIndices),'b');
		set(h,'facealpha',.1);
		axis square
		hold on
	end
	r = {};
	for i = 1:planner.n_e
		r{i} = zeros(2,t);
		r{i}(1,:) = planner.vars.p.value(1,i,1:t);
		r{i}(2,:) = planner.vars.p.value(2,i,1:t);
		scatter(r{i}(1,:),r{i}(2,:));
		hold on;
	end
	pause()
end
for i = 1:length(p.polygons)
	v = p.polygons{i}.v;
	hullIndices = convhull(v(1,:),v(2,:));
	h = patch(v(1,hullIndices),v(2,hullIndices),'b');
	set(h,'facealpha',.2);
	axis square
	hold on
end

figure(2)
r = {};
for i = 1:planner.n_e/2
	r{i} = zeros(2,t+1);
	r{i}(1,:) = planner.vars.p.value(1,i,1:t+1);
	r{i}(2,:) = planner.vars.p.value(2,i,1:t+1);
	scatter(r{i}(1,:),r{i}(2,:),'b');
	hold on;
end

for i = planner.n_e/2+1:planner.n_e
	r{i} = zeros(2,t+1);
	r{i}(1,:) = planner.vars.p.value(1,i,1:t+1);
	r{i}(2,:) = planner.vars.p.value(2,i,1:t+1);
	scatter(r{i}(1,:),r{i}(2,:),'r');
	hold on;
end

for i = 1:length(p.polygons)
	v = p.polygons{i}.v;
	hullIndices = convhull(v(1,:),v(2,:));
	h = patch(v(1,hullIndices),v(2,hullIndices),'b');
	set(h,'facealpha',.5);
	axis square
	hold on
end