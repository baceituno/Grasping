%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% MIT MCube Lab				   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run cage on a two polygon object
display('Clearing workspace')
clc; clear all; close all;

% reads a planar shape
nv = 4;
p = PlanarShape('PolyS');
samples = 9;
pushers = 4;
angles = linspace(-pi/8,pi/8,samples);

results = {};
k = 1;
for pushers = 4

	% sets-up the optimization program
	planner = MixedIntegerFullCagePlanningProblem(p,pushers,samples,angles);

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

	results{k} = planner;
	k = k + 1;
end

display('results')
r = [];
for i = 1:planner.n_pushers
	planner.vars.p.value(:,i)
	r = [r, planner.vars.p.value(:,i) + planner.vars.p_ref.value(:,(samples+1)/2)];
end

figure(1)
v = [];
for i = 1:length(p.polygons)
	v = p.polygons{i}.v;
	hullIndices = convhull(v(1,:),v(2,:));
	h = patch(v(1,hullIndices),v(2,hullIndices),'b');
	set(h,'facealpha',.1);
	axis square
	hold on
end
r = 0.03;
th = 0:pi/50:2*pi;
for j = 1:planner.n_pushers
	x = planner.vars.p.value(1,j);
	y = planner.vars.p.value(2,j);
	xunit = r * cos(th) + x;
	yunit = r * sin(th) + y;
	h = fill(xunit, yunit,'k');
	set(h,'facealpha',.9);
	axis square
	hold on
end

axis square

% draws the C-space slice
figure(2)
for j = 1:planner.n_pushers
	for i = 1:length(p.polygons)
		v = -p.polygons{i}.v + planner.vars.p.value(:,j);
		h = fill(v(1,:),v(2,:),'y');
		axis square
		set(h,'facealpha',.5);
		hold on
	end
end

axis square