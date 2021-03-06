%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% MIT MCube Lab				   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run Grasp planner on a box

display('Clearing workspace')
clc; clear all; close all;

display('checking for dependencies')
checkDependency('iris');
checkDependency('lcmgl');
path_handle = addpathTemporary(fileparts(mfilename('fullpath')));

differ = [];
improv = [];
for i = 1:20
	% adds the box polygonal regions
	safe_regions = createBox();
	planner = PlanGraspFromPolygon(safe_regions, 3, struct('lin_sides',20,'quad_approx',false));

	% parses the solution
	p = planner.vars.p.value;
	f = planner.vars.f_e.value;

	G = [];
	lambda = []
	for i = 1:planner.n_contacts
		G = [G round([eye(3); crossSkewSymMat(p(:,i))],3)];
		lambda = [lambda; f(:,i)];
	end

	p
	dk = G*lambda
	G
	rank(G)

	regions = planner.vars.region.value'*[1:length(safe_regions)]';
	normals = [];
	for i = 1:planner.n_contacts
		normals = [normals, safe_regions(regions(i)).normal/norm(safe_regions(regions(i)).normal)];
	end

	% runs the force adjustment
	optimal = ForceAdjustmentLP(G, normals);
	optimal = optimal.solve();

	% improvement percentage
	improvement = (optimal.vars.epsilon.value - planner.vars.epsilon.value)
	diff = mean(mean((optimal.vars.f_e.value - f).^2))
	cross_diff = mean(mean(cross(optimal.vars.f_e.value,f).^2))

	improv = [improv improvement];
	differ = [differ cross_diff];

	% adds the box polygonal regions
	safe_regions = createPyramid();
	planner = PlanGraspFromPolygon(safe_regions, 3, struct('lin_sides',20,'quad_approx',false));

	% parses the solution
	p = planner.vars.p.value;
	f = planner.vars.f_e.value;

	G = [];
	lambda = []
	for i = 1:planner.n_contacts
		G = [G round([eye(3); crossSkewSymMat(p(:,i))],3)];
		lambda = [lambda; f(:,i)];
	end

	p
	dk = G*lambda
	G
	rank(G)

	regions = planner.vars.region.value'*[1:length(safe_regions)]';
	normals = [];
	for i = 1:planner.n_contacts
		normals = [normals, safe_regions(regions(i)).normal/norm(safe_regions(regions(i)).normal)];
	end

	% runs the force adjustment
	optimal = ForceAdjustmentLP(G, normals);
	optimal = optimal.solve();

	% improvement percentage
	improvement = (optimal.vars.epsilon.value - planner.vars.epsilon.value)
	diff = mean(mean((optimal.vars.f_e.value - f).^2))
	cross_diff = mean(mean(cross(optimal.vars.f_e.value,f).^2))

	improv = [improv improvement];
	differ = [differ cross_diff];

	% adds the box polygonal regions
	safe_regions = creatOctahedron();
	planner = PlanGraspFromPolygon(safe_regions, 3, struct('lin_sides',20,'quad_approx',false));

	% parses the solution
	p = planner.vars.p.value;
	f = planner.vars.f_e.value;

	G = [];
	lambda = []
	for i = 1:planner.n_contacts
		G = [G round([eye(3); crossSkewSymMat(p(:,i))],3)];
		lambda = [lambda; f(:,i)];
	end

	p
	dk = G*lambda
	G
	rank(G)

	regions = planner.vars.region.value'*[1:length(safe_regions)]';
	normals = [];
	for i = 1:planner.n_contacts
		normals = [normals, safe_regions(regions(i)).normal/norm(safe_regions(regions(i)).normal)];
	end

	% runs the force adjustment
	optimal = ForceAdjustmentLP(G, normals);
	optimal = optimal.solve();

	% improvement percentage
	improvement = (optimal.vars.epsilon.value - planner.vars.epsilon.value)
	diff = mean(mean((optimal.vars.f_e.value - f).^2))
	cross_diff = mean(mean(cross(optimal.vars.f_e.value,f).^2))

	improv = [improv improvement];
	differ = [differ cross_diff];

	% adds the box polygonal regions
	safe_regions = createBall();
	planner = PlanGraspFromPolygon(safe_regions, 3, struct('lin_sides',20,'quad_approx',false));

	% parses the solution
	p = planner.vars.p.value;
	f = planner.vars.f_e.value;

	G = [];
	lambda = []
	for i = 1:planner.n_contacts
		G = [G round([eye(3); crossSkewSymMat(p(:,i))],3)];
		lambda = [lambda; f(:,i)];
	end

	p
	dk = G*lambda
	G
	rank(G)

	regions = planner.vars.region.value'*[1:length(safe_regions)]';
	normals = [];
	for i = 1:planner.n_contacts
		normals = [normals, safe_regions(regions(i)).normal/norm(safe_regions(regions(i)).normal)];
	end

	% runs the force adjustment
	optimal = ForceAdjustmentLP(G, normals);
	optimal = optimal.solve();

	% improvement percentage
	improvement = (optimal.vars.epsilon.value - planner.vars.epsilon.value)
	diff = mean(mean((optimal.vars.f_e.value - f).^2))
	cross_diff = mean(mean(cross(optimal.vars.f_e.value,f).^2))

	improv = [improv improvement];
	differ = [differ cross_diff];
end