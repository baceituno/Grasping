%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run Grasp planner on a box

display('Clearing workspace')
clc; clear all; close all;

display('checking for dependencies')
checkDependency('iris');
checkDependency('lcmgl');
path_handle = addpathTemporary(fileparts(mfilename('fullpath')));

% adds the box polygonal regions
box_size = [0.07;0.3;0.07];
[safe_regions, verts] = createBox(box_size);
% [safe_regions, verts] = createOctahedron(0.04);
% [safe_regions, verts] = createPyramid(0.04);

% use the drake visualizer
use_viz = false;

planner = PlanGraspFromPolygon(safe_regions, 3, struct());

% parses the solution
p = planner.vars.p.value;
f = planner.vars.f_e.value;

G = [];
lambda = [];
for i = 1:planner.n_contacts
	G = [G round([eye(3); crossSkewSymMat(p(:,i))],3)];
	lambda = [lambda; f(:,i)];
end

p
dk = G*lambda
rank(G)

regions = planner.vars.region.value'*[1:length(safe_regions)]';
normals = [];
friction_cones = cell(planner.n_contacts,1);
fc = cell(planner.n_contacts,1);


% computes the cone at flat ground
theta = linspace(0,2*pi,planner.num_edges+1);
theta = theta(1:end-1);
edges_0 = [planner.mu_object*cos(theta);planner.mu_object*sin(theta);ones(1,planner.num_edges)]; 

% computes the normals and the friction cones
for i = 1:planner.n_contacts
	normals = [normals, safe_regions(regions(i)).normal/norm(safe_regions(regions(i)).normal)];
	% computes the cone
	R_fc = rotateVectorToAlign([0;0;1],safe_regions(regions(i)).normal);
	friction_cones{i} = R_fc*edges_0;
	fc{i} = LinearizedFrictionCone(p(:,i),safe_regions(regions(i)).normal,planner.mu_object,R_fc*edges_0);
end

% runs the force adjustment
% optimal = ForceAdjustmentLP(G, normals);
% optimal = optimal.solve();

% % difference between vectors
% cross_diff = mean(mean(cross(optimal.vars.f_e.value,f).^2))

% computes the radius of the maximum epsilon ball in the wrench space
Qw = diag([10;10;10;500;500;500]);
epsilon = computeQ1LinFC(p,[0,0,0]',friction_cones,Qw)

% draws the resulting grasp
if use_viz
	% draws the object
	lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton,'shape');
	lcmgl.glColor3f(0,0,1);
	lcmgl.polyhedron(verts(1,:),verts(2,:),verts(3,:));
	lcmgl.switchBuffers();

	% draws the cones
	for i = 1:planner.n_contacts
		fc{i}.plot(use_viz, 0.07, sprintf('cone_%d',i));
	end
else
	figure(1)
	for i = 1:planner.n_contacts
		fc{i}.plot(use_viz, 0.07, sprintf('cone_%d',i));
		hold on
	end
	shape = alphaShape(verts(1,:)',verts(2,:)',verts(3,:)');
	plot(shape);
end