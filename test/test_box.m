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
box_size = [0.07;0.07;0.3];
[safe_regions, verts] = createBox(box_size);
% [safe_regions, verts] = createBall(0.04);

planner = PlanGraspFromPolygon(safe_regions, 3, struct('lin_sides',4,'quad_approx',false,... 
													   'palm_pos', [-0.03,0,0]'));

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
G
rank(G)

regions = planner.vars.region.value'*[1:length(safe_regions)]';
normals = [];
friction_cones = cell(planner.n_contacts,1);

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
end

% runs the force adjustment
optimal = ForceAdjustmentLP(G, normals);
optimal = optimal.solve();

% improvement percentage
improvement = (optimal.vars.epsilon.value - planner.vars.epsilon.value)
cross_diff = mean(mean(cross(optimal.vars.f_e.value,f).^2))

% computes the radius of the maximum epsilon ball in the wrench space
Qw = diag([10;10;10;500;500;500]);
epsilon = computeQ1LinFC(p,[0,0,0]',friction_cones,Qw)