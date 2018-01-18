%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function safe_regions = createOctagon()
	% creates a simple octagonagon-shape object
	safe_regions = iris.TerrainRegion.empty();

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [1,0,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [1;0;1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [1,0,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [1;0;-1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [-1,0,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [-1;0;1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [-1,0,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [-1;0;-1], normal);

	A = [1,0,0; -1,0,0]; b = [0.5;-0.5];
	normal = [1,0,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0.5;0;0], normal);

	A = [1,0,0; -1,0,0]; b = [-0.5;0.5];
	normal = [-1,0,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [-0.5;0;0], normal);

	A = [1,0,0; -1,0,0]; b = [1;1];
	normal = [0,0,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;0;0.5], normal);

	A = [1,0,0; -1,0,0]; b = [1;1];
	normal = [0,0,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;0;-0.5], normal);
end