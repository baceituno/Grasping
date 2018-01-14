%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function safe_regions = createBox()
	% creates a simple box
	safe_regions = iris.TerrainRegion.empty();

	A = [1,0,0; -1,0,0]; b = [1;-1];
	normal = [1.5,0,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [1;0;0], normal);

	A = [1,0,0; -1,0,0]; b = [-1;1];
	normal = [-1,0,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [-1;0;0], normal);

	A = [0,1,0; 0,-1,0]; b = [1;-1];
	normal = [0,1,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;1;0], normal);

	A = [0,1,0; 0,-1,0]; b = [-1;1];
	normal = [0,-1,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;-1;0], normal);

	A = [0,0,1; 0,0,-1]; b = [1;-1];
	normal = [0,0,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;0;1], normal);

	A = [0,0,1; 0,0,-1]; b = [-1;1];
	normal = [0,0,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;0;-1], normal);
end