%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [safe_regions, verts] = createOctahedron(radius)
	% creates a simple diamond-shape object
	safe_regions = iris.TerrainRegion.empty();

	verts = 3*radius*[1,0,0;0,1,0;-1,0,0;0,0,1;0,-1,0;0,0,-1]';

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [1,1,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[1;1;1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [1,1,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[1;1;-1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [1,-1,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[1;-1;1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [1,-1,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[1;-1;-1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [-1,1,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[-1;1;1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [-1,1,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[-1;1;-1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [-1,-1,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[-1;-1;1], normal);

	A = [0,0,1; 0,0,-1]; b = [0;0];
	normal = [-1,-1,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], radius*[-1;-1;-1], normal);
end