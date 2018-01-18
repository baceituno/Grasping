%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [safe_regions, verts] = createBox(box_size)
	% creates a simple box
	safe_regions = iris.TerrainRegion.empty();
	verts = repmat(box_size/2,1,8).*[1 1 1 1 -1 -1 -1 -1;1 1 -1 -1 1 1 -1 -1;1 -1 1 -1 1 -1 1 -1];

	box_size = box_size/2;
	A = [1,0,0; -1,0,0]; b = [box_size(1);-box_size(1)];
	normal = [1,0,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [box_size(1);0;0], normal);

	A = [1,0,0; -1,0,0]; b = [-box_size(1);box_size(1)];
	normal = [-1,0,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [-box_size(1);0;0], normal);

	A = [0,1,0; 0,-1,0]; b = [box_size(2);-box_size(2)];
	normal = [0,1,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;box_size(2);0], normal);

	A = [0,1,0; 0,-1,0]; b = [-box_size(2);box_size(2)];
	normal = [0,-1,0]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;-box_size(2);0], normal);

	A = [1,0,0; -1,0,0]; b = [box_size(1);box_size(1)];
	normal = [0,0,1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;0;box_size(3)], normal);

	A = [1,0,0; -1,0,0]; b = [box_size(1);box_size(1)];
	normal = [0,0,-1]';
	safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], [0;0;-box_size(3)], normal);
end