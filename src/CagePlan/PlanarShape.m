classdef PlanarShape
	properties
		polygons
		regions
		nv = 0
		G
	end

	methods
		function obj = PlanarShape(filename)
			% initializes a shape
			M = length(obj.polygons);
			for i = 1:M
				obj.nv += obj.polygons(i).iv;
			end
		end
	end
end