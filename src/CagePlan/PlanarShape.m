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
			% reads the file
			fid = fopen(filename,'r');
			np = fscanf(fid,'%d',1);
			obj.polygons = cell(np);

			iv = 0;
			
			% reads each polygon
			while true
			 	r = [fscanf(fid,'%d',1),fscanf(fid,'%d',1)];
			 	obj.polygons{r(1)} = struct('v', [], 'nv', r(2), 'iv', iv, 'center', []);
			 	iv = iv + r(2);
			 	v = [];
			 	for i = 1:r(2)
			 		verts = [fscanf(fid,'%d',1);fscanf(fid,'%d',1)];
			 		v = [v verts];
			 	end
			 	obj.polygons{r(1)}.v = v;
			 	obj.polygons{r(1)}.center = [mean(v(1,:)),mean(v(2,:))]';
			 	if (r(1) == np); break; end;
			end

			obj.nv = iv;

			% reads each region
			nr = fscanf(fid,'%d',1);
			obj.regions = cell(nr);

			while true
			 	r = [fscanf(fid,'%d',1),fscanf(fid,'%d',1)];
			 	obj.regions{r(1)} = struct('A', [], 'b', []);
			 	A = [];
			 	b = [];
			 	for i = 1:r(2)
			 		Ai = [fscanf(fid,'%d',1),fscanf(fid,'%d',1)];
			 		A = [A;Ai];
			 	end
			 	for i = 1:r(2)
			 		bi = fscanf(fid,'%d',1);
			 		b = [b;bi];
			 	end
			 	obj.regions{r(1)}.A = A;
			 	obj.regions{r(1)}.b = b;
			 	if (r(1) == nr); break; end;
			end

			% reads the G matrix
			n = fscanf(fid,'%d',1);
			m = fscanf(fid,'%d',1);

			obj.G = zeros(n,m);
			for i = 1:n
				for j = 1:m
					obj.G(i,j) = fscanf(fid,'%d',1);
				end
			end
		end
	end
end