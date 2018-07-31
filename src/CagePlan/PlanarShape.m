classdef PlanarShape
	properties
		polygons
		regions
		lines
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

			% % reads the line segments and determines the complementarity
			nl = fscanf(fid,'%d',1);
			obj.lines = cell(nl);

			for i = 1:nl
				obj.lines{i} = struct('v1', [], 'v2', [], 'opp', [], 'non_cod', [], 'angle', [], 'isCV', false);
				obj.lines{i}.v1 = [fscanf(fid,'%d',1);fscanf(fid,'%d',1)];
				obj.lines{i}.v2 = [fscanf(fid,'%d',1);fscanf(fid,'%d',1)];

				% computes the angle of the segment, given that the vertexes are counted counterclockwise
				dy = abs(obj.lines{i}.v2(2)-obj.lines{i}.v1(2));
				dx = abs(obj.lines{i}.v2(1)-obj.lines{i}.v1(1));

				if obj.lines{i}.v2(2) > obj.lines{i}.v1(2)
					if obj.lines{i}.v2(1) > obj.lines{i}.v1(1)
						obj.lines{i}.angle = -atan(dy/dx);
					else
						obj.lines{i}.angle = atan(dy/dx);
					end
				elseif obj.lines{i}.v2(2) == obj.lines{i}.v1(2)
					if obj.lines{i}.v2(1) > obj.lines{i}.v1(1)
						obj.lines{i}.angle = -pi/2;
					else
						obj.lines{i}.angle = pi/2;
					end
				else
					if obj.lines{i}.v2(1) > obj.lines{i}.v1(1)
						obj.lines{i}.angle = atan(dy/dx)-pi;
					else
						obj.lines{i}.angle = atan(dy/dx)+pi/2;
					end
				end

				% checks if the segment is a concave vertex
				if dx == 0 & dy == 0
					obj.lines{i}.isCV = true;
					obj.lines{i}.angle = 0;
				end
			end

			% determines the propierties of the faces
			% for i = 1:nl
			% 	% opposition of the faces
			% 	opposite = [];

			% 	if obj.lines{j}.isCV
			% 		opposite = 1:nl;
			% 	else
			% 		for j = 1:nl
			% 			if obj.lines{j}.isCV && j ~= i
			% 				opposite = [opposite, j];
			% 			end
			% 			if abs(obj.lines{i}.angle - obj.lines{j}.angle) == pi
			% 				opposite = [opposite, j];
			% 			end
			% 		end
			% 	end

			% 	obj.lines{i}.opp = opposite;

			% 	% co-directionaly of the faces
			% 	noncod = [];

			% 	if obj.lines{j}.isCV
			% 		noncod = 1:nl;
			% 	else
			% 		for j = 1:nl
			% 			if obj.lines{j}.isCV && j ~= i
			% 				noncod = [noncod, j];
			% 			end
			% 			dif_angs = abs(obj.lines{i}.angle - obj.lines{j}.angle);
			% 			if dif_angs > 0;
			% 				noncod = [opposite, j];
			% 			end
			% 		end
			% 	end

			% 	obj.lines{i}.non_cod = non_cod;
			% end
		end
	end
end