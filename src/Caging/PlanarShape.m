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
			np = fscanf(fid,'%f',1);
			obj.polygons = cell(np);

			iv = 0;
			
			% reads each polygon
			while true
			 	r = [fscanf(fid,'%f',1),fscanf(fid,'%f',1)];
			 	obj.polygons{r(1)} = struct('v', [], 'nv', r(2), 'iv', iv, 'center', []);
			 	iv = iv + r(2);
			 	v = [];
			 	for i = 1:r(2)
			 		verts = [fscanf(fid,'%f',1);fscanf(fid,'%f',1)];
			 		v = [v verts];
			 	end
			 	obj.polygons{r(1)}.v = v;
			 	obj.polygons{r(1)}.center = [mean(v(1,:)),mean(v(2,:))]';
			 	if (r(1) == np); break; end;
			end

			obj.nv = iv;

			% reads each region
			nr = fscanf(fid,'%f',1);
			obj.regions = cell(nr);

			while true
			 	r = [fscanf(fid,'%f',1),fscanf(fid,'%f',1)];
			 	obj.regions{r(1)} = struct('A', [], 'b', []);
			 	A = [];
			 	b = [];
			 	for i = 1:r(2)
			 		Ai = [fscanf(fid,'%f',1),fscanf(fid,'%f',1)];
			 		A = [A;Ai];
			 	end
			 	for i = 1:r(2)
			 		bi = fscanf(fid,'%f',1);
			 		b = [b;bi];
			 	end
			 	obj.regions{r(1)}.A = A;
			 	obj.regions{r(1)}.b = b;
			 	if (r(1) == nr); break; end;
			end

			% reads the G matrix
			n = fscanf(fid,'%f',1);
			m = fscanf(fid,'%f',1);

			obj.G = zeros(n,m);
			for i = 1:n
				for j = 1:m
					obj.G(i,j) = fscanf(fid,'%f',1);
				end
			end

			% % reads the line segments and determines the complementarity
			nl = fscanf(fid,'%f',1);
			obj.lines = cell(nl);

			for i = 1:nl
				obj.lines{i} = struct('v1', [], 'v2', [], 'opp', [], 'non_cod', [], 'angle', [], 'isCV', false);
				obj.lines{i}.v1 = [fscanf(fid,'%f',1);fscanf(fid,'%f',1)];
				obj.lines{i}.v2 = [fscanf(fid,'%f',1);fscanf(fid,'%f',1)];

				% computes the angle of the segment, given that the vertexes are counted counterclockwise
				dy = abs(obj.lines{i}.v2(2)-obj.lines{i}.v1(2));
				dx = abs(obj.lines{i}.v2(1)-obj.lines{i}.v1(1));

				if obj.lines{i}.v2(2) > obj.lines{i}.v1(2)
					if obj.lines{i}.v2(1) > obj.lines{i}.v1(1)
						obj.lines{i}.angle = -atan(dy/dx);
					elseif obj.lines{i}.v2(1) == obj.lines{i}.v1(1)
						obj.lines{i}.angle = 0;
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
			for i = 1:nl
				% facet normal
				if obj.lines{i}.isCV
					if i == 1
						deg = obj.lines{i+1}.angle + obj.lines{end}.angle;
						normal = [cos(deg/2);sin(deg/2)];
						normal = normal/norm(normal);
					elseif i == nl
						deg = obj.lines{i}.angle + obj.lines{1}.angle;
						normal = [cos(deg/2);sin(deg/2)];
						normal = normal/norm(normal);						
					else
						deg = obj.lines{i}.angle + obj.lines{i+1}.angle;
						normal = [cos(deg/2);sin(deg/2)];
						normal = normal/norm(normal);
					end
				else
					normal = [cos(obj.lines{i}.angle);sin(obj.lines{i}.angle)];
					normal = normal/norm(normal);
				end 

				obj.lines{i}.normal = normal;

				% opposition of the faces
				opposite = [];

				if obj.lines{i}.isCV
					if i == 1
						opposite = [i+2:nl-1];
					elseif i == nl
						opposite = [2:i-2];
					else
						opposite = [1:i-2,i+2:nl];
					end
				else
					for j = 1:nl
						if obj.lines{j}.isCV && j ~= i
							opposite = [opposite, j];
						end
						if abs(obj.lines{i}.angle - obj.lines{j}.angle) == pi
							opposite = [opposite, j];
						end
					end
				end

				obj.lines{i}.opp = opposite;

				% parallelity of the faces
				paral = [];

				if obj.lines{i}.isCV
					if i == 1
						kdx1 = nl;
						kdx2 = i+1;
					elseif i == nl
						kdx1 = i-i;
						kdx2 = 1;
					else
						kdx1 = i-1;
						kdx2 = i+1;
					end
					
					for j = 1:nl
						if obj.lines{j}.isCV
							idx1 = j - 1;
							idx2 = j + 1;

							if idx1 < 1; idx1 = nl; end;
							if idx2 > nl; idx2 = 1; end;

							dif_angs1 = abs(obj.lines{kdx1}.angle - obj.lines{idx1}.angle);
							dif_angs2 = abs(obj.lines{kdx1}.angle - obj.lines{idx2}.angle);

							dif_angs3 = abs(obj.lines{kdx2}.angle - obj.lines{idx1}.angle);
							dif_angs4 = abs(obj.lines{kdx2}.angle - obj.lines{idx2}.angle);

							if (dif_angs1 == 0) || (dif_angs1 == pi) || (dif_angs2 == 0) || (dif_angs2 == pi) || (dif_angs3 == 0) || (dif_angs3 == pi) || (dif_angs4 == 0) || (dif_angs4 == pi)
								paral = [paral, j];
							end
						else
							if kdx1 < 1; kdx1 = 1; end;
							if kdx2 > nl; kdx2 = nl; end;
							dif_angs1 = abs(obj.lines{kdx1}.angle - obj.lines{j}.angle);
							dif_angs2 = abs(obj.lines{kdx2}.angle - obj.lines{j}.angle);
							if (dif_angs1 == 0) || (dif_angs1 == pi) || (dif_angs2 == 0) || (dif_angs2 == pi)
								paral = [paral, j];
							end
						end
					end

				else
					for j = 1:nl
						if obj.lines{j}.isCV
							idx1 = j - 1;
							idx2 = j + 1;

							if idx1 < 1
								idx1 = nl;
							end

							if idx2 > nl
								idx2 = 1;
							end

							dif_angs1 = abs(obj.lines{i}.angle - obj.lines{idx1}.angle);
							dif_angs2 = abs(obj.lines{i}.angle - obj.lines{idx2}.angle);

							if (dif_angs1 == 0) || (dif_angs1 == pi) || (dif_angs2 == 0) || (dif_angs2 == pi)
								paral = [paral, j];
							end
						else
							dif_angs = abs(obj.lines{i}.angle - obj.lines{j}.angle);
							if (dif_angs == 0) || (dif_angs == pi)
								paral = [paral, j];
							end
						end
					end
				end

				obj.lines{i}.parallel_fac = paral;

				% co-directionaly of the faces
				noncod = [];

				if obj.lines{i}.isCV
					if i == 1
						noncod = [i+2:nl-1];
					elseif i == nl
						noncod = [2:i-2];
					else
						noncod = [1:i-2,i+2:nl];
					end
				else
					for j = 1:nl
						if obj.lines{j}.isCV && j ~= i
							if j == 1
								if i ~= nl && i ~= 2
									noncod = [noncod, j];
								end
							elseif j == nl
								if i ~= nl-1 && i ~= 1
									noncod = [noncod, j];
								end
							else
								if i ~= j-1 && i ~= j+1
									noncod = [noncod, j];
								end
							end
						end
						if ~obj.lines{j}.isCV
							dif_angs = abs(obj.lines{i}.angle - obj.lines{j}.angle);
							if (dif_angs > 0) && (dif_angs ~= pi);
								noncod = [noncod, j];
							end
						end
					end
				end

				obj.lines{i}.non_cod = noncod;

				% co-equality of the faces
				noneq = [];

				if obj.lines{i}.isCV
					if i == 1
						noneq = [i+2:nl-1];
					elseif i == nl
						noneq = [2:i-2];
					else
						noneq = [1:i-2,i+2:nl];
					end
				else
					for j = 1:nl
						if obj.lines{j}.isCV && j ~= i
							if j == 1
								if i ~= nl && i ~= 2
									noneq = [noneq, j];
								end
							elseif j == nl
								if i ~= nl-1 && i ~= 1
									noneq = [noneq, j];
								end
							else
								if i ~= j-1 && i ~= j+1
									noneq = [noneq, j];
								end
							end
						end
						if ~obj.lines{j}.isCV
							dif_angs = abs(obj.lines{i}.angle - obj.lines{j}.angle);
							if (dif_angs > 0);
								noneq = [noneq, j];
							end
						end
					end
				end

				obj.lines{i}.non_eq = noneq;
			end
		end
	end
end