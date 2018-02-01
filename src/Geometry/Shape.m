classdef Shape
	properties
		vertices
		faces
		normals
	end

	methods
		function obj = Shape(filename)
			% initializes a shape
			file = readObj(filename);

			% assignes the vertices
			obj.vertices = file.v;

			% assign the faces
			obj.faces = file.f;

			% assignes the normals
			obj.normals = file.vn;
		end

		function safe_regions = regions(obj)
			% initializes a cell for all the faces
			nf = size(obj.faces,1);
			safe_regions = cell(1,nf);

			% reads all of the faces
			for j = 1:nf
				% adds the vertices
				safe_regions{j}.vertices = [];
				safe_regions{j}.normal = zeros(3,1);

				% list the number of vertices
				safe_regions{j}.nv = length(obj.faces(j,:));

				% reads the vertices of the faces
				for i = 1:safe_regions{j}.nv
					safe_regions{j}.vertices = [safe_regions{j}.vertices,...
					 							obj.vertices(obj.faces(j,i),:)'];
					 % accumulates for the normal
					safe_regions{j}.normal += obj.vertices(obj.faces(j,i),:)';
				end

				% finds the average normal
				safe_regions{j}.normal /= safe_regions{j}.nv;

				% indicates that the region uses V-representation
				safe_regions{j}.isV = true;
			end
		end

		function safe_regions = rand_regions(obj, N)
			% initializes a cell for all the faces
			nf = size(obj.faces,1);
			safe_regions = cell(1,N);

			% generates random idx
			idx = randi([1,nf],N,1,'int16');

			% reads all of the faces
			for k = 1:N
				% reads the random idx
				j = idx(k);

				% adds the vertices
				safe_regions{j}.vertices = [];
				safe_regions{j}.normal = zeros(3,1);

				% list the number of vertices
				safe_regions{j}.nv = length(obj.faces(j,:));

				% reads the vertices of the faces
				for i = 1:safe_regions{j}.nv
					safe_regions{j}.vertices = [safe_regions{j}.vertices,...
					 							obj.vertices(obj.faces(j,i),:)'];
					 % accumulates for the normal
					safe_regions{j}.normal += obj.vertices(obj.faces(j,i),:)';
				end

				% finds the average normal
				safe_regions{j}.normal /= safe_regions{j}.nv;

				% indicates that the region uses V-representation
				safe_regions{j}.isV = true;
			end
		end
	end
end