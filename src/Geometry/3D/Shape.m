classdef Shape
	properties
		vertices
		faces
		normals
		CoM
	end

	methods
		function obj = Shape(filename)
			% initializes a shape
			file = readObj(filename);

			% assignes the vertices
			obj.vertices = file.v(:,1:3);

			% assign the faces
			obj.faces = file.f.v;

			% assignes the normals
			obj.normals = file.vn;

			% finds the CoM
			obj.CoM = [mean(obj.vertices(:,1));...
					   mean(obj.vertices(:,2));...
					   mean(obj.vertices(:,3))];

			if(size(obj.normals,1) == 0)
				error('model has no normals');
			end
		end

		function safe_regions = getRegions(obj)
			% initializes a cell for all the faces
			nf = size(obj.faces,1);

			% reads all of the faces
			for j = 1:nf
				% adds the vertices
				safe_regions(j).com = obj.CoM;
				safe_regions(j).vertices = [];
				safe_regions(j).normal = zeros(3,1);

				% list the number of vertices
				safe_regions(j).nv = length(obj.faces(j,:));

				% reads the vertices of the faces
				for i = 1:safe_regions(j).nv
					safe_regions(j).vertices = [safe_regions(j).vertices,...
					 							obj.vertices(obj.faces(j,i),:)'];
				end

				v1 = obj.vertices(obj.faces(j,2),:)' - obj.vertices(obj.faces(j,1),:)';
				v2 = obj.vertices(obj.faces(j,3),:)' - obj.vertices(obj.faces(j,2),:)';

				% accumulates for the normal
				safe_regions(j).normal = cross(v1,v2);

				% finds the average normal
				safe_regions(j).normal = safe_regions(j).normal/safe_regions(j).nv;

				% indicates that the region uses V-representation
				safe_regions(j).isV = true;
			end
		end

		function safe_regions = randRegions(obj, N)
			% initializes a cell for all the faces
			nf = size(obj.faces,1);

			% generates random idx
			idx = randi([1,nf],N,1,'int16');

			% reads all of the faces
			for k = 1:N
				% reads the random idx
				j = idx(k);
				safe_regions(k).com = obj.CoM;

				% adds the vertices
				safe_regions(k).vertices = [];
				safe_regions(k).normal = zeros(3,1);

				% list the number of vertices
				safe_regions(k).nv = length(obj.faces(j,:));

				% reads the vertices of the faces
				for i = 1:3
					safe_regions(k).vertices = [safe_regions(k).vertices,...
					 							obj.vertices(obj.faces(j,i),:)'];
				end

				% does the cross product for the normal
				v1 = obj.vertices(obj.faces(j,2),:)' - obj.vertices(obj.faces(j,1),:)';
				v2 = obj.vertices(obj.faces(j,3),:)' - obj.vertices(obj.faces(j,2),:)';

				% accumulates for the normal
				safe_regions(k).normal = cross(v1,v2);

				% finds the average normal
				safe_regions(k).normal = safe_regions(k).normal/safe_regions(k).nv;

				% indicates that the region uses V-representation
				safe_regions(k).isV = true;
			end
		end
	end
end