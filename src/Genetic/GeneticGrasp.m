classdef GeneticGrasp
properties
	safe_regions
	solution
	pos
	filename
end

methods
	function obj = GeneticGrasp(filename)
		% reads object
		obj.filename = filename;
		object = Shape(filename);
		obj.safe_regions = object.getRegions();
	end

	function obj = GeneticGraspRun(filename)
		% reads object
		obj.filename = filename;
		object = Shape(filename);
		obj.safe_regions = object.getRegions();

		% runs the GA
		obj = obj.SolveGrasp();

		% plots the solution
		obj = obj.plotSol();
	end

	function obj = SolveGrasp(obj)
		% computes the number of variables
		NVARS = 3*ceil(log2(length(obj.safe_regions))); 

		% scoring function
		fcn = @(x) obj.CostFunction(x);

		% calls for the genetic algorithns
		lb = zeros(NVARS,1); ub = ones(NVARS,1);
		obj.solution = ga(fcn,NVARS,[],[],[],[],lb,ub,[],[1:NVARS]);
	end

	function obj = plotSol(obj)
		% plots the solution of the grasping problem
		N = ceil(log2(length(obj.safe_regions))); 
		% reads the solution
		idx = [];
		normals = [];
		p = [];

		% read the assigned polygon
		for i = 1:3
			% index values
			i1 = (i-1)*N + 1;
			i2 = i*N;

			% lowest index must be one
			val = bi2de(obj.solution(i1:i2));
			if val == 0; val = 1; end;

			% higgers index must be N
			if val > length(obj.safe_regions); val = length(obj.safe_regions); end;

			idx = [idx val];
		end

		% gets the position of each polygon center
		for i = 1:3
			pos = [mean(obj.safe_regions(idx(i)).vertices(1,:)),... 
				   mean(obj.safe_regions(idx(i)).vertices(2,:)),... 
				   mean(obj.safe_regions(idx(i)).vertices(3,:))]';

			p  = [p pos];
		end

		% computes the cone at flat ground
		theta = linspace(0,2*pi,5);
		theta = theta(1:end-1);
		edges_0 = [cos(theta);sin(theta);ones(1,4)]; 

		% for each finger
		for i = 1:3
			% computes the cone
			R_fc = rotateVectorToAlign([0;0;1],obj.safe_regions(idx(i)).normal);
			fc{i} = LinearizedFrictionCone(p(:,i),obj.safe_regions(idx(i)).normal,1.0,R_fc*edges_0);
		end

		figure(1)
		for i = 1:3
			fc{i}.plot(false, 1, sprintf('cone_%d',i));
			hold on
		end

		b = Shape(obj.filename);
		verts = [];
		if size(verts,1) == 0
			verts = b.vertices';
		end

		shape = alphaShape(verts(1,:)',verts(2,:)',verts(3,:)');
		plot(shape);
	end

	function score = CostFunction(obj, x)
		% reads the index for each polygon
		N = ceil(log2(length(obj.safe_regions))); 
		idx = [];
		p = [];
		normals = [];

		% read the assigned polygon
		for i = 1:3
			% index values
			i1 = (i-1)*N + 1;
			i2 = i*N;

			% lowest index must be one
			val = bi2de(x(i1:i2));
			if val == 0; val = 1; end;

			% higgers index must be N
			if val > length(obj.safe_regions); val = length(obj.safe_regions); end;

			idx = [idx val];
		end

		% gets the position of each polygon center
		for i = 1:3
			pos = [mean(obj.safe_regions(idx(i)).vertices(1,:)),... 
				   mean(obj.safe_regions(idx(i)).vertices(2,:)),... 
				   mean(obj.safe_regions(idx(i)).vertices(3,:))]';

			p  = [p pos];
		end

		% computes the cone at flat ground
		theta = linspace(0,2*pi,5);
		theta = theta(1:end-1);
		edges_0 = [cos(theta);sin(theta);ones(1,4)]; 

		% for each finger
		for i = 1:3
			normals = [normals, obj.safe_regions(idx(i)).normal/norm(obj.safe_regions(idx(i)).normal)];
			% computes the cone
			R_fc = rotateVectorToAlign([0;0;1],obj.safe_regions(idx(i)).normal);
			friction_cones{i} = R_fc*edges_0;
		end

		% computes the score
		Qw = diag([10;10;10;500;500;500]);
		score = computeQ1LinFC(p,[0,0,0]',friction_cones,Qw);

		% inverts in order to fit the argmin formulation
		score = -score;
		obj.pos = p;
	end
end

end