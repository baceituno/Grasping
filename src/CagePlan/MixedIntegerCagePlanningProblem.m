classdef MixedIntegerCagePlanningProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (MIT MCube Lab)
%           and Hongkai Dai (Toyota Research Institute)
  properties
    n_pushers
    shape
  end

  methods
    function obj = MixedIntegerCagePlanningProblem(shape, n_pushers)
      % Constructs the optimization problem and declares the variables for each contact
      % @param shape: structure with elements:
      %               - shape.polygons : a disjunctive set of convex polygons which union 
      %                                  results in the shape, represented by vertices
      %               - shape.regions  : a set of convex regions that represent the complement 
      %                                  space of the shape
      %               - shape.nv       : Number of vertices in all the polygons of the shape
      % @param n_pushers: number of fingers of the gripper

      assert(nargin > 0);
      if nargin < 2; n_pushers = 3; end

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();
      obj.n_pushers = n_pushers;
      obj.shape = shape;

      % Pusher locations
      obj = obj.addVariable('p', 'C', [2, obj.n_pushers], -inf, inf);
    end

    function obj = addNoCollisionConstraint(obj)
      % constrains the pushers to lie within the complement of the object, 
      % represented through a set of convex polygonal regions using 
      % a logarithmic number of integer variables, as decribed in
      % Modeling Disjunctive Constraints with a Logarithmic Number of Binary
      % Variables and Constraints by J. Vielma and G. Nemhauser, 2011.

      % defines the region assignment variable
      nr = length(obj.shape.regions);
      obj = obj.addVariable('region', 'C', [nr, obj.n_pushers], 0, 1);

      % defines the log2 binary variables
      ny = ceil(log2(nr));
      obj = obj.addVariable('y', 'B', [ny, obj.n_pushers], 0, 1);

      % adds the sos1 constraint on lambda
      for i = 1:obj.n_pushers
        Aeq = sparse(1, obj.nv);
        beq = 1;
        Aeq(1, obj.vars.region.i(:,i)) = 1;
        obj = obj.addLinearConstraints([], [], Aeq, beq);        
      end

      % defines the gray coding on base 2
      codes = grayCodes(2,ny);
      codes = codes(1:nr,:);

      % adds the constraints on coding
      for i = 1:obj.n_pushers
        % for each digit
        for j = 1:ny
          Ai = sparse(2, obj.nv);
          bi = [0;1];
          % for each lambda
          for k = 1:nr
            if codes(k,j) == 1
              Ai(1, obj.vars.region.i(k,i)) = 1;
            elseif codes(k,j) == 0
              Ai(2, obj.vars.region.i(k,i)) = 1;
            else
              error('codes need to be 0 or 1');
            end  
          end
          Ai(1, obj.vars.y.i(j,i)) = -1;
          Ai(2, obj.vars.y.i(j,i)) = 1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end

      % big-M
      M = 10;

      % regions assignment constraint
      for r = 1:nr
        A = obj.shape.regions(r).A;
        b = obj.shape.regions(r).b;
        s = size(A, 1);
        % constrains for each finger
        for j = 1:obj.n_pushers
          Ai = sparse(size(A, 1), obj.nv);
          bi = zeros(size(A, 1), 1);

          Ai(:, obj.vars.p.i(:,j)) = A;
          Ai(:, obj.vars.region.i(r,j)) = M;
          bi(:) = b + M;
          
          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end
      
      % assigns each pusher to one regions
      Aeq = sparse(obj.n_pushers, obj.nv);
      beq = zeros(obj.n_pushers, 1);
      for j = 1:obj.n_pushers
        Aeq(j, obj.vars.region.i(:,j)) = 1;
      end
      obj = obj.addLinearConstraints([], [], Aeq, beq);
    end

    function obj = addCircleConstraint(obj)
      % Add mixed-integer constraints that require that 
      % the graph formed by the intersection of the 
      % C-space pushers forms a cyclic graph.

      % Defines the polygon interserction matrices
      M = length(obj.shape.polygons);
      obj = obj.addVariable('H','B',[obj.n_pushers,M,M], 0, 1);
      obj = obj.addVariable('G','B',[obj.n_pushers,M,M], 0, 1);

      % Defines the vertex multipliers
      % TODO: check if the size is correct
      obj = obj.addVariable('weight','C',[obj.n_pushers,obj.shape.nv], 0, inf);

      % Gets the G matrix
      G = object.shape.G;
      assert(G == G');

      % constrains the values of the G matrix
      for i = 1:M
        for j = 1:M
          Ai = sparse(1,obj.nv);
          bi = G(i,j);
          Ai(1,obj.vars.G.i(i,j)) = 1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end

      % big-K value
      K = 10;

      % requires that each pair of pushers intersect at least once
      for n = 1:obj.n_pushers
        Aeq = sparse(1,obj.nv);
        beq = 1;
        Aeq(1,obj.vars.H.i(n,:,:)) = 1;
        obj = obj.addLinearConstraints([], [], Aeq, beq);
      end

      % constrains the intersection between pushers
      for n = 1:obj.n_pushers
        idx_1 = n;
        idx_2 = mod(n+1,obj.n_pushers);
        for i = 1:M
          for j = 1:M
            % intersection of the pushers
            Ai = sparse(4,obj.nv);
            bi = K*ones(4,1);

            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;

            Ai(1:2,obj.vars.p.i(:,idx_1)) = eye(2);
            Ai(1:2,obj.vars.p.i(:,idx_2)) = -eye(2);

            Ai(3:4,obj.vars.p.i(:,idx_1)) = -eye(2);
            Ai(3:4,obj.vars.p.i(:,idx_2)) = eye(2);

            for k = 1:obj.shape.polygons(i).nv
              Ai(1:2,obj.vars.weight.i(n,obj.shape.polygons(i).iv+k)) = obj.shape.polygons(i).v(:,k);
              Ai(3:4,obj.vars.weight.i(n,obj.shape.polygons(i).iv+k)) = -obj.shape.polygons(i).v(:,k);           
            end

            for l = 1:obj.shape.polygons(j).nv
              Ai(1:2,obj.vars.weight.i(n,obj.shape.polygons(j).iv+l)) = -obj.shape.polygons(j).v(:,l);           
              Ai(3:4,obj.vars.weight.i(n,obj.shape.polygons(j).iv+l)) = obj.shape.polygons(j).v(:,l); 
            end

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % sum of the weights
            Aeq = sparse(2,obj.nv);
            beq = zeros(2,1);

            Aeq(:,obj.vars.H.i(n,i,j)) = -1;

            range_1 = obj.shape.polygons(i).iv:(obj.shape.polygons(i).iv + obj.shape.polygons(i).nv);
            Aeq(1,obj.vars.weight.i(n,range_1)) = 1;

            range_2 = obj.shape.polygons(j).iv:(obj.shape.polygons(j).iv + obj.shape.polygons(j).nv);
            Aeq(2,obj.vars.weight.i(n,range_2)) = 1;

            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end
        end 
      end

      % constrains the existance of a closed circle
      for n = 1:obj.n_pushers
        idx_1 = n;
        idx_2 = mod(n+1,obj.n_pushers);
        for i = 1:M
          for j = 1:M
            Ai = sparse(2,obj.nv);
            bi = [0;1];

            % for two pushers the graph must go in and out of the 
            % C-space pushers
            if(obj.n_pushers > 2)
              l_range = 1:M;
            else
              l_range = [1:i-1,i+1:M];
            end

            Ai(1,obj.vars.H.i(idx_2,j,l_range)) = -1;
            Ai(1,obj.vars.G.i(idx_2,:)) = -1;
            Ai(1,obj.vars.H.i(idx_1,i,j)) = 1;

            Ai(2,obj.vars.H.i(idx_2,j,l_range)) = 1;
            Ai(2,obj.vars.G.i(idx_2,:)) = 1;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end

        for p = 1:M
          for q = 1:M
            Ai = sparse(2,obj.nv);
            bi = [1;0];

            r_range = [1:p-1, p+1:M];

            Ai(1,obj.vars.H.i(idx_2,q,:)) = 1;
            Ai(1,obj.vars.G.i(idx_2,q,r_range)) = 1;

            Ai(2,obj.vars.H.i(idx_2,q,:)) = -1;
            Ai(2,obj.vars.G.i(idx_2,q,r_range)) = -1;
            Ai(2,obj.vars.G.i(idx_2,p,q)) = 1;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end
    end

    function obj = addEnclosingConstraint(obj)
      % constrains the origin to lie inside the circle formed by the 
      % pushers in the C-space of the object by using

      % Defines the polygon interserction matrices
      M = length(obj.shape.polygons);
      obj = obj.addVariable('F','B',[obj.n_pushers,M], 0, 1);
      
      % Defins teh XOR slack variable
      obj = obj.addVariable('c','C',[obj.n_pushers,M-1], 0, 1);
      
      % performs the XOR sequentially
      obj = obj.addXORConstraint(obj.vars.c.i(1,1),obj.vars.F.i(1,1),obj.vars.F.i(1,2));
      k = 1;
      l = 1;

      for i = 1:obj.n_pushers
        for j = 1:M-1
          obj = obj.addXORConstraint(obj.vars.c.i(i,j),obj.vars.c.i(k,l),obj.vars.F.i(i,j+1));
          k = i;
          l = j;
        end
      end

      % requires that F adds to an odd value
      Aeq = sparse(1,obj.nv);
      beq = 1;
      Aeq(1,obj.vars.c.i(end,end)) = 1;
      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % TODO: add MI intersection constraints
      b = [0;1];
      obj = obj.addVariable('beta','C',[obj.n_pushers,M], 0, inf);
      obj = obj.addVariable('lambda','C',[obj.n_pushers,M], 0, 1);
    end

    function obj = addXORConstraint(obj,c,b1,b2)
      % constrains that c = b1 XOR b2 as:
      % c <= b1 + b2
      % c >= b1 - b2
      % c >= b2 - b1
      % c <= 2 - b2 - b1
      % with c in [0,1]

      Ai = sparse(4,obj.nv);
      bi = [zeros(1,3), 2]';

      Ai(1,b1) = -1;
      Ai(1,b2) = -1;
      Ai(1,c) = 1;

      Ai(2,b1) = -1;
      Ai(2,b2) = 1;
      Ai(2,c) = -1;

      Ai(3,b1) = 1;
      Ai(3,b2) = -1;
      Ai(3,c) = -1;

      Ai(4,b1) = 1;
      Ai(4,b2) = 1;
      Ai(4,c) = 1;

      obj = obj.addLinearConstraints(Ai, bi, [], []);      
    end

    function obj = addCostFunction(obj)
      % Minimizes a cost function defined by the user
      % TODO:
    end
    % end of methods
  end
end