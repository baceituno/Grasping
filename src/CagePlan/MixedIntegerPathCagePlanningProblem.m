classdef MixedIntegerPathCagePlanningProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (MIT MCube Lab)
%           and Hongkai Dai (Toyota Research Institute)
  properties
    n_pushers
    shape
    r = 0.1 % min pusher separation
    nT
    x_G
    dt = 0.1
  end

  methods
    function obj = MixedIntegerPathCagePlanningProblem(shape, n_pushers, nT, xG)
      % Constructs the optimization problem and declares the variables for each contact
      % @param shape: structure with elements:
      %               - shape.polygons : a disjunctive set of convex polygons which union 
      %                                  results in the shape, represented by vertices
      %                                  - shape.polygons{i}.nv = number of vertices of the poligon
      %                                  - shape.polygons{i}.v = a 2xnv array with all the vertices
      %                                  - shape.polygons{i}.iv = sum(shape.polygons(k).nv) from k = 1 to i-1
      %                                  - shape.polygons{i}.center = center of the polygon      
      %               - shape.regions  : a set of convex regions that represent the complement 
      %                                  space of the shape
      %               - shape.nv       : Number of vertices in all the polygons of the shape
      % @param n_pushers: number of fingers of the gripper
      % @param nT: number of time-steps
      % @param xG: goal position of the object

      assert(nargin > 0);
      if nargin < 2
        n_pushers = 3; 
      end
      if nargin < 3
        nT = 10; 
        xG = [1,2]; 
      end

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();
      obj.n_pushers = n_pushers;
      obj.shape = shape;
      obj.x_G = xG;
      obj.nT = nT

      % Pusher locations
      obj = obj.addVariable('p', 'C', [2, obj.n_pushers], -inf, inf);

      obj = obj.addVariable('r_obj', 'C', [2, nT], -inf, inf);
      obj = obj.addVariable('p_hand', 'C', [2, nT], -inf, inf);
      obj = obj.addVariable('dp_hand', 'C', [2, nT], -inf, inf);
      obj = obj.addVariable('ddp_hand', 'C', [2, nT], -inf, inf);
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
        A = obj.shape.regions{r}.A;
        b = obj.shape.regions{r}.b;
        % constrains for each finger
        for j = 1:obj.n_pushers
          Ai = sparse(size(A, 1), obj.nv);
          bi = zeros(size(A, 1), 1);
          Ai(:, obj.vars.p.i(1:2,j)) = A;
          Ai(:, obj.vars.region.i(r,j)) = M;
          bi(:) = b + M;
          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end
      
      % assigns each pusher to one regions
      Aeq = sparse(obj.n_pushers,obj.nv);
      beq = ones(obj.n_pushers,1);
      for j = 1:obj.n_pushers
        Aeq(j, obj.vars.region.i(:,j)) = 1;
      end
      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % ensures that the object does not penetrate the ground
      M = length(obj.shape.polygons);
      for j = 1:obj.nT
        for i = 1:M
          for k = obj.shape.polygons{i}.nv
            Ai = sparse(2,obj.nv);
            bi = obj.shape.polygons{i}.v(:,k);
            Ai(:,obj.vars.r_obj.i(:,j)) = -eye(2)
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
      end      

      % ensures that the pushers keep a minimum separation
      % obj = obj.addVariable('d', 'B', [2, obj.n_pushers], 0, 1);

      % for n = 1:obj.n_pushers
      %   idx_1 = n;
      %   idx_2 = mod(n+1,obj.n_pushers);
      %   if(idx_2 == 0); idx_2 = obj.n_pushers; end;

      %   Ai = sparse(4,obj.nv);
      %   bi = -obj.r + [0,M,M,2*M]';
        
      %   Ai(1,obj.vars.d.i(:,idx_1)) = -M;
      %   Ai(1,obj.vars.p.i(:,idx_1)) = [1,1];
      %   Ai(1,obj.vars.p.i(:,idx_2)) = [-1,-1];

      %   Ai(2,obj.vars.d.i(1,idx_1)) = -M;
      %   Ai(2,obj.vars.d.i(1,idx_1)) = M;
      %   Ai(2,obj.vars.p.i(:,idx_1)) = [-1,1];
      %   Ai(2,obj.vars.p.i(:,idx_2)) = [1,-1];
        
      %   Ai(3,obj.vars.d.i(2,idx_1)) = -M;
      %   Ai(3,obj.vars.d.i(2,idx_1)) = M;
      %   Ai(3,obj.vars.p.i(:,idx_1)) = [1,-1];
      %   Ai(3,obj.vars.p.i(:,idx_2)) = [-1,1];

      %   Ai(4,obj.vars.d.i(:,idx_1)) = M;
      %   Ai(4,obj.vars.p.i(:,idx_1)) = [-1,-1];
      %   Ai(4,obj.vars.p.i(:,idx_2)) = [1,1];
        
      %   obj = obj.addLinearConstraints(Ai, bi, [], []);
      % end
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
      obj = obj.addVariable('weight','C',[obj.n_pushers,2,obj.shape.nv], 0, 1);

      % constrains the values of the G matrix
      for n = 1:obj.n_pushers
        for i = 1:M
          for j = 1:M
            Ai = sparse(1,obj.nv);
            bi = obj.shape.G(i,j);
            Ai(1,obj.vars.G.i(n,i,j)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
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
        if(idx_2 == 0); idx_2 = obj.n_pushers; end;

        % intersection constraint
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

            for k = 1:obj.shape.polygons{i}.nv
              Ai(1:2,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k)) = obj.shape.polygons{i}.v(:,k);
              Ai(3:4,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k)) = -obj.shape.polygons{i}.v(:,k);           
            end

            for l = 1:obj.shape.polygons{j}.nv
              Ai(1:2,obj.vars.weight.i(idx_1,2,obj.shape.polygons{j}.iv+l)) = -obj.shape.polygons{j}.v(:,l);           
              Ai(3:4,obj.vars.weight.i(idx_1,2,obj.shape.polygons{j}.iv+l)) = obj.shape.polygons{j}.v(:,l); 
            end

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end

        % weights must add to 1 if intersecting, 0 otherwise
        for i = 1:M
          Aeq = sparse(2,obj.nv);
          beq = zeros(2,1);

          range_1 = (obj.shape.polygons{i}.iv+1):(obj.shape.polygons{i}.iv + obj.shape.polygons{i}.nv);

          Aeq(1,obj.vars.H.i(idx_1,i,:)) = -1;
          Aeq(1,obj.vars.weight.i(idx_1,1,range_1)) = 1;

          Aeq(2,obj.vars.H.i(idx_1,:,i)) = -1;
          Aeq(2,obj.vars.weight.i(idx_1,2,range_1)) = 1;
          
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end 
      end

      % constrains the existance of a closed circle
      for n = 1:obj.n_pushers
        idx_1 = n;
        idx_2 = mod(n+1,obj.n_pushers);
        if(idx_2 == 0); idx_2 = obj.n_pushers; end;

        for i = 1:M
          for j = 1:M
            Ai = sparse(2,obj.nv);
            bi = [K+1,K-1]';

            % for two pushers the graph must go in and out of the 
            % C-space pushers
            if(obj.n_pushers > 2)
              l_range = 1:M;
            else
              l_range = [1:i-1,i+1:M];
            end

            Ai(1,obj.vars.H.i(idx_2,j,l_range)) = 1;
            Ai(1,obj.vars.G.i(idx_2,j,:)) = 1;

            Ai(2,obj.vars.H.i(idx_2,j,l_range)) = -1;
            Ai(2,obj.vars.G.i(idx_2,j,:)) = -1;
            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end

        for p = 1:M
          for q = 1:M
            Ai = sparse(2,obj.nv);
            bi = [K+1,K-1]';

            r_range = [1:p-1,p+1:M];

            Ai(1,obj.vars.H.i(idx_2,q,:)) = 1;
            Ai(1,obj.vars.G.i(idx_2,q,r_range)) = 1;

            Ai(2,obj.vars.H.i(idx_2,q,:)) = -1;
            Ai(2,obj.vars.G.i(idx_2,q,r_range)) = -1;
            Ai(:,obj.vars.G.i(idx_2,p,q)) = K;

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
      obj = obj.addVariable('F','B',[obj.n_pushers,M,6], 0, 1);
      
      % Defins the XOR slack variable
      obj = obj.addVariable('c','C',[obj.n_pushers*M-1,1], 0, 1);

      % adds MI intersection constraints
      obj = obj.addVariable('beta','C',[obj.n_pushers,M], 0, inf);
      obj = obj.addVariable('lambda','C',[obj.n_pushers,M], 0, 1);
      
      % big K
      K = 10;

      for i = 1:obj.n_pushers
        idx_1 = i;
        idx_2 = mod(i+1,obj.n_pushers);
        if(idx_2 == 0); idx_2 = obj.n_pushers; end
        for j = 1:M
          % requires that each F adds to one
          Aeq = sparse(1,obj.nv);
          beq = 1;
          Aeq(1,obj.vars.F.i(i,j,[1:4,6])) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);

          Ai = sparse(1,obj.nv);
          bi = 0;
          Ai(1,obj.vars.F.i(i,j,1:4)) = 1;
          Ai(1,obj.vars.H.i(i,j,:)) = -1;
          Ai(1,obj.vars.G.i(i,j,:)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Aeq = sparse(1,obj.nv);
          beq = 1;
          Aeq(1,obj.vars.F.i(i,j,6)) = 1;
          Aeq(1,obj.vars.H.i(i,j,:)) = 1;
          Aeq(1,obj.vars.G.i(i,j,:)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end
      end

      % performs XOR sequentially
      n = 1;

      for i = 1:obj.n_pushers
        for j = 1:M
          if j == 1 && i == 1
            obj = obj.addXORConstraint(obj.vars.c.i(1,1),obj.vars.F.i(1,1,1),obj.vars.F.i(1,2,1));
          else 
            if j < M
              obj = obj.addXORConstraint(obj.vars.c.i(n),obj.vars.c.i(n-1,1),obj.vars.F.i(i,j+1,1));
            else
              if i < obj.n_pushers
                obj = obj.addXORConstraint(obj.vars.c.i(n),obj.vars.c.i(n-1,1),obj.vars.F.i(i+1,1,1));
              end
            end
          end
          n = n + 1;
        end
      end

      % requires that F adds to an odd value
      Aeq = sparse(1,obj.nv);
      beq = 1;
      Aeq(1,obj.vars.c.i(end,end)) = 1;
      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % small difference to avoid including edges
      dx = 0.01;

      for n = 1:obj.n_pushers
        idx_1 = n;
        idx_2 = mod(n+1,obj.n_pushers);
        if(idx_2 == 0); idx_2 = obj.n_pushers; end
        for i = 1:M
          for j = 1:M
            % adds the constraint when the 
            % polygon is not intersecting
            Ai = sparse(4,obj.nv);
            bi = 2*K*ones(4,1);

            Ai(1:2,obj.vars.p.i(2,idx_1)) = -1;
            Ai(3:4,obj.vars.p.i(1,idx_1)) = [1;-1];
            Ai(:,obj.vars.G.i(idx_1,i,j)) = K;
            Ai(:,obj.vars.F.i(idx_1,i,1)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(2);
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(2);
            if obj.shape.polygons{i}.center(1) > obj.shape.polygons{j}.center(1)
              bi(3,1) = bi(3,1) - obj.shape.polygons{j}.center(1) - dx;
              bi(4,1) = bi(4,1) + obj.shape.polygons{i}.center(1) - dx;
            else 
              bi(3,1) = bi(3,1) - obj.shape.polygons{i}.center(1) - dx;
              bi(4,1) = bi(4,1) + obj.shape.polygons{j}.center(1) - dx;
            end

            obj = obj.addLinearConstraints(Ai, bi, [], []);
            
            % if the ray passes above the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1:2,obj.vars.p.i(2,idx_1)) = 1;
            Ai(1:2,obj.vars.G.i(idx_1,i,j)) = K;
            Ai(1:2,obj.vars.F.i(idx_1,i,2)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(2) - dx;
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(2) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the right of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1)) = -1;
            Ai(2,obj.vars.p.i(1,idx_1)) = -1;
            Ai(1:2,obj.vars.G.i(idx_1,i,j)) = K;
            Ai(1:2,obj.vars.F.i(idx_1,i,3)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the left of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1)) = 1;
            Ai(2,obj.vars.p.i(1,idx_1)) = 1;
            Ai(1:2,obj.vars.G.i(idx_1,i,j)) = K;
            Ai(1:2,obj.vars.F.i(idx_1,i,4)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % adds the constraint constrains when the 
            % polygon is intersecting with other pusher
            
            % the ray intersects the segment
            Ai = sparse(4,obj.nv);
            bi = 3*K*ones(4,1);

            Ai(1,obj.vars.p.i(2,idx_1)) = -1;
            Ai(2,obj.vars.p.i(2,idx_2)) = -1;
            Ai(3,obj.vars.p.i(1,idx_1)) = 1;
            Ai(4,obj.vars.p.i(1,idx_2)) = -1;
            Ai(:,obj.vars.F.i(idx_1,i,1)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;
            Ai(:,obj.vars.F.i(idx_1,i,5)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(2);
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(2);
            bi(3,1) = bi(3,1) - obj.shape.polygons{i}.center(1) - dx;
            bi(4,1) = bi(4,1) + obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            Ai = sparse(4,obj.nv);
            bi = 2*K*ones(4,1);

            Ai(1,obj.vars.p.i(2,idx_1)) = -1;
            Ai(2,obj.vars.p.i(2,idx_2)) = -1;
            Ai(3,obj.vars.p.i(1,idx_1)) = -1;
            Ai(4,obj.vars.p.i(1,idx_2)) = 1;
            Ai(:,obj.vars.F.i(idx_1,i,1)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;
            Ai(:,obj.vars.F.i(idx_1,i,5)) = -K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(2);
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(2);
            bi(3,1) = bi(3,1) + obj.shape.polygons{i}.center(1) - dx;
            bi(4,1) = bi(4,1) - obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
 
            % the ray passes above the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(2,idx_1)) = 1;
            Ai(2,obj.vars.p.i(2,idx_2)) = 1;
            Ai(:,obj.vars.F.i(idx_1,i,2)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(2) - dx;
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(2) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the right of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1)) = -1;
            Ai(2,obj.vars.p.i(1,idx_2)) = -1;
            Ai(:,obj.vars.F.i(idx_1,i,3)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the left of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1)) = 1;
            Ai(2,obj.vars.p.i(1,idx_2)) = 1;
            Ai(:,obj.vars.F.i(idx_1,i,4)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end
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

    function obj = addIntegrationConstraints(obj)
      % defines the initial conditions as zero
      Aeq = sparse(6, obj.nv);
      beq = zeros(6, 1);

      Aeq(1:2, obj.vars.p_hand.i(:,1)) = eye(2);
      Aeq(3:4, obj.vars.dp_hand.i(:,1)) = eye(2);
      Aeq(5:6, obj.vars.ddp_hand.i(:,1)) = eye(2);

      obj = obj.addLinearConstraints([],[],Aeq,beq);

      % direct transcription of the hand motion
      for j = 2:obj.nT
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, obj.vars.p_hand.i(:,j)) = eye(2);
        Aeq(:, obj.vars.p_hand.i(:,j-1)) = -eye(2);
        Aeq(:, obj.vars.dp_hand.i(:,j)) = -obj.dt*eye(2);

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, obj.vars.dp_hand.i(:,j)) = eye(2);
        Aeq(:, obj.vars.dp_hand.i(:,j-1)) = -eye(2);
        Aeq(:, obj.vars.ddp_hand.i(:,j)) = -obj.dt*eye(2);

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end
    end

    function obj = addCostFunction(obj)
      % Minimizes the separation between pushers
      for j = 1:obj.n_pushers-1
        Qi = sparse(obj.nv,obj.nv);
        Qi(obj.vars.p.i(:,j),obj.vars.p.i(:,j)) = eye(2);
        Qi(obj.vars.p.i(1,j),obj.vars.p.i(1,j+1)) = -2;
        Qi(obj.vars.p.i(2,j),obj.vars.p.i(2,j+1)) = -2;
        Qi(obj.vars.p.i(:,j+1),obj.vars.p.i(:,j+1)) = eye(2);
        obj = obj.addCost(Qi,[],[]);
      end

      % minimizes acceleration
      for j = 1:obj.nT
        Qi = sparse(obj.nv,obj.nv);
        Qi(obj.vars.ddp_hand.i(:,j),obj.vars.ddp_hand.i(:,j)) = eye(2);
        obj = obj.addCost(Qi,[],[]);
      end

      % minimizes distance to goal
      Qi = sparse([],[],[],obj.nv,obj.nv,2);
      ci = sparse(obj.nv, 1);
      Qi(obj.vars.p_hand.i(:,obj.nT), obj.vars.p_hand.i(:,obj.nT)) = 1;
      ci(obj.vars.p_hand.i(:,obj.nT)) = -2*x_G;
      objcon_i = x_G'*x_G;
      obj = obj.addCost(Qi, ci, objcon_i);
    end
    % end of methods
  end
end