classdef MixedIntegerFullCagePlanningProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (MIT MCube Lab)
%           and Hongkai Dai (Toyota Research Institute)
  properties
    n_pushers
    shape
    r = 0.1 % min pusher separation
    n_samples
  end

  methods
    function obj = MixedIntegerFullCagePlanningProblem(shape, n_pushers, n_samples)
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
      %               - shape.lines    : a set of line segments covering the shape
      %               - shape.nv       : Number of vertices in all the polygons of the shape
      % @param n_pushers: number of fingers of the gripper

      assert(nargin > 0);
      assert(mod(n_samples,2) ~= 0);
      if nargin < 2; n_pushers = 3; end

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();
      obj.n_pushers = n_pushers;
      obj.shape = shape;
      obj.n_samples = n_samples;

      % Pusher formation
      obj = obj.addVariable('p', 'C', [2, obj.n_pushers], -inf, inf);

      % Hand location
      obj = obj.addVariable('p_ref', 'C', [2, obj.n_samples], -inf, inf);

      % Critical Slice
      obj = obj.addVariable('Th', 'B', [1, obj.n_samples], 0, 1);
    end

    function obj = addSliceConstraints(obj)
      % constrains that circle enclosing constraints are not added
      % after the object reaches a limit orientation (in both orientations)

      % for each slice with negative orientation
      for i = 1:(obj.n_samples-1)/2
        Ai = sparse(1,obj.nv);
        bi = zeros(1,1);

        % Th_i < Th_i+1
        Ai(1,obj.vars.Th.i(1,i)) = 1;
        Ai(1,obj.vars.Th.i(1,i+1)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end

      % for each slice with positive orientation
      for i = (obj.n_samples+1)/2:obj.n_samples
        Ai = sparse(1,obj.nv);
        bi = zeros(1,1);

        % Th_i-1 < Th_i
        Ai(1,obj.vars.Th.i(1,i-1)) = -1;
        Ai(1,obj.vars.Th.i(1,i)) = 1;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end
    end

    function obj = addLimitOrientationConstraints(obj)
      % Constrains the conditions required for the object to reach a limit orientation,
      % depending on the number of pushers. A limit orientation is that for which the
      % object is fully impobilized by the pushers, with a zero-area free-space.


      % defines the line assignment variable
      nl = length(obj.shape.lines);
      obj = obj.addVariable('line', 'C', [nl, obj.n_pushers, obj.n_samples], 0, 1);
      obj = obj.addVariable('lambda_l', 'C', [nl, obj.n_pushers, obj.n_samples], 0, 1);

      % defines the log2 binary variables
      ny = ceil(log2(nl));
      obj = obj.addVariable('y2', 'B', [ny, obj.n_pushers, obj.n_samples], 0, 1);

      % defines the gray coding on base 2
      codes = grayCodes(2,ny);
      codes = codes(1:nl,:);

      % adds the constraints on coding
      for w = 1:obj.n_samples 
        for i = 1:obj.n_pushers
          % for each digit
          for j = 1:ny
            Ai = sparse(2, obj.nv);
            bi = [0;1];
            % for each lambda
            for k = 1:nl
              if codes(k,j) == 1
                Ai(1, obj.vars.line.i(k,i,w)) = 1;
              elseif codes(k,j) == 0
                Ai(2, obj.vars.line.i(k,i,w)) = 1;
              else
                error('codes need to be 0 or 1');
              end  
            end
            Ai(1, obj.vars.y2.i(j,i,w)) = -1;
            Ai(2, obj.vars.y2.i(j,i,w)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end

      % big K 
      K = 100;

      % line assignment
      for s = 1:obj.n_samples
        for i = 1:obj.n_pushers
          % each pusher can be at most in one line
          Ai = sparse(1, obj.nv);
          bi = ones(1, 1);
          Ai(1,obj.vars.line.i(:,i,s)) = 1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          for l = 1:nl
            Ai = sparse(4, obj.nv);
            bi = K*ones(4, 1);

            Ai(1:2,obj.vars.p.(:,i)) = eye(2);
            Ai(1:2,obj.vars.p_ref.(:,s)) = eye(2);

            Ai(1:2,obj.vars.lambda_l.i(l,i,s)) = -obj.shape.lines{l}.v1;
            Ai(1:2,obj.vars.lambda_l.i(l,i,s)) = obj.shape.lines{l}.v2;
            bi(1:2) = bi(1:2) + obj.shape.lines{l}.v2;

            Ai(3:4,obj.vars.p.(:,i)) = -eye(2);
            Ai(3:4,obj.vars.p_ref.(:,s)) = -eye(2);

            Ai(3:4,obj.vars.lambda_l.i(l,i,s)) = obj.shape.lines{l}.v1;
            Ai(3:4,obj.vars.lambda_l.i(l,i,s)) = -obj.shape.lines{l}.v2;
            bi(3:4) = bi(3:4) - obj.shape.lines{l}.v2; 
            
            Ai(:,obj.vars.line.i(l,i,s)) = K;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end

      % depending on the number of fingers determines what is a critical orientation
      if n_pushers == 2
        for s = 1:obj.n_samples
          % both fingers must be in contact to have a critical orientation
          Ai = sparse(2, obj.nv);
          bi = zeros(2, 1);

          Ai(:,obj.vars.Th.(1,s)) = 1;
          if s < (obj.n_samples + 1)/2
            Ai(:,obj.vars.Th.(1,s+1)) = -1;
          else 
            Ai(:,obj.vars.Th.(1,s-1)) = -1;
          end
          Ai(1,obj.vars.line.i(:,1,s)) = -1;
          Ai(2,obj.vars.line.i(:,2,s)) = -1;

          obj = obj.addLinearConstraints(Ai, bi, [], []);   

          % only contact with opposite faces counts
          for i = 1:obj.n_pushers
            idx_1 = i;
            idx_2 = mod(i+1,obj.n_pushers);
            if(idx_2 == 0); idx_2 = obj.n_pushers; end;
            for l = 1:nl
              Ai = sparse(1, obj.nv);
              bi = zeros(1, 1);

              Ai(1,obj.vars.line.i(l,idx_1,s)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.opp,idx_2,s)) = -1;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end 
        end
      else
        % there are three types of critical orientations
        obj = obj.addVariable('crit_type', 'B', [2, obj.n_samples], 0, 1);

        for s = 1:obj.n_samples
          Ai = sparse(1, obj.nv);
          bi = zeros(1,1);          

          Ai(1,obj.vars.Th.i(1,s)) = 1;
          Ai(1,obj.vars.crit_type.i(:,s)) = -1;

          if s < (obj.n_samples + 1)/2
            Ai(1,obj.vars.Th.(1,s+1)) = -1;
          else 
            Ai(1,obj.vars.Th.(1,s-1)) = -1;
          end

          obj = obj.addLinearConstraints(Ai, bi, [], []);

          % For 2 pushers in contact:
          % only contact with opposite faces counts
          for i = 1:obj.n_pushers
            idx_1 = i;
            idx_comp = [1:idx_1-1,idx+1:obj.n_pushers];

            for l = 1:nl
              Ai = sparse(1, obj.nv);
              bi = K;

              Ai(1,obj.vars.line.i(l,idx_1,s)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.opp,idx_comp,s)) = -1;

              Ai(1,obj.vars.crit_type.i(1,s)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end 

          % For 3 pushers or more in contact:
          % only contact with non-codirectional faces counts
          for i = 1:obj.n_pushers
            for l = 1:nl
              Ai = sparse(1, obj.nv);
              bi = K*ones(1,1);

              Ai(1,obj.vars.line.i(l,i,s)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.non_cod,:,s)) = -1/2;

              Ai(1,obj.vars.crit_type.i(2,s)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
      end
    end

    function rotation = Rotate(obj,sample)
      % determines the rotation matrix for the sample
      ang = -pi + 2*pi*(sample-1)/obj.n_samples;
      Rotate = [cos(ang), -sin(ang); sin(ang), cos(ang)];
    end

    function obj = addNoCollisionConstraint(obj)
      % constrains the pushers to lie within the complement of the object, 
      % represented through a set of convex polygonal regions using 
      % a logarithmic number of integer variables, as decribed in
      % Modeling Disjunctive Constraints with a Logarithmic Number of Binary
      % Variables and Constraints by J. Vielma and G. Nemhauser, 2011.

      % defines the region assignment variable
      nr = length(obj.shape.regions);
      obj = obj.addVariable('region', 'C', [nr, obj.n_pushers, obj.n_samples], 0, 1);

      % defines the log2 binary variables
      ny = ceil(log2(nr));
      obj = obj.addVariable('y', 'B', [ny, obj.n_pushers, obj.n_samples], 0, 1);

      % defines the gray coding on base 2
      codes = grayCodes(2,ny);
      codes = codes(1:nr,:);

      % adds the constraints on coding
      for w = 1:obj.n_samples 
        for i = 1:obj.n_pushers
          % for each digit
          for j = 1:ny
            Ai = sparse(2, obj.nv);
            bi = [0;1];
            % for each lambda
            for k = 1:nr
              if codes(k,j) == 1
                Ai(1, obj.vars.region.i(k,i,w)) = 1;
              elseif codes(k,j) == 0
                Ai(2, obj.vars.region.i(k,i,w)) = 1;
              else
                error('codes need to be 0 or 1');
              end  
            end
            Ai(1, obj.vars.y.i(j,i,w)) = -1;
            Ai(2, obj.vars.y.i(j,i,w)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end

      % big-M
      M = 10;

      % regions assignment constraint
      for i = 1:obj.n_samples
        for r = 1:nr
          A = obj.shape.regions{r}.A;
          b = obj.shape.regions{r}.b;
          % constrains for each finger
          for j = 1:obj.n_pushers
            Ai = sparse(size(A, 1), obj.nv);
            bi = zeros(size(A, 1), 1);
            Ai(:, obj.vars.p.i(1:2,j)) = A;
            Ai(:, obj.vars.p_ref.i(1:2,i)) = A;
            Ai(:, obj.vars.region.i(r,j,i)) = M;
            bi(:) = b + M;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end

      % assigns each pusher to one region when Th = 1 (sos1 constraint)
      for i = 1:obj.n_samples
        Ai = sparse(2*obj.n_pushers,obj.nv);
        bi = ones(2*obj.n_pushers,1);
        for j = 1:obj.n_pushers
          Ai(2*j-1, obj.vars.region.i(:,j,i)) = 1;
          Ai(2*j-1, obj.vars.Th.i(1,i)) = -M;
          bi(2*j-1,1) = 1;

          Ai(2*j, obj.vars.region.i(:,j,i)) = -1;
          Ai(2*j, obj.vars.Th.i(1,i)) = -M;
          bi(2*j,1) = -1;
        end
        obj = obj.addLinearConstraints([], [], Aeq, beq);
      end
    end

    function obj = addCircleConstraint(obj)
      % Add mixed-integer constraints that require that 
      % the graph formed by the intersection of the 
      % C-space pushers forms a cyclic graph.

      % Defines the polygon interserction matrices
      M = length(obj.shape.polygons);
      obj = obj.addVariable('H','B',[obj.n_pushers,M,M,obj.n_samples], 0, 1);
      obj = obj.addVariable('G','B',[obj.n_pushers,M,M,obj.n_samples], 0, 1);

      % Defines the vertex multipliers
      obj = obj.addVariable('weight','C',[obj.n_pushers,2,obj.shape.nv,obj.n_samples], 0, 1);

      % constrains the values of the G matrix
      for s = 1:obj.n_samples
        for n = 1:obj.n_pushers
          for i = 1:M
            for j = 1:M
              Ai = sparse(1,obj.nv);
              bi = obj.shape.G(i,j);
              Ai(1,obj.vars.G.i(n,i,j,s)) = 1;
              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
      end

      % big-K value
      K = 10;

      % requires that each pair of pushers intersect at least once
      for s = 1:obj.n_samples
        for n = 1:obj.n_pushers
          Ai = sparse(2,obj.nv);
          bi = [1;-1];
          Aeq(1,obj.vars.H.i(n,:,:,s)) = 1;
          Aeq(-1,obj.vars.H.i(n,:,:,s)) = -1;
          Aeq(:,obj.vars.Th.i(1,s)) = -K;
          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end

      % constrains the intersection between pushers
      for s = 1:obj.n_samples
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

              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;

              Ai(1:2,obj.vars.p.i(:,idx_1)) = eye(2);
              Ai(1:2,obj.vars.p.i(:,idx_2)) = -eye(2);

              Ai(3:4,obj.vars.p.i(:,idx_1)) = -eye(2);
              Ai(3:4,obj.vars.p.i(:,idx_2)) = eye(2);

              for k = 1:obj.shape.polygons{i}.nv
                Ai(1:2,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,s)) = obj.Rotate(s)*obj.shape.polygons{i}.v(:,k);
                Ai(3:4,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,s)) = -obj.Rotate(s)*obj.shape.polygons{i}.v(:,k);           
              end

              for l = 1:obj.shape.polygons{j}.nv
                Ai(1:2,obj.vars.weight.i(idx_1,2,obj.shape.polygons{j}.iv+l,s)) = -obj.Rotate(s)*obj.shape.polygons{j}.v(:,l);           
                Ai(3:4,obj.vars.weight.i(idx_1,2,obj.shape.polygons{j}.iv+l,s)) = obj.Rotate(s)*obj.shape.polygons{j}.v(:,l); 
              end

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end

          % weights must add to 1 if intersecting, 0 otherwise
          for i = 1:M
            Aeq = sparse(2,obj.nv);
            beq = zeros(2,1);

            range_1 = (obj.shape.polygons{i}.iv+1):(obj.shape.polygons{i}.iv + obj.shape.polygons{i}.nv);

            Aeq(1,obj.vars.H.i(idx_1,i,:,s)) = -1;
            Aeq(1,obj.vars.weight.i(idx_1,1,range_1,s)) = 1;

            Aeq(2,obj.vars.H.i(idx_1,:,i,s)) = -1;
            Aeq(2,obj.vars.weight.i(idx_1,2,range_1,s)) = 1;
            
            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end 
        end
      end

      % constrains the existance of a closed circle
      for s = 1:obj.n_samples
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

              Ai(1,obj.vars.H.i(idx_2,j,l_range,s)) = 1;
              Ai(1,obj.vars.G.i(idx_2,j,:,s)) = 1;

              Ai(2,obj.vars.H.i(idx_2,j,l_range,s)) = -1;
              Ai(2,obj.vars.G.i(idx_2,j,:,s)) = -1;
              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end

          for p = 1:M
            for q = 1:M
              Ai = sparse(2,obj.nv);
              bi = [K+1,K-1]';

              r_range = [1:p-1,p+1:M];

              Ai(1,obj.vars.H.i(idx_2,q,:,s)) = 1;
              Ai(1,obj.vars.G.i(idx_2,q,r_range,s)) = 1;

              Ai(2,obj.vars.H.i(idx_2,q,:,s)) = -1;
              Ai(2,obj.vars.G.i(idx_2,q,r_range,s)) = -1;
              Ai(:,obj.vars.G.i(idx_2,p,q,s)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
      end 
    end

    function obj = addContinuousBoundaryVariationConstraints(obj)
      % constrains that the geometric boundary of the circle must change continuously
      % this is done by requiring that the intersection polygong between C-space pushers
      % is not broken between adjacents samples

      obj = obj.addVariable('weight2','C',[obj.n_pushers,obj.shape.nv,obj.n_samples], 0, 1);
      obj = obj.addVariable('weight3','C',[obj.n_pushers,obj.shape.nv,obj.n_samples], 0, 1);

      % number of polygons
      M = length(obj.shape.polygons);

      % big-K number
      K = 100;

      for s = 1:obj.n_samples-1
        for n = 1:obj.n_pushers
          idx_1 = n;
          idx_2 = mod(n+1,obj.n_pushers);
          if(idx_2 == 0); idx_2 = obj.n_pushers; end

          % weight must add to 1
          for i = 1:M
            Aeq = sparse(2,obj.nv);
            beq = zeros(2,1);

            range_1 = (obj.shape.polygons{i}.iv+1):(obj.shape.polygons{i}.iv + obj.shape.polygons{i}.nv);

            Aeq(1,obj.vars.H.i(idx_1,i,:,s)) = -1;
            Aeq(1,obj.vars.weight2.i(idx_1,range_1,s)) = 1;
            
            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end

          % intersection must not break when rotating from one slice to the next
          for i = 1:M
            for j = 1:M
              Ai = sparse(4,obj.nv);
              bi = K*ones(4,1);

              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;

              % requires that the intersection point is in the polygon during adjacent samples
              Ai(1:2,obj.vars.p.i(:,idx_1)) = eye(2);
              Ai(1:2,obj.vars.p.i(:,idx_2)) = -eye(2);

              Ai(3:4,obj.vars.p.i(:,idx_1)) = -eye(2);
              Ai(3:4,obj.vars.p.i(:,idx_2)) = eye(2);

              for k = 1:obj.shape.polygons{i}.nv
                Ai(1:2,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,s)) = obj.Rotate(s+1)*obj.shape.polygons{i}.v(:,k);
                Ai(3:4,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,s)) = -obj.Rotate(s+1)*obj.shape.polygons{i}.v(:,k);           
              end

              for l = 1:obj.shape.polygons{j}.nv
                Ai(1:2,obj.vars.weight2.i(idx_1,obj.shape.polygons{j}.iv+l,s)) = -obj.Rotate(s+1)*obj.shape.polygons{j}.v(:,l);           
                Ai(3:4,obj.vars.weight2.i(idx_1,obj.shape.polygons{j}.iv+l,s)) = obj.Rotate(s+1)*obj.shape.polygons{j}.v(:,l); 
              end

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              Ai = sparse(4,obj.nv);
              bi = K*ones(4,1);

              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;

              % requires that the arc formed by the intersection point is in the polygon during adjacent samples
              Ai(1:2,obj.vars.p.i(:,idx_1)) = eye(2);
              Ai(1:2,obj.vars.p.i(:,idx_2)) = -eye(2);

              Ai(3:4,obj.vars.p.i(:,idx_1)) = -eye(2);
              Ai(3:4,obj.vars.p.i(:,idx_2)) = eye(2);

              for k = 1:obj.shape.polygons{i}.nv
                Ai(1:2,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,s)) = (eye(2) + PlanarRotMat(-pi/2)*tan(delta/2))*obj.Rotate(s+1)*obj.shape.polygons{i}.v(:,k);
                Ai(3:4,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,s)) = -(eye(2) + PlanarRotMat(-pi/2)*tan(delta/2))*obj.Rotate(s+1)*obj.shape.polygons{i}.v(:,k);           
              end

              for l = 1:obj.shape.polygons{j}.nv
                Ai(1:2,obj.vars.weight3.i(idx_1,obj.shape.polygons{j}.iv+l,s)) = -obj.Rotate(s+1)*obj.shape.polygons{j}.v(:,l);           
                Ai(3:4,obj.vars.weight3.i(idx_1,obj.shape.polygons{j}.iv+l,s)) = obj.Rotate(s+1)*obj.shape.polygons{j}.v(:,l); 
              end

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
      end

    end
    

    function obj = addEnclosingConstraint(obj)
      % constrains the origin to lie inside the circle formed by the 
      % pushers in the C-space of the object by using

      % Defines the polygon interserction matrices
      M = length(obj.shape.polygons);
      obj = obj.addVariable('F','B',[obj.n_pushers,M,6,obj.n_samples], 0, 1);
      
      % Defins the XOR slack variable
      obj = obj.addVariable('c','C',[obj.n_pushers*M-1,1,obj.n_samples], 0, 1);

      % adds MI intersection constraints
      obj = obj.addVariable('beta','C',[obj.n_pushers,M,obj.n_samples], 0, inf);
      obj = obj.addVariable('lambda','C',[obj.n_pushers,M,obj.n_samples], 0, 1);
      
      % big K
      K = 10;

      for s = 1:obj.n_samples
        for i = 1:obj.n_pushers
          idx_1 = i;
          idx_2 = mod(i+1,obj.n_pushers);
          if(idx_2 == 0); idx_2 = obj.n_pushers; end
          for j = 1:M
            % requires that each F adds to one
            Aeq = sparse(1,obj.nv);
            beq = 1;
            Aeq(1,obj.vars.F.i(i,j,[1:4,6],s)) = 1;
            obj = obj.addLinearConstraints([], [], Aeq, beq);

            Ai = sparse(1,obj.nv);
            bi = 0;
            Ai(1,obj.vars.F.i(i,j,1:4,s)) = 1;
            Ai(1,obj.vars.H.i(i,j,:,s)) = -1;
            Ai(1,obj.vars.G.i(i,j,:,s)) = -1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            Aeq = sparse(1,obj.nv);
            beq = 1;
            Aeq(1,obj.vars.F.i(i,j,6,s)) = 1;
            Aeq(1,obj.vars.H.i(i,j,:,s)) = 1;
            Aeq(1,obj.vars.G.i(i,j,:,s)) = 1;
            obj = obj.addLinearConstraints([], [], Aeq, beq);
          end
        end
      end

      % performs XOR sequentially
      for s = 1:obj.n_samples
        n = 1;
        for i = 1:obj.n_pushers
          for j = 1:M
            if j == 1 && i == 1
              obj = obj.addXORConstraint(obj.vars.c.i(1,1,s),obj.vars.F.i(1,1,1,s),obj.vars.F.i(1,2,1,s));
            else 
              if j < M
                obj = obj.addXORConstraint(obj.vars.c.i(n,1,s),obj.vars.c.i(n-1,1,s),obj.vars.F.i(i,j+1,1,s));
              else
                if i < obj.n_pushers
                  obj = obj.addXORConstraint(obj.vars.c.i(n,1,s),obj.vars.c.i(n-1,1,s),obj.vars.F.i(i+1,1,1,s));
                end
              end
            end
            n = n + 1;
          end
        end
      end

      % requires that F adds to an odd value when Th = 1
      for s = 1:obj,n_samples
        Ai = sparse(2,obj.nv);
        bi = [1;-1];
        Ai(1,obj.vars.c.i(end,end,s)) = 1;
        Ai(2,obj.vars.c.i(end,end,s)) = -1;
        Ai(:,obj.vars.Th.i(s,1)) = -K;
        obj = obj.addLinearConstraints([], [], Aeq, beq);
      end

      % small difference to avoid including edges
      dx = 0.01;

      for s = 1:obj.n_samples
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
              Ai(1:2,obj.vars.p_ref.i(2,s)) = -1;
              Ai(3:4,obj.vars.p.i(1,idx_1)) = [1;-1];
              Ai(3:4,obj.vars.p_ref.i(1,s)) = [1;-1];

              Ai(:,obj.vars.G.i(idx_1,i,j,s)) = K;
              Ai(:,obj.vars.F.i(idx_1,i,1,s)) = K;
              bi(1,1) = bi(1,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(2);
              bi(2,1) = bi(2,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(2);
              if obj.Rotate(s)*obj.shape.polygons{i}.center(1) > obj.Rotate(s)*obj.shape.polygons{j}.center(1)
                bi(3,1) = bi(3,1) - obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;
                bi(4,1) = bi(4,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              else 
                bi(3,1) = bi(3,1) - obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
                bi(4,1) = bi(4,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;
              end

              obj = obj.addLinearConstraints(Ai, bi, [], []);
              
              % if the ray passes above the segment
              Ai = sparse(2,obj.nv);
              bi = 2*K*ones(2,1);

              Ai(1:2,obj.vars.p.i(2,idx_1)) = 1;
              Ai(1:2,obj.vars.p_ref.i(2,s)) = 1;
              Ai(1:2,obj.vars.G.i(idx_1,i,j,s)) = K;
              Ai(1:2,obj.vars.F.i(idx_1,i,2,s)) = K;
              bi(1,1) = bi(1,1) - obj.Rotate(s)*obj.shape.polygons{i}.center(2) - dx;
              bi(2,1) = bi(2,1) - obj.Rotate(s)*obj.shape.polygons{j}.center(2) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              % if the ray passes to the right of the segment
              Ai = sparse(2,obj.nv);
              bi = 2*K*ones(2,1);

              Ai(1,obj.vars.p.i(1,idx_1)) = -1;
              Ai(1,obj.vars.p_ref.i(1,s)) = -1;
              Ai(2,obj.vars.p.i(1,idx_1)) = -1;
              Ai(2,obj.vars.p_ref.i(1,s)) = -1;
              Ai(1:2,obj.vars.G.i(idx_1,i,j,s)) = K;
              Ai(1:2,obj.vars.F.i(idx_1,i,3,s)) = K;
              bi(1,1) = bi(1,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              bi(2,1) = bi(2,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              % if the ray passes to the left of the segment
              Ai = sparse(2,obj.nv);
              bi = 2*K*ones(2,1);

              Ai(1,obj.vars.p.i(1,idx_1)) = 1;
              Ai(1,obj.vars.p_ref.i(1,s)) = 1;
              Ai(2,obj.vars.p.i(1,idx_1)) = 1;
              Ai(2,obj.vars.p_ref.i(1,s)) = 1;
              Ai(1:2,obj.vars.G.i(idx_1,i,j,s)) = K;
              Ai(1:2,obj.vars.F.i(idx_1,i,4,s)) = K;
              bi(1,1) = bi(1,1) - obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              bi(2,1) = bi(2,1) - obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              % adds the constraint constrains when the 
              % polygon is intersecting with other pusher
              
              % the ray intersects the segment
              Ai = sparse(4,obj.nv);
              bi = 3*K*ones(4,1);

              Ai(1,obj.vars.p.i(2,idx_1)) = -1;
              Ai(1,obj.vars.p_ref.i(2,s)) = -1;
              Ai(2,obj.vars.p.i(2,idx_2)) = -1;
              Ai(2,obj.vars.p_ref.i(2,s)) = -1;
              Ai(3,obj.vars.p.i(1,idx_1)) = 1;
              Ai(3,obj.vars.p_ref.i(1,s)) = 1;
              Ai(4,obj.vars.p.i(1,idx_2)) = -1;
              Ai(4,obj.vars.p_ref.i(1,s)) = -1;
              Ai(:,obj.vars.F.i(idx_1,i,1,s)) = K;
              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;
              Ai(:,obj.vars.F.i(idx_1,i,5,s)) = K;
              bi(1,1) = bi(1,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(2);
              bi(2,1) = bi(2,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(2);
              bi(3,1) = bi(3,1) - obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              bi(4,1) = bi(4,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              Ai = sparse(4,obj.nv);
              bi = 2*K*ones(4,1);

              Ai(1,obj.vars.p.i(2,idx_1)) = -1;
              Ai(1,obj.vars.p_ref.i(2,s)) = -1;
              Ai(2,obj.vars.p.i(2,idx_2)) = -1;
              Ai(2,obj.vars.p_ref.i(2,s)) = -1;
              Ai(3,obj.vars.p.i(1,idx_1)) = -1;
              Ai(3,obj.vars.p_ref.i(1,s)) = -1;
              Ai(4,obj.vars.p.i(1,idx_2)) = 1;
              Ai(4,obj.vars.p_ref.i(1,s)) = 1;
              Ai(:,obj.vars.F.i(idx_1,i,1,s)) = K;
              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;
              Ai(:,obj.vars.F.i(idx_1,i,5,s)) = -K;
              bi(1,1) = bi(1,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(2);
              bi(2,1) = bi(2,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(2);
              bi(3,1) = bi(3,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              bi(4,1) = bi(4,1) - obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
   
              % the ray passes above the segment
              Ai = sparse(2,obj.nv);
              bi = 2*K*ones(2,1);

              Ai(1,obj.vars.p.i(2,idx_1)) = 1;
              Ai(1,obj.vars.p_ref.i(2,s)) = 1;
              Ai(2,obj.vars.p.i(2,idx_2)) = 1;
              Ai(2,obj.vars.p_ref.i(2,s)) = 1;
              Ai(:,obj.vars.F.i(idx_1,i,2,s)) = K;
              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;
              bi(1,1) = bi(1,1) - obj.Rotate(s)*obj.shape.polygons{i}.center(2) - dx;
              bi(2,1) = bi(2,1) - obj.Rotate(s)*obj.shape.polygons{j}.center(2) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              % if the ray passes to the right of the segment
              Ai = sparse(2,obj.nv);
              bi = 2*K*ones(2,1);

              Ai(1,obj.vars.p.i(1,idx_1)) = -1;
              Ai(2,obj.vars.p.i(1,idx_2)) = -1;
              Ai(:,obj.vars.p_ref.i(1,s)) = -1;
              Ai(:,obj.vars.F.i(idx_1,i,3,s)) = K;
              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;
              bi(1,1) = bi(1,1) + obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              bi(2,1) = bi(2,1) + obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);

              % if the ray passes to the left of the segment
              Ai = sparse(2,obj.nv);
              bi = 2*K*ones(2,1);

              Ai(1,obj.vars.p.i(1,idx_1)) = 1;
              Ai(2,obj.vars.p.i(1,idx_2)) = 1;
              Ai(:,obj.vars.p_ref.i(1,s)) = 1;
              Ai(:,obj.vars.F.i(idx_1,i,4,s)) = K;
              Ai(:,obj.vars.H.i(idx_1,i,j,s)) = K;
              bi(1,1) = bi(1,1) - obj.Rotate(s)*obj.shape.polygons{i}.center(1) - dx;
              bi(2,1) = bi(2,1) - obj.Rotate(s)*obj.shape.polygons{j}.center(1) - dx;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
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
    end
    % end of methods
  end
end