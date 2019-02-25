classdef MixedIntegerCertifiedGraspingProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (MIT MCube Lab)
  properties
    n_e
    shape
    r = 0.1 % min pusher separation
    nT
    dt = 0.1
    mu = 0.5
    theta_range = 3
    theta
    tau
    n_joints
    M
    nl
    nr
  end

  methods
    function obj = MixedIntegerCertifiedGraspingProblem(shape, n_e, nT, th)
      % Constructs the optimization problem and declares the variables for each finger
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
      % @param n_e: number of fingers of the gripper
      % @param nT: number of time-steps

      assert(nargin > 0);
      if nargin < 2; n_e = 3; end
      if nargin < 3; nT = 5; end
      if nargin < 4
        th = [pi/8,pi/12,pi/24,pi/48,0];
      end

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();
      obj.n_e = n_e;
      obj.shape = shape;
      obj.nT = nT;

      obj.theta = zeros(obj.theta_range,obj.nT);
      for t = 1:obj.nT
        obj.theta(:,t) = linspace(-th(t),th(t),obj.theta_range);
      end

      % Pusher locations
      obj = obj.addVariable('p', 'C', [2,obj.n_e,obj.nT+1], -inf, inf);
      obj = obj.addVariable('dp', 'C', [2,obj.n_e,obj.nT+1], -inf, inf);
      obj = obj.addVariable('ddp', 'C', [2,obj.n_e,obj.nT+1], -inf, inf);

      % loop existance
      obj.M = length(obj.shape.polygons);
      obj = obj.addVariable('H','B',[obj.n_e,obj.M,obj.M,obj.nT], 0, 1);
      obj = obj.addVariable('G','B',[obj.n_e,obj.M,obj.M,obj.nT], 0, 1);
      obj = obj.addVariable('weight','C',[obj.n_e,2,obj.shape.nv,obj.nT,obj.theta_range], 0, 1);

      % loop enclosing
      obj = obj.addVariable('F','B',[obj.n_e,obj.M,6], 0, 1);
      obj = obj.addVariable('c','C',[obj.n_e*obj.M-1,1], 0, 1);

      % limit orientation constraints
      obj.nl = length(obj.shape.lines);
      
      obj = obj.addVariable('line', 'B', [obj.nl, obj.n_e, obj.nT, 2], 0, 1);
      obj = obj.addVariable('line_f', 'B', [obj.nl, obj.n_e], 0, 1);

      obj = obj.addVariable('rho_l', 'C', [obj.nl, obj.n_e, obj.nT, 2], 0, 0.99);
      obj = obj.addVariable('rho_l_f', 'C', [obj.nl, obj.n_e], 0, 0.99);

      % continuous boundary variation
      obj = obj.addVariable('weight2','C',[obj.n_e,obj.shape.nv,obj.nT], 0, 1);
      
      % non-penetration
      obj.nr = length(obj.shape.regions);
      obj = obj.addVariable('region', 'B', [obj.nr,obj.n_e,obj.nT], 0, 1);

      % Grasp Observability
      obj = obj.addVariable('M_np', 'B', [1, obj.n_e*(obj.n_e-1)/2], 0, 1);
      obj = obj.addVariable('M_p', 'B', [1, obj.n_e*(obj.n_e-1)/2], 0, 1);

      obj = obj.addVariable('p_int', 'C', [2, obj.n_e*(obj.n_e-1)/2], -inf, inf);
      obj = obj.addVariable('rho_int', 'C', [2, obj.n_e*(obj.n_e-1)/2], -inf, inf);

      obj = obj.addVariable('s_int', 'C', [obj.n_e, obj.n_e*(obj.n_e-1)/2], -inf, inf);
      % obj = obj.addVariable('z_int', 'B', [obj.n_e, obj.n_e*(obj.n_e-1)/2], -inf, inf);
    end

    function obj = addNoCollisionConstraint(obj)
      % constrains the pushers to lie within the complement of the object, 
      % represented through a set of convex polygonal regions using 
      % a logarithmic number of integer variables, as decribed in
      % Modeling Disjunctive Constraints with a Logarithmic Number of Binary
      % Variables and Constraints by J. Vielma and G. Nemhauser, 2011.

      % defines the log2 binary variables
      % ny = ceil(log2(obj.nr));
      % obj = obj.addVariable('y', 'B', [ny, obj.n_e,obj.nT], 0, 1);

      % % defines the gray coding on base 2
      % codes = grayCodes(2,ny);
      % codes = codes(1:obj.nr,:);

      % big-M
      M = 100;

      for t = 1:obj.nT
        % adds the constraints on coding
        % for i = 1:obj.n_e
        %   % for each digit
        %   for j = 1:ny
        %     Ai = sparse(2, obj.nv);
        %     bi = [0;1];
        %     % for each lambda
        %     for k = 1:obj.nr
        %       if codes(k,j) == 1
        %         Ai(1, obj.vars.region.i(k,i,t)) = 1;
        %       elseif codes(k,j) == 0
        %         Ai(2, obj.vars.region.i(k,i,t)) = 1;
        %       else
        %         error('codes need to be 0 or 1');
        %       end  
        %     end
        %     Ai(1, obj.vars.y.i(j,i,t)) = -1;
        %     Ai(2, obj.vars.y.i(j,i,t)) = 1;
        %     obj = obj.addLinearConstraints(Ai, bi, [], []);
        %   end
        % end

        dth = 0;
        theta = obj.theta;

        % regions assignment constraint
        for r = 1:obj.nr
          A = obj.shape.regions{r}.A;
          b = obj.shape.regions{r}.b;
          % constrains for each finger
          for j = 1:obj.n_e
            for th = 1:obj.theta_range
              if t < obj.nT
                R_theta = [cos(theta(th,t)),-sin(theta(th,t));sin(theta(th,t)),cos(theta(th,t))];
                R_theta_1 = [cos(theta(th,t+1)),-sin(theta(th,t+1));sin(theta(th,t+1)),cos(theta(th,t+1))];
              else
                R_theta = [cos(theta(th,t)),-sin(theta(th,t));sin(theta(th,t)),cos(theta(th,t))];
                R_theta_1 = eye(2);
              end
              
              Ai = sparse(size(A, 1), obj.nv);
              bi = zeros(size(A, 1), 1);

              Ai(:, obj.vars.p.i(:,j,t)) = A*(R_theta');
              Ai(:, obj.vars.region.i(r,j,t)) = M;
              bi(:) = b + M;
              obj = obj.addLinearConstraints(Ai, bi, [], []);

              Ai = sparse(size(A, 1), obj.nv);
              bi = zeros(size(A, 1), 1);

              Ai(:, obj.vars.p.i(:,j,t+1)) = A*(R_theta_1');
              Ai(:, obj.vars.region.i(r,j,t)) = M;
              bi(:) = b + M;
              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
        
        % assigns each pusher to one of the region
        for th = 1
          Aeq = sparse(obj.n_e,obj.nv);
          beq = ones(obj.n_e,1);
          for j = 1:obj.n_e
            Aeq(j, obj.vars.region.i(:,j,t)) = 1;
          end
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end
      end
    end

    function obj = addCircleConstraint(obj)
      % Add mixed-integer constraints that require that 
      % the graph formed by the intersection of the 
      % C-space pushers forms a cyclic graph.

      % constrains for all times
      for t = 1:obj.nT
        theta_range = linspace(-pi/(18*t),pi/(18*t),5);
        % constrains the values of the G matrix
        for n = 1:obj.n_e
          for i = 1:obj.M
            for j = 1:obj.M
              Ai = sparse(1,obj.nv);
              bi = obj.shape.G(i,j);
              Ai(1,obj.vars.G.i(n,i,j,t)) = 1;
              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end

        % big-K value
        K = 10;
        % requires that each pair of pushers intersect at least once
        for n = 1:obj.n_e
          Aeq = sparse(1,obj.nv);
          beq = 1;
          Aeq(1,obj.vars.H.i(n,:,:,t)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end

        % constrains the intersection between pushers
        for n = 1:obj.n_e
          idx_1 = n;
          idx_2 = mod(n+1,obj.n_e);
          if(idx_2 == 0); idx_2 = obj.n_e; end;

          % intersection constraint
          for i = 1:obj.M
            for j = 1:obj.M
              for th = 1:obj.theta_range
                % intersection of the pushers
                Ai = sparse(4,obj.nv);
                bi = K*ones(4,1);

                R = [cos(obj.theta(th,t)), sin(obj.theta(th,t)); -sin(obj.theta(th,t)), cos(obj.theta(th,t))];

                Ai(:,obj.vars.H.i(idx_1,i,j,t)) = K;

                Ai(1:2,obj.vars.p.i(:,idx_1,t)) = eye(2);
                Ai(1:2,obj.vars.p.i(:,idx_2,t)) = -eye(2);

                Ai(3:4,obj.vars.p.i(:,idx_1,t)) = -eye(2);
                Ai(3:4,obj.vars.p.i(:,idx_2,t)) = eye(2);

                for k = 1:obj.shape.polygons{i}.nv
                  Ai(1:2,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,t,th)) = -R*obj.shape.polygons{i}.v(:,k);
                  Ai(3:4,obj.vars.weight.i(idx_1,1,obj.shape.polygons{i}.iv+k,t,th)) = R*obj.shape.polygons{i}.v(:,k);           
                end

                for l = 1:obj.shape.polygons{j}.nv
                  Ai(1:2,obj.vars.weight.i(idx_1,2,obj.shape.polygons{j}.iv+l,t,th)) = R*obj.shape.polygons{j}.v(:,l);           
                  Ai(3:4,obj.vars.weight.i(idx_1,2,obj.shape.polygons{j}.iv+l,t,th)) = -R*obj.shape.polygons{j}.v(:,l); 
                end

                obj = obj.addLinearConstraints(Ai, bi, [], []);
              end
            end
          end

          % weights must add to 1 if intersecting, 0 otherwise
          for i = 1:obj.M
            for th = 1:obj.theta_range
              Aeq = sparse(2,obj.nv);
              beq = zeros(2,1);

              range_1 = (obj.shape.polygons{i}.iv+1):(obj.shape.polygons{i}.iv + obj.shape.polygons{i}.nv);

              Aeq(1,obj.vars.H.i(idx_1,i,:,t)) = -1;
              Aeq(1,obj.vars.weight.i(idx_1,1,range_1,t,th)) = 1;

              Aeq(2,obj.vars.H.i(idx_1,:,i,t)) = -1;
              Aeq(2,obj.vars.weight.i(idx_1,2,range_1,t,th)) = 1;
              
              obj = obj.addLinearConstraints([], [], Aeq, beq);
            end
          end
        end

        % constrains the existance of a closed circle
        for n = 1:obj.n_e
          idx_1 = n;
          idx_2 = mod(n+1,obj.n_e);
          if(idx_2 == 0); idx_2 = obj.n_e; end;

          for i = 1:obj.M
            for j = 1:obj.M
              Ai = sparse(2,obj.nv);
              bi = [K+1,K-1]';

              % for two pushers the graph must go in and out of the 
              % C-space pushers
              if(obj.n_e > 2)
                l_range = 1:obj.M;
              else
                l_range = [1:i-1,i+1:obj.M];
              end

              Ai(1,obj.vars.H.i(idx_2,j,l_range,t)) = 1;
              Ai(1,obj.vars.G.i(idx_2,j,:,t)) = 1;

              Ai(2,obj.vars.H.i(idx_2,j,l_range,t)) = -1;
              Ai(2,obj.vars.G.i(idx_2,j,:,t)) = -1;
              Ai(:,obj.vars.H.i(idx_1,i,j,t)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end

          for p = 1:obj.M
            for q = 1:obj.M
              Ai = sparse(2,obj.nv);
              bi = [K+1,K-1]';

              r_range = [1:p-1,p+1:obj.M];

              Ai(1,obj.vars.H.i(idx_2,q,:,t)) = 1;
              Ai(1,obj.vars.G.i(idx_2,q,r_range,t)) = 1;

              Ai(2,obj.vars.H.i(idx_2,q,:,t)) = -1;
              Ai(2,obj.vars.G.i(idx_2,q,r_range,t)) = -1;
              Ai(:,obj.vars.G.i(idx_2,p,q,t)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
      end
    end

    function obj = addEnclosingConstraint(obj)
      % constrains the origin to lie inside the circle formed by the 
      % pushers in the C-space of the object by using
      
      % big K
      K = 100;

      for i = 1:obj.n_e
        idx_1 = i;
        idx_2 = mod(i+1,obj.n_e);
        if(idx_2 == 0); idx_2 = obj.n_e; end
        for j = 1:obj.M
          % requires that each F adds to one
          Aeq = sparse(1,obj.nv);
          beq = 1;
          Aeq(1,obj.vars.F.i(i,j,[1:4,6])) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);

          Ai = sparse(1,obj.nv);
          bi = 0;
          Ai(1,obj.vars.F.i(i,j,1:4)) = 1;
          Ai(1,obj.vars.H.i(i,j,:,1)) = -1;
          Ai(1,obj.vars.G.i(i,j,:,1)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Aeq = sparse(1,obj.nv);
          beq = 1;
          Aeq(1,obj.vars.F.i(i,j,6)) = 1;
          Aeq(1,obj.vars.H.i(i,j,:,1)) = 1;
          Aeq(1,obj.vars.G.i(i,j,:,1)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end
      end

      % performs XOR sequentially
      n = 1;

      for i = 1:obj.n_e
        for j = 1:obj.M
          if j == 1 && i == 1
            obj = obj.addXORConstraint(obj.vars.c.i(1,1),obj.vars.F.i(1,1,1),obj.vars.F.i(1,2,1));
          else 
            if j < obj.M
              obj = obj.addXORConstraint(obj.vars.c.i(n),obj.vars.c.i(n-1,1),obj.vars.F.i(i,j+1,1));
            else
              if i < obj.n_e
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
      dx = 0.001;

      for n = 1:obj.n_e
        idx_1 = n;
        idx_2 = mod(n+1,obj.n_e);
        if(idx_2 == 0); idx_2 = obj.n_e; end
        for i = 1:obj.M
          for j = 1:obj.M
            % adds the constraint when the
            % polygon is not intersecting
            Ai = sparse(4,obj.nv);
            bi = 2*K*ones(4,1);

            Ai(1:2,obj.vars.p.i(2,idx_1,1)) = -1;
            Ai(3:4,obj.vars.p.i(1,idx_1,1)) = [1;-1];
            Ai(:,obj.vars.G.i(idx_1,i,j,1)) = K;
            Ai(:,obj.vars.F.i(idx_1,i,1)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(2);
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(2);
            if obj.shape.polygons{i}.center(1) > obj.shape.polygons{j}.center(1)
              bi(3,1) = bi(3,1) + obj.shape.polygons{j}.center(1) - dx;
              bi(4,1) = bi(4,1) - obj.shape.polygons{i}.center(1) - dx;
            else 
              bi(3,1) = bi(3,1) + obj.shape.polygons{i}.center(1) - dx;
              bi(4,1) = bi(4,1) - obj.shape.polygons{j}.center(1) - dx;
            end

            obj = obj.addLinearConstraints(Ai, bi, [], []);
            
            % if the ray passes above the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1:2,obj.vars.p.i(2,idx_1,1)) = 1;
            Ai(1:2,obj.vars.G.i(idx_1,i,j,1)) = K;
            Ai(1:2,obj.vars.F.i(idx_1,i,2)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(2) - dx;
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(2) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the right of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1,1)) = -1;
            Ai(2,obj.vars.p.i(1,idx_1,1)) = -1;
            Ai(1:2,obj.vars.G.i(idx_1,i,j,1)) = K;
            Ai(1:2,obj.vars.F.i(idx_1,i,3)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the left of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1,1)) = 1;
            Ai(2,obj.vars.p.i(1,idx_1,1)) = 1;
            Ai(1:2,obj.vars.G.i(idx_1,i,j,1)) = K;
            Ai(1:2,obj.vars.F.i(idx_1,i,4)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % adds the constraint constrains when the 
            % polygon is intersecting with other pusher
            
            % the ray intersects the segment
            Ai = sparse(4,obj.nv);
            bi = 3*K*ones(4,1);

            Ai(1,obj.vars.p.i(2,idx_1,1)) = -1;
            Ai(2,obj.vars.p.i(2,idx_2,1)) = -1;
            Ai(3,obj.vars.p.i(1,idx_1,1)) = 1;
            Ai(4,obj.vars.p.i(1,idx_2,1)) = -1;
            Ai(:,obj.vars.F.i(idx_1,i,1)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j,1)) = K;
            Ai(:,obj.vars.F.i(idx_1,i,5)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(2);
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(2);
            bi(3,1) = bi(3,1) + obj.shape.polygons{i}.center(1) - dx;
            bi(4,1) = bi(4,1) - obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            Ai = sparse(4,obj.nv);
            bi = 2*K*ones(4,1);

            Ai(1,obj.vars.p.i(2,idx_1,1)) = -1;
            Ai(2,obj.vars.p.i(2,idx_2,1)) = -1;
            Ai(3,obj.vars.p.i(1,idx_1,1)) = -1;
            Ai(4,obj.vars.p.i(1,idx_2,1)) = 1;
            Ai(:,obj.vars.F.i(idx_1,i,1)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j,1)) = K;
            Ai(:,obj.vars.F.i(idx_1,i,5)) = -K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(2);
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(2);
            bi(3,1) = bi(3,1) - obj.shape.polygons{i}.center(1) - dx;
            bi(4,1) = bi(4,1) + obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
 
            % the ray passes above the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(2,idx_1,1)) = 1;
            Ai(2,obj.vars.p.i(2,idx_2,1)) = 1;
            Ai(:,obj.vars.F.i(idx_1,i,2)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j,1)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(2) - dx;
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(2) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the right of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1,1)) = -1;
            Ai(2,obj.vars.p.i(1,idx_2,1)) = -1;
            Ai(:,obj.vars.F.i(idx_1,i,3)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j,1)) = K;
            bi(1,1) = bi(1,1) - obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) - obj.shape.polygons{j}.center(1) - dx;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % if the ray passes to the left of the segment
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(1,obj.vars.p.i(1,idx_1,1)) = 1;
            Ai(2,obj.vars.p.i(1,idx_2,1)) = 1;
            Ai(:,obj.vars.F.i(idx_1,i,4)) = K;
            Ai(:,obj.vars.H.i(idx_1,i,j,1)) = K;
            bi(1,1) = bi(1,1) + obj.shape.polygons{i}.center(1) - dx;
            bi(2,1) = bi(2,1) + obj.shape.polygons{j}.center(1) - dx;

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
      % Performs a linear integration over the dynamics
      % of the hand motion. In this case, Direct Transcription
      % is used to linearly interpolate between timesteps.

      for i = 1:obj.n_e
        % direct transcription of the motion
        for j = 1:obj.nT
          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, obj.vars.p.i(:,i,j+1)) = eye(2);
          Aeq(:, obj.vars.p.i(:,i,j)) = -eye(2);
          Aeq(:, obj.vars.dp.i(:,i,j)) = -obj.dt*eye(2);

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, obj.vars.dp.i(:,i,j+1)) = eye(2);
          Aeq(:, obj.vars.dp.i(:,i,j)) = -eye(2);
          Aeq(:, obj.vars.ddp.i(:,i,j)) = -obj.dt*eye(2);

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          % upper bounds
          Ai = sparse(4, obj.nv);
          bi = 30*ones(4,1);

          Ai(1:2,obj.vars.ddp.i(:,i,j)) = eye(2);
          Ai(3:4,obj.vars.ddp.i(:,i,j)) = -eye(2);

          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end
    end

    function obj = addCspaceContractionConstraints(obj)
      % Constrains the conditions required for the object to reach a limit orientation,
      % depending on the number of pushers. A limit orientation is that for which the
      % object is fully impobilized by the pushers, with a zero-area free-space component.

      % big K 
      K = 100;
      nl = obj.nl;

      % line assignment
      for t = 1:obj.nT
        R_1 = [cos(obj.theta(1,t)),-sin(obj.theta(1,t)); sin(obj.theta(1,t)),cos(obj.theta(1,t))];
        R_2 = [cos(obj.theta(end,t)),-sin(obj.theta(end,t)); sin(obj.theta(end,t)),cos(obj.theta(end,t))];

        for i = 1:obj.n_e
          % each pusher can be at most in one line
          Ai = sparse(2, obj.nv);
          bi = ones(2, 1);
          Ai(1,obj.vars.line.i(:,i,t,1)) = 1;
          Ai(2,obj.vars.line.i(:,i,t,2)) = 1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          for l = 1:nl
            % multipliers are only activated when in contact
            Ai = sparse(2, obj.nv);
            bi = zeros(2,1);

            Ai(1,obj.vars.line.i(l,i,t,1)) = -1;
            Ai(1,obj.vars.rho_l.i(l,i,t,1)) = 1;

            Ai(2,obj.vars.line.i(l,i,t,2)) = -1;
            Ai(2,obj.vars.rho_l.i(l,i,t,2)) = 1;
            
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % cannot make contact with borders
            Ai = sparse(2, obj.nv);
            bi = zeros(2,1);

            Ai(1,obj.vars.line.i(l,i,t,1)) = 0.001;
            Ai(1,obj.vars.rho_l.i(l,i,t,1)) = -1;

            Ai(2,obj.vars.line.i(l,i,t,2)) = 0.001;
            Ai(2,obj.vars.rho_l.i(l,i,t,2)) = -1;
            
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            % line assignment via big-M formulation
            Ai = sparse(8,obj.nv);
            bi = K*ones(8,1);

            v11 = R_1*obj.shape.lines{l}.v1;
            v21 = R_1*obj.shape.lines{l}.v2;

            Ai(1:4,obj.vars.line.i(l,i,t,1)) = K;

            Ai(1:2,obj.vars.p.i(:,i,t)) = eye(2);
            Ai(1:2,obj.vars.rho_l.i(l,i,t,1)) = v21-v11;
            bi(1:2) = bi(1:2) + v21;

            Ai(3:4,obj.vars.p.i(:,i,t)) = -eye(2);
            Ai(3:4,obj.vars.rho_l.i(l,i,t,1)) = v11-v21;
            bi(3:4) = bi(3:4) - v21;

            v12 = R_2*obj.shape.lines{l}.v1;
            v22 = R_2*obj.shape.lines{l}.v2;

            Ai(5:8,obj.vars.line.i(l,i,t,2)) = K;

            Ai(5:6,obj.vars.p.i(:,i,t)) = eye(2);
            Ai(5:6,obj.vars.rho_l.i(l,i,t,2)) = v22-v12;
            bi(5:6) = bi(5:6) + v22;

            Ai(7:8,obj.vars.p.i(:,i,t)) = -eye(2);
            Ai(7:8,obj.vars.rho_l.i(l,i,t,2)) = v12-v22;
            bi(7:8) = bi(7:8) - v22;

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end

      for i = 1:obj.n_e
        for l = 1:obj.nl
          Aeq = sparse(1, obj.nv);
          beq = zeros(1,1);

          Aeq(1,obj.vars.line.i(l,i,obj.nT,1)) = 1;
          Aeq(1,obj.vars.line.i(l,i,obj.nT,2)) = -1;
          
          % obj = obj.addLinearConstraints([], [], Aeq, beq);
        end
      end

      % depending on the number of fingers determines what is a critical orientation
      if obj.n_e == 2
        for t = 1:obj.nT
          % both fingers must be in contact to have a critical orientation
          Aeq = sparse(2, obj.nv);
          beq = ones(2, 1);

          Aeq(1,obj.vars.line.i(:,1,t,1)) = 1;
          Aeq(2,obj.vars.line.i(:,2,t,2)) = 1;

          obj = obj.addLinearConstraints([], [], Aeq, beq);   

          % only contact with opposite faces counts
          for i = 1:obj.n_e
            idx_1 = i;
            idx_2 = mod(i+1,obj.n_e);
            if(idx_2 == 0); idx_2 = obj.n_e; end;
            
            for l = 1:nl
              Ai = sparse(2,obj.nv);
              bi = zeros(2,1);

              Ai(1,obj.vars.line.i(l,idx_1,t,1)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.opp,idx_2,t,1)) = -1;

              Ai(2,obj.vars.line.i(l,idx_1,t,2)) = 1;
              Ai(2,obj.vars.line.i(obj.shape.lines{l}.opp,idx_2,t,2)) = -1;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end 
        end
      else
        % there are three types of critical orientations
        obj = obj.addVariable('crit_type', 'B', [3,obj.nT,2], 0, 1);

        for t = 1:obj.nT
          Aeq = sparse(2, obj.nv);
          beq = ones(2,1);          

          Aeq(1,obj.vars.crit_type.i(:,t,1)) = 1;
          Aeq(2,obj.vars.crit_type.i(:,t,2)) = 1;

          obj = obj.addLinearConstraints([], [], Aeq, beq);

          Ai = sparse(2, obj.nv);
          bi = -2*ones(2,1);

          Ai(1,obj.vars.line.i(:,:,t,1)) = -1;
          Ai(2,obj.vars.line.i(:,:,t,2)) = -1;

          obj = obj.addLinearConstraints(Ai, bi, [], []);

          % For 2 pushers in contact:
          % only contact with opposite faces counts
          for i = 1:obj.n_e
            idx_1 = i;
            idx_comp = [1:idx_1-1,idx_1+1:obj.n_e];

            for l = 1:nl
              Ai = sparse(2, obj.nv);
              bi = K*ones(2,1);

              Ai(1,obj.vars.line.i(l,idx_1,t,1)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.opp,:,t,1)) = -1;
              Ai(1,obj.vars.crit_type.i(1,t,1)) = K;

              Ai(2,obj.vars.line.i(l,idx_1,t,2)) = 1;
              Ai(2,obj.vars.line.i(obj.shape.lines{l}.opp,:,t,2)) = -1;
              Ai(2,obj.vars.crit_type.i(1,t,2)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end 

          % For 3 pushers or more in contact:
          % only contact with non-codirectional facets counts
          for i = 1:obj.n_e
            for l = 1:nl
              Ai = sparse(2, obj.nv);
              bi = K*ones(2,1);

              Ai(1,obj.vars.line.i(l,i,t,1)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.non_cod,:,t,1)) = -1/2;
              Ai(1,obj.vars.crit_type.i(2,t,1)) = K;

              Ai(2,obj.vars.line.i(l,i,t,2)) = 1;
              Ai(2,obj.vars.line.i(obj.shape.lines{l}.non_cod,:,t,2)) = -1/2;
              Ai(2,obj.vars.crit_type.i(2,t,2)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end

          % For 3 pushers or more in contact:
          % only contact with non-equal facets counts
          for i = 1:obj.n_e
            for l = 1:nl
              Ai = sparse(2, obj.nv);
              bi = K*ones(2,1);

              Ai(1,obj.vars.line.i(l,i,t,1)) = 1;
              Ai(1,obj.vars.line.i(obj.shape.lines{l}.non_eq,:,t,1)) = -1/3;
              Ai(1,obj.vars.crit_type.i(3,t,1)) = K;

              Ai(2,obj.vars.line.i(l,i,t,2)) = 1;
              Ai(2,obj.vars.line.i(obj.shape.lines{l}.non_eq,:,t,2)) = -1/3;
              Ai(2,obj.vars.crit_type.i(3,t,2)) = K;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
        end
      end
    end

    function obj = addTerminalConstraints(obj,mode)
      if nargin < 2
        mode = 1:2;
      end

      % constrains the final time-step to be a form-closure grasp
      K = 100;
      nl = obj.nl;
      t = obj.nT+1;

      Ai = sparse(1, obj.nv);
      bi = -2*ones(1,1);
      Ai(1,obj.vars.line_f.i(:,:)) = -1;
      obj = obj.addLinearConstraints(Ai, bi, [], []);

      % Ai = sparse(4, obj.nv);
      % bi = -ones(4,1);
      % Ai(1,obj.vars.line_f.i(5,4)) = -1;
      % Ai(2,obj.vars.line_f.i(6,3)) = -1;
      % Ai(3,obj.vars.line_f.i(11,1)) = -1;
      % Ai(4,obj.vars.line_f.i(11,2)) = -1;
      % obj = obj.addLinearConstraints(Ai, bi, [], []);

      obj = obj.addVariable('assigned', 'C', [obj.nl,1], 0, 1);

      for i = 1:obj.nl
        for e = 1:obj.n_e
          Ai = sparse(2,obj.nv);
          bi = [K+1;K-1];

          Ai(1,obj.vars.assigned.i(i,1)) = 1;
          Ai(2,obj.vars.assigned.i(i,1)) = -1;
          Ai(:,obj.vars.line_f.i(i,e)) = K;

          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end

        Ai = sparse(2,obj.nv);
        bi = zeros(2,1);

        Ai(1,obj.vars.assigned.i(i,1)) = 1;
        Ai(2,obj.vars.assigned.i(i,1)) = -1;
        Ai(:,obj.vars.line_f.i(i,:)) = -K;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end

      % line assignment for the final case
      for i = 1:obj.n_e
        % one assignment per finger
        Ai = sparse(1, obj.nv);
        bi = ones(1, 1);
        Ai(1,obj.vars.line_f.i(:,i)) = 1;
        obj = obj.addLinearConstraints(Ai, bi, [], []);

        for l = 1:nl
          % line assignment
          Ai = sparse(2, obj.nv);
          bi = zeros(2,1);

          Ai(1,obj.vars.line_f.i(l,i)) = -1;
          Ai(1,obj.vars.rho_l_f.i(l,i)) = 1;

          Ai(2,obj.vars.line_f.i(l,i)) = 0.001;
          Ai(2,obj.vars.rho_l_f.i(l,i)) = -1;

          obj = obj.addLinearConstraints(Ai, bi, [], []);

          % line assignment via big-M formulation
          Ai = sparse(4,obj.nv);
          bi = K*ones(4,1);

          v1 = obj.shape.lines{l}.v1;
          v2 = obj.shape.lines{l}.v2;

          Ai(1:4,obj.vars.line_f.i(l,i)) = K;

          Ai(1:2,obj.vars.p.i(:,i,t)) = eye(2);
          Ai(1:2,obj.vars.rho_l_f.i(l,i)) = v2-v1;
          bi(1:2) = bi(1:2) + v2;

          Ai(3:4,obj.vars.p.i(:,i,t)) = -eye(2);
          Ai(3:4,obj.vars.rho_l_f.i(l,i)) = v1-v2;
          bi(3:4) = bi(3:4) - v2;

          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end

      % types of limit orientation that lead to form closure
      obj = obj.addVariable('crit_type_f', 'B', [3,1], 0, 1);

      Aeq = sparse(1, obj.nv);
      beq = ones(1,1);          
      Aeq(1,obj.vars.crit_type_f.i(mode,1)) = 1;
      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % determines the number of contacts active
      for l = 1:nl
        Ai = sparse(1, obj.nv);
        bi = K*ones(1,1);

        Ai(1,obj.vars.assigned.i(l,1)) = 1;
        Ai(1,obj.vars.assigned.i(obj.shape.lines{l}.opp,1)) = -1;
        Ai(1,obj.vars.crit_type_f.i(3,1)) = K;

        obj = obj.addLinearConstraints(Ai, bi, [], []);

        Ai = sparse(1, obj.nv);
        bi = K*ones(1,1);

        Ai(1,obj.vars.assigned.i(l,1)) = 1;
        Ai(1,obj.vars.assigned.i(obj.shape.lines{l}.non_cod,1)) = -1/2;
        Ai(1,obj.vars.crit_type_f.i(1,1)) = K;

        obj = obj.addLinearConstraints(Ai, bi, [], []);

        Ai = sparse(1, obj.nv);
        bi = K*ones(1,1);

        Ai(1,obj.vars.assigned.i(l,1)) = 1;
        Ai(1,obj.vars.assigned.i(obj.shape.lines{l}.non_eq,1)) = -1/3;
        Ai(1,obj.vars.crit_type_f.i(2,1)) = K;

        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end
    end

    function obj = addObservabilityConstraints(obj)
      % constrains the existance of possitive linear independance
      % between contact normals with respect to each contact location

      % checks for pairs of non-parallel normals
      idx = 1;
      K = 100;
      for i = 1:obj.n_e-1
        for j = i+1:obj.n_e

          Ai = sparse(4,obj.nv);
          bi = zeros(4,1);
          Ai(:,obj.vars.M_np.i(1,idx)) = -K;
          Ai(:,obj.vars.M_p.i(1,idx)) = -K;
          Ai(1:2,obj.vars.p_int.i(:,idx)) = eye(2);
          Ai(3:4,obj.vars.p_int.i(:,idx)) = -eye(2);
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Ai = sparse(1,obj.nv);
          bi = zeros(1,1);
          Ai(1,obj.vars.M_np.i(1,idx)) = 1;
          Ai(1,obj.vars.line_f.i(:,i)) = -K;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          for l = 1:obj.nl
            Ai = sparse(1,obj.nv);
            bi = 0.5;

            Ai(1,obj.vars.M_np.i(1,idx)) = -1;
            Ai(1,obj.vars.line_f.i(l,i)) = 1/2;
            Ai(1,obj.vars.line_f.i(obj.shape.lines{l}.non_cod,j)) = 1/2;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            for l_2 = 1:obj.nl
              Ai = sparse(8,obj.nv);
              bi = 3*K*ones(8,1);

              Ai(:,obj.vars.M_np.i(1,idx)) = K;
              Ai(:,obj.vars.line_f.i(l,i)) = K;
              Ai(:,obj.vars.line_f.i(l_2,j)) = K;

              Ai(1:2,obj.vars.p_int.i(:,idx)) = -eye(2);
              Ai(1:2,obj.vars.p.i(:,i,end)) = eye(2);
              Ai(1:2,obj.vars.rho_int.i(1,idx)) = obj.shape.lines{l}.normal;

              Ai(5:6,obj.vars.p_int.i(:,idx)) = eye(2);
              Ai(5:6,obj.vars.p.i(:,i,end)) = -eye(2);
              Ai(5:6,obj.vars.rho_int.i(1,idx)) = -obj.shape.lines{l}.normal;

              Ai(3:4,obj.vars.p_int.i(:,idx)) = -eye(2);
              Ai(3:4,obj.vars.p.i(:,j,end)) = eye(2);
              Ai(3:4,obj.vars.rho_int.i(2,idx)) = obj.shape.lines{l_2}.normal;

              Ai(7:8,obj.vars.p_int.i(:,idx)) = eye(2);
              Ai(7:8,obj.vars.p.i(:,j,end)) = -eye(2);
              Ai(7:8,obj.vars.rho_int.i(2,idx)) = -obj.shape.lines{l_2}.normal;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
          idx = idx + 1;
        end
      end
      
      % checks for pairs of opposite normals 
      idx = 1;
      for i = 1:obj.n_e-1
        for j = i+1:obj.n_e
          Ai = sparse(1,obj.nv);
          bi = zeros(1,1);
          Ai(1,obj.vars.M_p.i(1,idx)) = 1;
          Ai(1,obj.vars.line_f.i(:,i)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          for l = 1:obj.nl
            Ai = sparse(1,obj.nv);
            bi = zeros(1,1);

            Ai(1,obj.vars.M_p.i(1,idx)) = 1;
            Ai(1,obj.vars.line_f.i(l,i)) = -1/2;
            Ai(1,obj.vars.line_f.i(obj.shape.lines{l}.parallel_fac,j)) = -1/2;

            obj = obj.addLinearConstraints(Ai, bi, [], []);

            for l_2 = 1:obj.nl
              Ai = sparse(8,obj.nv);
              bi = 3*K*ones(8,1);

              Ai(:,obj.vars.M_p.i(1,idx)) = K;
              Ai(:,obj.vars.line_f.i(l,i)) = K;
              Ai(:,obj.vars.line_f.i(l_2,j)) = K;

              Ai(1:2,obj.vars.p_int.i(:,idx)) = -eye(2);
              Ai(1:2,obj.vars.p.i(:,i,end)) = eye(2);
              Ai(1:2,obj.vars.rho_int.i(1,idx)) = obj.shape.lines{l}.normal;

              Ai(5:6,obj.vars.p_int.i(:,idx)) = eye(2);
              Ai(5:6,obj.vars.p.i(:,i,end)) = -eye(2);
              Ai(5:6,obj.vars.rho_int.i(1,idx)) = -obj.shape.lines{l}.normal;

              Ai(3:4,obj.vars.p_int.i(:,idx)) = -eye(2);
              Ai(3:4,obj.vars.p.i(:,j,end)) = eye(2);
              Ai(3:4,obj.vars.rho_int.i(2,idx)) = obj.shape.lines{l_2}.normal;

              Ai(7:8,obj.vars.p_int.i(:,idx)) = eye(2);
              Ai(7:8,obj.vars.p.i(:,j,end)) = -eye(2);
              Ai(7:8,obj.vars.rho_int.i(2,idx)) = -obj.shape.lines{l_2}.normal;

              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
          end
          idx = idx + 1;
        end
      end

      % computes if there is a non-congruent point
      for i = 1:obj.n_e
        for k = 1:(obj.n_e-1)*obj.n_e/2
          Ai = sparse(2,obj.nv);
          bi = zeros(2,1);

          Ai(:,obj.vars.M_np.i(1,k)) = -K;
          Ai(:,obj.vars.M_p.i(1,k)) = -K;
          Ai(1,obj.vars.s_int.i(i,k)) = 1;
          Ai(2,obj.vars.s_int.i(i,k)) = -1;

          obj = obj.addLinearConstraints(Ai, bi, [], []);
          
          for l = 1:obj.nl
            Ai = sparse(2,obj.nv);
            bi = 2*K*ones(2,1);

            Ai(:,obj.vars.line_f.i(l,i)) = K;
            Ai(:,obj.vars.M_np.i(1,k)) = K;
            Ai(:,obj.vars.M_p.i(1,k)) = K;

            Ai(1,obj.vars.s_int.i(i,k)) = 1;
            Ai(1,obj.vars.p_int.i(:,k)) = -[obj.shape.lines{l}.normal(2),-obj.shape.lines{l}.normal(1)];
            Ai(1,obj.vars.p.i(:,i,end)) = [obj.shape.lines{l}.normal(2),-obj.shape.lines{l}.normal(1)];

            Ai(2,obj.vars.s_int.i(i,k)) = -1;
            Ai(2,obj.vars.p_int.i(:,k)) = [obj.shape.lines{l}.normal(2),-obj.shape.lines{l}.normal(1)];
            Ai(2,obj.vars.p.i(:,i,end)) = -[obj.shape.lines{l}.normal(2),-obj.shape.lines{l}.normal(1)];

            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end

      % constrains that at least three of the contact normals are non-coincident
      Ai = sparse(1,obj.nv);
      bi = -0.001;
      Ai(1,obj.vars.s_int.i(:,:)) = 1;
      obj = obj.addLinearConstraints(Ai, bi, [], []);
    end

    function obj = addCostFunction(obj)
      % Minimizes the separation between pushers in the final configuration
      for t = 1:obj.nT
        for j = 1:obj.n_e  
          Qi = sparse(obj.nv,obj.nv);
          Qi(obj.vars.ddp.i(:,j,t),obj.vars.ddp.i(:,j,t)) = eye(2);
          obj = obj.addCost(Qi,[],[]);
        end

        for j = 1:obj.n_e-1
          Qi = sparse(obj.nv,obj.nv);
          Qi(obj.vars.p.i(:,j,t),obj.vars.p.i(:,j,t)) = eye(2);
          Qi(obj.vars.p.i(1,j,t),obj.vars.p.i(1,j+1,t)) = -2;
          Qi(obj.vars.p.i(2,j,t),obj.vars.p.i(2,j+1,t)) = -2;
          Qi(obj.vars.p.i(:,j+1,t),obj.vars.p.i(:,j+1,t)) = eye(2);
          % obj = obj.addCost(Qi,[],[]);
        end
      end
    end

    function obj = addYumiKinConstraints(obj)
      % Adds the kinematics of Yumi as a constraint
      assert(obj.n_e==4);

      obj = obj.addVariable('p1_hand', 'C', [2,obj.nT+1], -inf, inf);
      obj = obj.addVariable('p2_hand', 'C', [2,obj.nT+1], -inf, inf);

      Ai = sparse(2, obj.nv);
      bi = [-0.6,-0.6]';

      Ai(1, obj.vars.p.i(1,1,1)) = 1;
      Ai(1, obj.vars.p.i(1,1,end)) = -1;

      Ai(2, obj.vars.p.i(1,4,1)) = 1;
      Ai(2, obj.vars.p.i(1,4,end)) = -1;

      % obj = obj.addLinearConstraints(Ai,bi,[],[]);

      % relative position between hands must remain constant
      for i = 1:obj.nT+1
        % for hand 1 (left)
        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);

        Aeq(1, obj.vars.p.i(2,2,i)) = 1;
        Aeq(1, obj.vars.p.i(2,1,i)) = -1;

        Aeq(2, obj.vars.p1_hand.i(2,i)) = 1;
        Aeq(2, obj.vars.p.i(2,2,i)) = -1;

        Aeq(3, obj.vars.p1_hand.i(1,i)) = 2;
        Aeq(3, obj.vars.p.i(1,1,i)) = -1;
        Aeq(3, obj.vars.p.i(1,2,i)) = -1;

        Ai = sparse(3, obj.nv);
        bi = zeros(3, 1);

        Ai(1, obj.vars.p.i(1,1,i)) = 1;
        Ai(1, obj.vars.p.i(1,2,i)) = -1;

        % maximum finger separation is 5 cm
        Ai(2, obj.vars.p.i(1,1,i)) = -1;
        Ai(2, obj.vars.p.i(1,2,i)) = 1;
        bi(2,1) = 5;

        Ai(3, obj.vars.p.i(1,1,i)) = 1;
        Ai(3, obj.vars.p.i(1,2,i)) = -1;
        bi(3,1) = -2;

        obj = obj.addLinearConstraints(Ai,bi,Aeq,beq);

        % for hand 2 (right)
        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);

        Aeq(1, obj.vars.p2_hand.i(2,i)) = 1;
        Aeq(1, obj.vars.p.i(2,3,i)) = -1;

        Aeq(2, obj.vars.p2_hand.i(2,i)) = 1;
        Aeq(2, obj.vars.p.i(2,4,i)) = -1;

        Aeq(3, obj.vars.p2_hand.i(1,i)) = 2;
        Aeq(3, obj.vars.p.i(1,3,i)) = -1;
        Aeq(3, obj.vars.p.i(1,4,i)) = -1;

        Ai = sparse(3, obj.nv);
        bi = zeros(3, 1);

        Ai(1, obj.vars.p.i(1,4,i)) = 1;
        Ai(1, obj.vars.p.i(1,3,i)) = -1;

        % maximum finger separation is 5 cm
        Ai(2, obj.vars.p.i(1,4,i)) = -1;
        Ai(2, obj.vars.p.i(1,3,i)) = 1;
        bi(2,1) = 5;

        Ai(3, obj.vars.p.i(1,4,i)) = 1;
        Ai(3, obj.vars.p.i(1,3,i)) = -1;
        bi(3,1) = -2;

        obj = obj.addLinearConstraints(Ai,bi,Aeq,beq);

        % vetical separation between hands
        Ai = sparse(2, obj.nv);
        bi = [-1.5;4.5];

        Ai(1, obj.vars.p1_hand.i(2,i)) = -1;
        Ai(1, obj.vars.p2_hand.i(2,i)) = 1;

        Ai(2, obj.vars.p1_hand.i(2,i)) = 1;
        Ai(2, obj.vars.p2_hand.i(2,i)) = -1;

        obj = obj.addLinearConstraints(Ai,bi,[],[]);
      end
    end

    % end of methods
  end
end
%EOF