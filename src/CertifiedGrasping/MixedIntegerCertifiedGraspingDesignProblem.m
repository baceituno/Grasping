classdef MixedIntegerCertifiedGraspingDesignProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (MIT MCube Lab)
  properties
    n_e
    shape
    p
    p_hand
    r = 0.1 % min pusher separation
    nT
    dt = 0.1
    theta_range
    theta
    M
    nl
    nr
    mult
  end

  methods
    function obj = MixedIntegerCertifiedGraspingDesignProblem(previous, mult)
      % Constructs the optimization problem and declares the variables for each finger
      % @param previous: previously solved MIP for Certified Grasping design
      % @param mult: number of interpolation points in between contacts

      assert(nargin > 0);
      if nargin < 2
        mult = 3; 
      end

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();
      obj.n_e = previous.n_e;
      obj.shape = previous.shape;
      obj.nT = previous.nT;
      obj.mult = mult;
      obj.p = previous.vars.p.value;
      obj.p_hand = previous.vars.p_hand.value;
      obj.region_old = previous.vars.region.value;

      % Pusher locations
      obj = obj.addVariable('p_int', 'C', [2,(obj.n_e)*mult,obj.nT], -inf, inf);

      % non-penetration
      obj.nr = length(obj.shape.regions);
      obj = obj.addVariable('region', 'B', [obj.nr,obj.n_e*mult,obj.nT], 0, 1);

      % sets-up the limit orientations rage
      obj.theta = previous.theta;
      obj.theta_range = previous.theta_range;

      obj = obj.addNoCollisionConstraint();
      obj = obj.addHandConstraint();
      obj = obj.addCostFunction();
    end

    function obj = addNoCollisionConstraint(obj)
      % constrains the pushers to lie within the complement of the object, 
      % represented through a set of convex polygonal regions

      % big-M
      M = 100;

      for t = 1:obj.nT
        theta = obj.theta;

        % regions assignment constraint
        for r = 1:obj.nr
          A = obj.shape.regions{r}.A;
          b = obj.shape.regions{r}.b;

          % constrains for each point in the first hand
          for j = 1:obj.n_e/2

            A = sparse(1,obj.nv);
            bi = -obj.region_old(r,j,t);
            Ai(1,obj.vars.region.i(r,(j-1)*obj.mult+1,t)) = -1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            if j < obj.n_e/2
              A = sparse(1,obj.nv);
              bi = -obj.region_old(r,j+1,t);
              Ai(1,obj.vars.region.i(r,j*obj.mult,t)) = -1;
              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end
            

            for k = 1:obj.mult
              idx = (j-1)*obj.mult+k;
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

                Ai(:, obj.vars.p_int.i(:,idx,t)) = A*(R_theta');
                Ai(:, obj.vars.region.i(r,idx,t)) = M;
                bi(:) = b + M;
                obj = obj.addLinearConstraints(Ai, bi, [], []);

                Ai = sparse(size(A, 1), obj.nv);
                bi = zeros(size(A, 1), 1);

                Ai(:, obj.vars.p_int.i(:,idx,t+1)) = A*(R_theta_1');
                Ai(:, obj.vars.region.i(r,idx,t)) = M;
                bi(:) = b + M;
                obj = obj.addLinearConstraints(Ai, bi, [], []);

                % next intermediate point must also satistfy this
                if j < obj.n_e/2
                  Ai = sparse(size(A, 1), obj.nv);
                  bi = zeros(size(A, 1), 1);

                  Ai(:, obj.vars.p_int.i(:,idx+1,t)) = A*(R_theta');
                  Ai(:, obj.vars.region.i(r,idx,t)) = M;
                  bi(:) = b + M;
                  obj = obj.addLinearConstraints(Ai, bi, [], []);
                end
              end
            end
          end

          % constraints on the Ssecond hand
          for j = 1:obj.n_e/2+1:obj.n_e

            A = sparse(1,obj.nv);
            bi = -obj.region_old(r,j,t);
            Ai(1,obj.vars.region.i(r,(j-1)*obj.mult+1,t)) = -1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            if j < obj.n_e
              A = sparse(1,obj.nv);
              bi = -obj.region_old(r,j+1,t);
              Ai(1,obj.vars.region.i(r,j*obj.mult,t)) = -1;
              obj = obj.addLinearConstraints(Ai, bi, [], []);
            end

            for k = 1:obj.mult
              idx = (j-1)*obj.mult+k;
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

                Ai(:, obj.vars.p_int.i(:,idx,t)) = A*(R_theta');
                Ai(:, obj.vars.region.i(r,idx,t)) = M;
                bi(:) = b + M;
                obj = obj.addLinearConstraints(Ai, bi, [], []);

                Ai = sparse(size(A, 1), obj.nv);
                bi = zeros(size(A, 1), 1);

                Ai(:, obj.vars.p_int.i(:,idx,t+1)) = A*(R_theta_1');
                Ai(:, obj.vars.region.i(r,idx,t)) = M;
                bi(:) = b + M;
                obj = obj.addLinearConstraints(Ai, bi, [], []);

                % next intermediate point must also satistfy this
                if j < obj.n_e
                  Ai = sparse(size(A, 1), obj.nv);
                  bi = zeros(size(A, 1), 1);

                  Ai(:, obj.vars.p_int.i(:,idx+1,t)) = A*(R_theta');
                  Ai(:, obj.vars.region.i(r,idx,t)) = M;
                  bi(:) = b + M;
                  obj = obj.addLinearConstraints(Ai, bi, [], []);
                end
              end
            end
          end
        end
        
        % assigns each pusher to one regions
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

    function obj = addCostFunction(obj)
      % Minimizes the separation between pushers in the final configuration
      for t = 1:obj.nT
        for j = 1:obj.n_e/2-1
          Qi = sparse(obj.nv,obj.nv);
          Qi(obj.vars.p_int.i(:,j,t),obj.vars.p_int.i(:,j,t)) = eye(2);
          Qi(obj.vars.p_int.i(1,j,t),obj.vars.p_int.i(1,j+1,t)) = -1;
          Qi(obj.vars.p_int.i(2,j,t),obj.vars.p_int.i(2,j+1,t)) = -1;
          Qi(obj.vars.p_int.i(:,j+1,t),obj.vars.p_int.i(:,j+1,t)) = eye(2);
          obj = obj.addCost(Qi,[],[]);
        end
        for j = obj.n_e/2+1:obj.n_e-1
          Qi = sparse(obj.nv,obj.nv);
          Qi(obj.vars.p_int.i(:,j,t),obj.vars.p_int.i(:,j,t)) = eye(2);
          Qi(obj.vars.p_int.i(1,j,t),obj.vars.p_int.i(1,j+1,t)) = -1;
          Qi(obj.vars.p_int.i(2,j,t),obj.vars.p_int.i(2,j+1,t)) = -1;
          Qi(obj.vars.p_int.i(:,j+1,t),obj.vars.p_int.i(:,j+1,t)) = eye(2);
          obj = obj.addCost(Qi,[],[]);
        end
      end
    end

    function obj = addHandConstraints(obj)
      for t = 2:obj.nT
        for i = 1:obj.n_e/2
          for k = 1:obj.mult
            Aeq = sparse(2,obj.nv);
            beq = obj.p_hand(:,1,t)-obj.p_hand(:,1,1);

            idx = (i-1)*mult+k;
            Aeq(:,obj.vars.p_int.i(:,idx,t)) = eye(2);
            Aeq(:,obj.vars.p_int.i(:,idx,1)) = -eye(2);

            obj = obj.addLinearConstraints([],[],Aeq,beq);
          end
        end

        for i = obj.n_e/2+1:obj.n_e
          for k = 1:obj.mult
            Aeq = sparse(2,obj.nv);
            beq = obj.p_hand(:,2,t)-obj.p_hand(:,2,1);

            idx = (i-1)*mult+k;
            Aeq(:,obj.vars.p_int.i(:,idx,t)) = eye(2);
            Aeq(:,obj.vars.p_int.i(:,idx,1)) = -eye(2);

            obj = obj.addLinearConstraints([],[],Aeq,beq);
          end
        end
      end
    end

    % end of methods
  end
end