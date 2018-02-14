classdef MixedIntegerGraspPlanningProblem < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (Mechatronics Group, USB C Laboratory)
%            and Hongkai Dai (Toyota Research Institute)
  properties
    n_contacts
    safe_regions

    tau_max = 1;
    mu_object = 1.0;
    num_edges = 4;

    com = zeros(3,1);

    q_cws = 1;
    q_u = 1;
  end

  methods
    function obj = MixedIntegerGraspPlanningProblem(safe_regions, n_contacts)
      % Constructs the optimization problem and declares the variables for each contact
      % @param n_contacts: number of fingers of the gripper
      assert(nargin > 0);
      if nargin < 2; n_contacts = 3; end

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();
      obj.n_contacts = n_contacts;
      obj.safe_regions = safe_regions;
      if isfield(safe_regions(1),'isV'); obj.com = obj.safe_regions(1).com; end

      % Contact locations
      obj = obj.addVariable('p', 'C', [3, obj.n_contacts], -inf, inf);
      obj = obj.addVariable('l', 'C', [3, obj.n_contacts], -inf, inf);
      % Contact forces
      obj = obj.addVariable('f_e', 'C', [3, obj.n_contacts],-inf, inf);

      % contact surface normal and force cones
      obj = obj.addVariable('alpha', 'C', [1, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('epsilon', 'C', [1, 1],0.1, inf);

      % friction cone edge multipliers
      obj = obj.addVariable('lambda_e', 'C', [obj.num_edges, obj.n_contacts],0, inf);

      % convex decomposition of bilinear terms
      obj = obj.addVariable('u_plus', 'C', [3, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('u_min', 'C', [3, obj.n_contacts],-inf, inf);
    end

    function obj = addConvexRegions(obj)
      % Add mixed-integer constraints that require that 
      % each contact lie within one of the safe regions described by A and b.
      % such that for contact i H(i,:) implies that A*fi < b
      % where H is a binary matrix with sum(H(i)) == 1
      nr = length(obj.safe_regions);
      obj = obj.addVariable('region', 'C', [nr, obj.n_contacts], 0, 1);

      % defines the log2 binary variables
      ny = ceil(log2(nr));
      obj = obj.addVariable('y', 'B', [ny, obj.n_contacts], 0, 1);

      % adds the sos1 constraint on lambda
      for i = 1:obj.n_contacts
        Aeq = sparse(1, obj.nv);
        beq = 1;
        Aeq(1, obj.vars.region.i(:,i)) = 1;
        obj = obj.addLinearConstraints([], [], Aeq, beq);        
      end

      % defines the gray coding on base 2
      codes = grayCodes(2,ny);
      codes = codes(1:nr,:);

      % adds the constraints on coding
      for i = 1:obj.n_contacts
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

      Ai = sparse((obj.n_contacts) * sum(cellfun(@(x) size(x, 1) + 2, {obj.safe_regions.A})), obj.nv);
      bi = zeros(size(Ai, 1), 1);
      offset_ineq = 0;
      Aeq = sparse(obj.n_contacts, obj.nv);
      beq = ones(obj.n_contacts, 1);
      offset_eq = 0;

      for r = 1:nr
        A = obj.safe_regions(r).A;
        b = obj.safe_regions(r).b;

        Ar = [A(:,1:2), sparse(size(A, 1), 1);
              obj.safe_regions(r).normal';
              -obj.safe_regions(r).normal'];
        br = [b;
              obj.safe_regions(r).normal' * obj.safe_regions(r).point;
              -obj.safe_regions(r).normal' * obj.safe_regions(r).point];
        s = size(Ar, 1);
        M = 10;
        for j = 1:obj.n_contacts
          Ai(offset_ineq + (1:s), obj.vars.p.i(:,j)) = Ar;
          Ai(offset_ineq + (1:s), obj.vars.region.i(r,j)) = M;
          bi(offset_ineq + (1:s)) = br + M;
          offset_ineq = offset_ineq + s;
        end
      end
      assert(offset_ineq == size(Ai, 1));
      for j = 1:obj.n_contacts
        Aeq(offset_eq + 1, obj.vars.region.i(:,j)) = 1;
        offset_eq = offset_eq + 1;
      end
      assert(offset_eq == size(Aeq, 1));
      obj = obj.addLinearConstraints(Ai, bi, Aeq, beq);
    end

    function obj = addVregions(obj, logvars)
      % Add mixed-integer constraints that require that 
      % each contact lie within one of the regions described by vertices
      % using a logarithmic number of integer variables, as decribed in
      % Modeling Disjunctive Constraints with a Logarithmic Number of Binary
      % Variables and Constraints by J. Vielma and G. Nemhauser, 2011.

      if logvars
        % defines the lambdas
        nr = length(obj.safe_regions);
        obj = obj.addVariable('region', 'C', [nr, obj.n_contacts], 0, 1);

        % defines the log2 binary variables
        ny = ceil(log2(nr));
        obj = obj.addVariable('y', 'B', [ny, obj.n_contacts], 0, 1);

        % adds the sos1 constraint on lambda
        for i = 1:obj.n_contacts
          Aeq = sparse(1, obj.nv);
          beq = 1;
          Aeq(1, obj.vars.region.i(:,i)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);        
        end

        % defines the gray coding on base 2
        codes = grayCodes(2,ny);
        codes = codes(1:nr,:);

        % adds the constraints on coding
        for i = 1:obj.n_contacts
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
      else
        nr = length(obj.safe_regions);
        obj = obj.addVariable('region', 'B', [nr, obj.n_contacts], 0, 1);
      end
        

      % constrains each contact to a single region
      use_3 = 0;
      if use_3 == 1
        disp('with big-M');
        % uses all weights
        % uses three weights only
        obj = obj.addVariable('w', 'C', [3, obj.n_contacts], 0, inf);
        
        % we define a big M
        M = 10;

        % for each contact and region
        for i = 1:obj.n_contacts
          idx = zeros(3,1);
          for j = 1:nr
            Ai = sparse(6, obj.nv);
            bi = M*ones(6,1);
            % the contact must be in a linear combination
            % of the vertices.
            Ai(1:3, obj.vars.p.i(:,i)) = eye(3);
            Ai(4:6, obj.vars.p.i(:,i)) = -eye(3);
            for k = 1:3
              Ai(1:3, obj.vars.w.i(k,i)) = -obj.safe_regions(j).vertices(:,k);
              Ai(4:6, obj.vars.w.i(k,i)) = obj.safe_regions(j).vertices(:,k);
            end
            % big-M constraint
            Ai(:, obj.vars.region.i(j,i)) = M;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end

        % all weights must add to one
        for i = 1:obj.n_contacts
          Aeq = sparse(1, obj.nv);
          beq = 1;
          Aeq(1, obj.vars.w.i(:,i)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);        
        end
      else
        disp('without big-M');
        % uses a weight for each vertex
        obj = obj.addVariable('w', 'C', [3*nr, obj.n_contacts], 0, inf);
        
        % makes all the weights add up to 1
        for i = 1:obj.n_contacts
            Aeq = sparse(1, obj.nv);
            beq = 1;
            Aeq(1, obj.vars.w.i(:,i)) = 1;
            obj = obj.addLinearConstraints([], [], Aeq, beq);        
        end

        % for each region
        for i = 1:obj.n_contacts
          % Çosntrains the contact
          Aeq = sparse(3, obj.nv);
          beq = zeros(3,1);
          
          Aeq(1:3, obj.vars.p.i(:,i)) = eye(3);
          for j = 1:nr
            idx = (j-1)*3 + 1;
            Aeq(1:3, obj.vars.w.i(idx:idx+2,i)) = -obj.safe_regions(j).vertices(:,1:3);
          end
          obj = obj.addLinearConstraints([], [], Aeq, beq);

          for j = 1:nr
            % adds requires the sum of the weights to be
            % at most lambda
            Ai = sparse(1, obj.nv);
            bi = 0;

            idx = (j-1)*3 + 1;
            Ai(1, obj.vars.w.i(idx:idx+2,i)) = 1;
            Ai(1, obj.vars.region.i(j,i)) = -1;
            
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end
      end
    end

    function obj = constrainRegions(obj)
      % constrains each contact to lie on a separate face
      nr = length(obj.safe_regions);
      for i = 1:nr
        Ai = sparse(1, obj.nv);
        bi = 1;
        Ai(1, obj.vars.region.i(i,:)) = 1;
        obj = obj.addLinearConstraints(Ai, bi, [], []);        
      end
    end

    function obj = addForceClosureConstraints(obj)
      % Constrains force closure over the object:
      % 1) Zero aggregated force
      % 2) Zero aggregated torque

      % constrains zero linear force
      Aeq = sparse(3, obj.nv);
      beq = zeros(3,1);
      for j = 1:obj.n_contacts
        Aeq(:, obj.vars.f_e.i(:,j)) = eye(3);
      end
      obj = obj.addLinearConstraints([],[],Aeq,beq);

      % gets the distance to CoM
      for j = 1:obj.n_contacts
        Aeq = sparse(3, obj.nv);
        beq = obj.com;
        Aeq(:, obj.vars.p.i(:,j)) = eye(3);
        Aeq(:, obj.vars.l.i(:,j)) = -eye(3);
        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % obtains the angular segment
      Aeq = sparse(3, obj.nv);
      beq = zeros(3, 1);
      for j = 1:obj.n_contacts
        Aeq(:, obj.vars.u_plus.i(:,j)) = eye(3);
        Aeq(:, obj.vars.u_min.i(:,j)) = -eye(3);        
      end
      obj = obj.addLinearConstraints([],[],Aeq,beq);
    end

    function obj = addKinematicConstraints(obj, palm_pos)
      % Constrains the fingers to lie on separate polytope volumes
      obj = obj.addVariable('r', 'C', [3, 1], -inf, inf);

      % models compliant actuation on z axis:
      Ai = sparse(1,obj.nv);
      bi = zeros(1,1);

      Ai(1,obj.vars.p.i(3,1)) = 1;
      Ai(1,obj.vars.p.i(3,3)) = -1;

      obj = obj.addLinearConstraints(Ai,bi,[],[]);

      % forces must point in a single direction
      Ai = sparse(3,obj.nv);
      bi = zeros(3,1);

      Ai(1,obj.vars.f_e.i(3,3)) = -1;
      bi(1,1) = -0.5;
      Ai(2,obj.vars.f_e.i(3,2)) = 1;
      bi(2,1) = 0.5;
      Ai(3,obj.vars.f_e.i(3,1)) = 1;
      bi(3,1) = 0.5;
      % obj = obj.addLinearConstraints(Ai,bi,[],[]);

      % all fingers in the same x
      Aeq = sparse(3,obj.nv);
      beq = zeros(3,1);

      Aeq(1,obj.vars.p.i(1,1)) = 1;
      Aeq(1,obj.vars.p.i(1,2)) = -1;

      Aeq(2,obj.vars.p.i(1,1)) = 1;
      Aeq(2,obj.vars.p.i(1,3)) = -1;

      obj = obj.addLinearConstraints([],[],Aeq,beq);

      % for finger 1
      Ai = sparse(6,obj.nv);
      bi = zeros(6,1);

      Ai(1:3,obj.vars.p.i(:,1)) = eye(3);
      Ai(1:3,obj.vars.r.i(:,1)) = -eye(3);
      bi(1:3) = [-palm_pos(1) + 10;-palm_pos(2) + 0.036; -palm_pos(3) + 10];
      Ai(4:6,obj.vars.p.i(:,1)) = -eye(3);
      Ai(4:6,obj.vars.r.i(:,1)) = eye(3);
      bi(4:6) = [palm_pos(1) + 10;palm_pos(2) - 0.036; -palm_pos(3) + 10];

      obj = obj.addLinearConstraints(Ai,bi,[],[]);

      % for finger 2
      Ai = sparse(6,obj.nv);
      bi = zeros(6,1);

      Ai(1:3,obj.vars.p.i(:,2)) = eye(3);
      Ai(1:3,obj.vars.r.i(:,1)) = -eye(3);
      bi(1:3) = [-palm_pos(1) + 10;-palm_pos(2) - 0.036; -palm_pos(3) + 10];
      Ai(4:6,obj.vars.p.i(:,2)) = -eye(3);
      Ai(4:6,obj.vars.r.i(:,1)) = eye(3);
      bi(4:6) = [palm_pos(1) + 10;palm_pos(2) + 0.036; -palm_pos(3) + 10];

      obj = obj.addLinearConstraints(Ai,bi,[],[]);

      % for finger 3
      Ai = sparse(6,obj.nv);
      bi = zeros(6,1);

      Ai(1:3,obj.vars.p.i(:,3)) = eye(3);
      Ai(1:3,obj.vars.r.i(:,1)) = -eye(3);
      bi(1:3) = [-palm_pos(1) + 10;-palm_pos(2); -palm_pos(3) + 10];
      Ai(4:6,obj.vars.p.i(:,3)) = -eye(3);
      Ai(4:6,obj.vars.r.i(:,1)) = eye(3);
      bi(4:6) = [palm_pos(1) + 10;palm_pos(2); -palm_pos(3) + 10];

      obj = obj.addLinearConstraints(Ai,bi,[],[]);
    end

    function obj = addLinConvexDecompositionofBilinearTerms(obj)
      % computes the angular momentum at each timestep
      % using a linear approximation of the cross
      % product for the angular momentum contibution
      % of each contact force:
      % l x f = (U+ - U-)/4 
      % with the following tight upper bound:
      % U+ > |l_bar+f_bar|, U- > |l_bar-f_bar|
      % where we add a quadratic cost at U+ and U-

      % introduces the quadratic approximation of bilinear terms
      % convex difference dissagregates for contact torques
      obj = obj.addVariable('a_p_d', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('a_m_d', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('b_p_e', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('b_m_e', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('c_p_f', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('c_m_f', 'C', [2, obj.n_contacts],-inf, inf);

      % defines the approximation constraints
      As = [1,1;-1,1;1,-1;-1,-1];
      bs = [1;1;1;1];
      sides = 4;

      % defines the descomposition elements
      for j = 1:obj.n_contacts
        % for each leg
        apd_idx = obj.vars.a_p_d.i(:,j);
        amd_idx = obj.vars.a_m_d.i(:,j);

        bpe_idx = obj.vars.b_p_e.i(:,j);
        bme_idx = obj.vars.b_m_e.i(:,j);

        cpf_idx = obj.vars.c_p_f.i(:,j);
        cmf_idx = obj.vars.c_m_f.i(:,j);

        l_1_idx = obj.vars.l.i(1,j);
        l_2_idx = obj.vars.l.i(2,j);
        l_3_idx = obj.vars.l.i(3,j);

        f_1_idx = obj.vars.f_e.i(1,j);
        f_2_idx = obj.vars.f_e.i(2,j);
        f_3_idx = obj.vars.f_e.i(3,j);        

        % defines the leg a and d
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, apd_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = -1;
        Aeq(1, f_2_idx) = 1;

        Aeq(2, l_2_idx) = 1;
        Aeq(2, f_3_idx) = 1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, amd_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = -1;
        Aeq(1, f_2_idx) = -1;

        Aeq(2, l_2_idx) = 1;
        Aeq(2, f_3_idx) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % defines the leg b and e
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, bpe_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = 1;
        Aeq(1, f_1_idx) = 1;

        Aeq(2, l_1_idx) = -1;
        Aeq(2, f_3_idx) = 1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, bme_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = 1;
        Aeq(1, f_1_idx) = -1;

        Aeq(2, l_1_idx) = -1;
        Aeq(2, f_3_idx) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % defines the leg c and f
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, cpf_idx) = -eye(2);
        
        Aeq(1, l_2_idx) = -1;
        Aeq(1, f_1_idx) = 1;

        Aeq(2, l_1_idx) = 1;
        Aeq(2, f_2_idx) = 1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, cmf_idx) = -eye(2);
        
        Aeq(1, l_2_idx) = -1;
        Aeq(1, f_1_idx) = -1;

        Aeq(2, l_1_idx) = 1;
        Aeq(2, f_2_idx) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % defines a linear approximation of the convex decomposition 
      for j = 1:obj.n_contacts
        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_p_d.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus.i(1,j)) = -bs;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_p_e.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_p_f.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_m_d.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min.i(1,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_m_e.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_m_f.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);
      end

      % minimizes the centroidal angular momentum of the motion
      for j = 1:obj.n_contacts
        Qi = sparse(obj.nv, obj.nv);
        Qi(obj.vars.u_min.i(:,j), obj.vars.u_min.i(:,j)) = eye(3);
        Qi(obj.vars.u_plus.i(:,j), obj.vars.u_plus.i(:,j)) = eye(3);
        obj = obj.addCost(obj.q_u*Qi, [], []);

        c = sparse(obj.nv, 1);
        c(obj.vars.u_min.i(:,j),1) = 1;
        obj = obj.addCost([], obj.q_u*c, []);

        c = sparse(obj.nv, 1);
        c(obj.vars.u_plus.i(:,j),1) = 1;
        obj = obj.addCost([], obj.q_u*c, []);
      end
    end

    function obj = addQuadConvexDecompositionofBilinearTerms(obj)
      % computes the angular momentum at each timestep
      % using a quadratic approximation of the cross
      % product for the angular momentum contibution
      % of each contact force:
      % l x f = (U+ - U-)/4 
      % with the following tight upper bound:
      % U+ >= (l_bar+f_bar)^2, U- >= (l_bar-f_bar)^2 
      % where we add a quadratic cost at U+ and U-

      % introduces the quadratic approximation of bilinear terms
      % convex difference dissagregates for contact torques
      obj = obj.addVariable('a_p_d', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('a_m_d', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('b_p_e', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('b_m_e', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('c_p_f', 'C', [2, obj.n_contacts],-inf, inf);
      obj = obj.addVariable('c_m_f', 'C', [2, obj.n_contacts],-inf, inf);

      % defines the descomposition elements
      for j = 1:obj.n_contacts
        % for each leg
        apd_idx = obj.vars.a_p_d.i(:,j);
        amd_idx = obj.vars.a_m_d.i(:,j);

        bpe_idx = obj.vars.b_p_e.i(:,j);
        bme_idx = obj.vars.b_m_e.i(:,j);

        cpf_idx = obj.vars.c_p_f.i(:,j);
        cmf_idx = obj.vars.c_m_f.i(:,j);

        l_1_idx = obj.vars.l.i(1,j);
        l_2_idx = obj.vars.l.i(2,j);
        l_3_idx = obj.vars.l.i(3,j);

        f_1_idx = obj.vars.f_e.i(1,j);
        f_2_idx = obj.vars.f_e.i(2,j);
        f_3_idx = obj.vars.f_e.i(3,j);        

        % defines the leg a and d
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, apd_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = -1;
        Aeq(1, f_2_idx) = 1;

        Aeq(2, l_2_idx) = 1;
        Aeq(2, f_3_idx) = 1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, amd_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = -1;
        Aeq(1, f_2_idx) = -1;

        Aeq(2, l_2_idx) = 1;
        Aeq(2, f_3_idx) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % defines the leg b and e
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, bpe_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = 1;
        Aeq(1, f_1_idx) = 1;

        Aeq(2, l_1_idx) = -1;
        Aeq(2, f_3_idx) = 1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, bme_idx) = -eye(2);
        
        Aeq(1, l_3_idx) = 1;
        Aeq(1, f_1_idx) = -1;

        Aeq(2, l_1_idx) = -1;
        Aeq(2, f_3_idx) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % defines the leg c and f
        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, cpf_idx) = -eye(2);
        
        Aeq(1, l_2_idx) = -1;
        Aeq(1, f_1_idx) = 1;

        Aeq(2, l_1_idx) = 1;
        Aeq(2, f_2_idx) = 1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2, obj.nv);
        beq = zeros(2, 1);

        Aeq(:, cmf_idx) = -eye(2);
        
        Aeq(1, l_2_idx) = -1;
        Aeq(1, f_1_idx) = -1;

        Aeq(2, l_1_idx) = 1;
        Aeq(2, f_2_idx) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % defines the convex part of the cross product
      for j = 1:obj.n_contacts
        % defines the variable U_plus
        quadcon = struct('Qc', sparse(obj.nv, obj.nv), 'q', zeros(obj.nv, 1), 'rhs', 0);

        quadcon.Qc(obj.vars.a_p_d.i(1,j), obj.vars.a_p_d.i(1,j)) = 1;
        quadcon.Qc(obj.vars.a_p_d.i(2,j), obj.vars.a_p_d.i(2,j)) = 1;
        quadcon.q(obj.vars.u_plus.i(1,j)) = -1;

        obj = obj.addQuadcon(quadcon);

        quadcon = struct('Qc', sparse(obj.nv, obj.nv), 'q', zeros(obj.nv, 1), 'rhs', 0);

        quadcon.Qc(obj.vars.b_p_e.i(1,j), obj.vars.b_p_e.i(1,j)) = 1;
        quadcon.Qc(obj.vars.b_p_e.i(2,j), obj.vars.b_p_e.i(2,j)) = 1;
        quadcon.q(obj.vars.u_plus.i(2,j)) = -1;

        obj = obj.addQuadcon(quadcon);

        quadcon = struct('Qc', sparse(obj.nv, obj.nv), 'q', zeros(obj.nv, 1), 'rhs', 0);

        quadcon.Qc(obj.vars.c_p_f.i(1,j), obj.vars.c_p_f.i(1,j)) = 1;
        quadcon.Qc(obj.vars.c_p_f.i(2,j), obj.vars.c_p_f.i(2,j)) = 1;
        quadcon.q(obj.vars.u_plus.i(3,j)) = -1;

        obj = obj.addQuadcon(quadcon);
      end
      
      % defines the concave part of the cross product
      for j = 1:obj.n_contacts 
        % defines the variable U_minus
        quadcon = struct('Qc', sparse(obj.nv, obj.nv), 'q', zeros(obj.nv, 1), 'rhs', 0);

        quadcon.Qc(obj.vars.a_m_d.i(1,j), obj.vars.a_m_d.i(1,j)) = 1;
        quadcon.Qc(obj.vars.a_m_d.i(2,j), obj.vars.a_m_d.i(2,j)) = 1;
        quadcon.q(obj.vars.u_min.i(1,j)) = -1;

        obj = obj.addQuadcon(quadcon);

        quadcon = struct('Qc', sparse(obj.nv, obj.nv), 'q', zeros(obj.nv, 1), 'rhs', 0);

        quadcon.Qc(obj.vars.b_m_e.i(1,j), obj.vars.b_m_e.i(1,j)) = 1;
        quadcon.Qc(obj.vars.b_m_e.i(2,j), obj.vars.b_m_e.i(2,j)) = 1;
        quadcon.q(obj.vars.u_min.i(2,j)) = -1;

        obj = obj.addQuadcon(quadcon);

        quadcon = struct('Qc', sparse(obj.nv, obj.nv), 'q', zeros(obj.nv, 1), 'rhs', 0);

        quadcon.Qc(obj.vars.c_m_f.i(1,j), obj.vars.c_m_f.i(1,j)) = 1;
        quadcon.Qc(obj.vars.c_m_f.i(2,j), obj.vars.c_m_f.i(2,j)) = 1;
        quadcon.q(obj.vars.u_min.i(3,j)) = -1;

        obj = obj.addQuadcon(quadcon);
      end

      % minimizes the centroidal angular momentum of the motion
      for j = 1:obj.n_contacts
        Qi = sparse(obj.nv, obj.nv);
        Qi(obj.vars.u_min.i(:,j), obj.vars.u_min.i(:,j)) = eye(3);
        Qi(obj.vars.u_plus.i(:,j), obj.vars.u_plus.i(:,j)) = eye(3);
        obj = obj.addCost(obj.q_u*Qi, [], []);
      end
    end

    function obj = addFrictionConesConstraints(obj)
      % Constrains the contact force to lie within
      % a friction cone
      if isempty(obj.safe_regions)
        obj.safe_regions = obj.seed_plan.obj.safe_regions;
      end
      nr = length(obj.safe_regions);
    
      % adds a reward to the cone robustness
      c = sparse(obj.nv, 1);
      c(obj.vars.epsilon.i(:,:),1) = 1;
      obj = obj.addCost([], -obj.q_cws*c, []);

      % computes the cone at flat ground
      theta = linspace(0,2*pi,obj.num_edges+1);
      theta = theta(1:end-1);
      edges_0 = [obj.mu_object*cos(theta);obj.mu_object*sin(theta);ones(1,obj.num_edges)]; 

      region_edges = cell(nr,1);

      % big M
      M = 10*obj.tau_max;
      l = 1;

      % computes all the cones of the regions
      for j = 1:nr
        R_fc = rotateVectorToAlign([0;0;1],obj.safe_regions(j).normal);
        region_edges{j} = R_fc*edges_0;
      end

      % for each region
      for r = 1:nr
        % for each contact
        for i = 1:obj.n_contacts
          % for current timestep
          f_idx = obj.vars.f_e.i(:,i);

          % constrains that the force stays in the cone of the region 
          % if the contact is assigned to it
          Ai = sparse(6,obj.nv);
          bi = M*ones(6,1);

          % H_{r,i} => f_{i,j} in FC_e
          Ai(1:3,f_idx) = eye(3);
          Ai(1:3,obj.vars.alpha.i(l,i)) = -obj.safe_regions(r).normal;
          Ai(1:3,obj.vars.region.i(r,i)) = M;

          Ai(4:6,f_idx) = -eye(3);
          Ai(4:6,obj.vars.alpha.i(l,i)) = obj.safe_regions(r).normal;
          Ai(4:6,obj.vars.region.i(r,i)) = M;

          % for each edge adds a positive weight
          for e = 1:obj.num_edges
            Ai(1:3,obj.vars.lambda_e.i(e,i)) = -region_edges{r}(1:3,e);
            Ai(4:6,obj.vars.lambda_e.i(e,i)) = region_edges{r}(1:3,e);
          end

          obj = obj.addLinearConstraints(Ai, bi, [], []);
        end
      end

      % uses the minimum alpha as e-margin
      for j = 1:obj.n_contacts
        Ai = sparse(1, obj.nv);
        bi = zeros(1,1);
        
        Ai(:, obj.vars.alpha.i(1,j)) = -1;
        Ai(:, obj.vars.epsilon.i(1,1)) = 1;
        
        obj = obj.addLinearConstraints(Ai,bi,[],[]);
      end

      % adds bounds to the normal force components
      for j = 1:obj.n_contacts
        for r = 1:nr 
          Ai = sparse(1, obj.nv);
          bi = obj.tau_max + M;
          
          Ai(1, obj.vars.f_e.i(:,j)) = obj.safe_regions(r).normal(:);
          Ai(1, obj.vars.region.i(r,j)) = M;
          
          obj = obj.addLinearConstraints(Ai,bi,[],[]);
        end
      end
    end
    % end of methods
  end
end