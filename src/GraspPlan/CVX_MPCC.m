classdef Quad_MixedIntegerMotionPlanningProblem < Quad_MixedIntegerConvexProgram
% Structure for Quadruped motion planning via Mixed Integer Convex Optimization
% Currently implemented with numerical variables for kinematic and dynamic constraints
% using tight upper bound relaxation for angular dynamics and mixed-{0,1} variables
% for safe region assignment and body-yaw orientation

  properties
    % robot model
    Quadruped;
    
    % user varaibles
    n_legs = 4;
    nsteps;
    ts;
    seed_plan;
    L_leg = 0.487;
    offset = [0.717, -0.717, pi - 0.717, -pi + 0.717];

    % regions parameters
    dz = 0.15;
    
    % trajectory parameters
    dt = 0.05;
    
    % conditioning parameters
    normalization = 1;
    l_norm = 1;
    f_norm = 1;

    % robot parameters
    r_v = [0.1;0.1;0.01];
    d_lim = 0.1;
    l_bnd = 0.3;
    mass = 84.896;
    I = [4.1927 0.0279 -0.2036; 0.0279 11.4264 -0.0045; -0.2036 -0.0045 12.501];
    g = 9.8;
    com_height = 0.5268;
    l_max = 0.45;
    swing_ht = 0.15;
    
    % optimization parameters
    q_dko = 1e1;
    q_cws = 1e-1;
    q_a = 1e-2;
    q_f = 1e-6;

    t_max = [90, 60, 90];

    % stability parameters
    num_edges = 4;
    mu_ground = 0.7;

    weights;
    max_distance = 500;
    pose_indices = [1,2,3,6];
  end

  methods
    function obj = CVX_MPCC()
      % CoM position and derivatives
      obj = obj.addVariable('r', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('r_nominal', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('dr', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('ddr', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('dddr', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('ddddr', 'C', [3, obj.ts],-inf, inf);
      
      % angular momentum
      obj = obj.addVariable('dk', 'C', [3, obj.ts],-inf, inf);

      % ground reaction forces
      obj = obj.addVariable('f_e_1', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('f_e_2', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('f_e_3', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('f_e_4', 'C', [3, obj.ts],-inf, inf);

      % contact surface normal and force cones
      obj = obj.addVariable('alpha', 'C', [4, obj.ts],0, inf);

      obj = obj.addVariable('lambda_e_1', 'C', [obj.num_edges, obj.ts],0, inf);
      obj = obj.addVariable('lambda_e_2', 'C', [obj.num_edges, obj.ts],0, inf);
      obj = obj.addVariable('lambda_e_3', 'C', [obj.num_edges, obj.ts],0, inf);
      obj = obj.addVariable('lambda_e_4', 'C', [obj.num_edges, obj.ts],0, inf);

      % end effector positions
      obj = obj.addVariable('p_1', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('p_2', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('p_3', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('p_4', 'C', [3, obj.ts], -inf, inf);

      obj = obj.addVariable('dp_1', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('dp_2', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('dp_3', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('dp_4', 'C', [3, obj.ts], -inf, inf);

      obj = obj.addVariable('ddp_1', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('ddp_2', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('ddp_3', 'C', [3, obj.ts], -inf, inf);
      obj = obj.addVariable('ddp_4', 'C', [3, obj.ts], -inf, inf);

      obj = obj.addVariable('l_1', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('l_2', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('l_3', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('l_4', 'C', [3, obj.ts],-inf, inf);

      obj = obj.addVariable('u_plus_1', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('u_min_1', 'C', [3, obj.ts],-inf, inf);

      obj = obj.addVariable('u_plus_2', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('u_min_2', 'C', [3, obj.ts],-inf, inf);

      obj = obj.addVariable('u_plus_3', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('u_min_3', 'C', [3, obj.ts],-inf, inf);

      obj = obj.addVariable('u_plus_4', 'C', [3, obj.ts],-inf, inf);
      obj = obj.addVariable('u_min_4', 'C', [3, obj.ts],-inf, inf);
    end

    function obj = addQuadraticGoalObjective(obj, goal_pos)
      % For each index j in step_indices, add a caudratic cost of the form:
      % weights(j) * (footstep(j) - xgoal(6 - N +j))' * Qfoal * (footstep(j) - xgoal(6 - N +j)
      w_goal = diag([1e2 1e2 1e2]);
      % adds a cost to the final CoM position
      Qi = sparse([], [], [], obj.nv, obj.nv, 4);
      ci = sparse(obj.nv, 1);
      Qi(obj.vars.r.i(:,end), obj.vars.r.i(:,end)) = w_goal;
      ci(obj.vars.r.i(:,end)) = -2 * w_goal * goal_pos;
      objcon_i = goal_pos' * w_goal * goal_pos;
      obj = obj.addCost(Qi, ci, objcon_i);
    end

    function obj = addIntegrationConstraints(obj)
      % Performs an euler integration over the ODE 
      for j = 2:obj.ts
        Aeq = sparse(12, obj.nv);
        beq = zeros(12, 1);

        Aeq(1:3, obj.vars.r.i(:,j)) = eye(3);
        Aeq(1:3, obj.vars.r.i(:,j-1)) = -eye(3);
        Aeq(1:3, obj.vars.dr.i(:,j)) = -obj.dt*eye(3)/2;
        Aeq(1:3, obj.vars.dr.i(:,j-1)) = -obj.dt*eye(3)/2;

        Aeq(4:6, obj.vars.dr.i(:,j)) = eye(3);
        Aeq(4:6, obj.vars.dr.i(:,j-1)) = -eye(3);
        Aeq(4:6, obj.vars.ddr.i(:,j)) = -obj.dt*eye(3);

        Aeq(7:9, obj.vars.th.i(:,j)) = obj.I;
        Aeq(7:9, obj.vars.th.i(:,j-1)) = -obj.I;
        Aeq(7:9, obj.vars.k.i(:,j)) = -obj.dt*eye(3);

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(6, obj.nv);
        beq = zeros(6, 1);

        Aeq(1:3, obj.vars.ddr.i(:,j)) = eye(3);
        Aeq(1:3, obj.vars.ddr.i(:,j-1)) = -eye(3);
        Aeq(1:3, obj.vars.dddr.i(:,j)) = -obj.dt*eye(3);

        Aeq(4:6, obj.vars.dddr.i(:,j)) = eye(3);
        Aeq(4:6, obj.vars.dddr.i(:,j-1)) = -eye(3);
        Aeq(4:6, obj.vars.ddddr.i(:,j)) = -obj.dt*eye(3);

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(12, obj.nv);
        beq = zeros(12, 1);

        Aeq(1:3, obj.vars.p_1.i(:,j)) = eye(3);
        Aeq(1:3, obj.vars.p_1.i(:,j-1)) = -eye(3);
        Aeq(1:3, obj.vars.dp_1.i(:,j)) = -obj.dt*eye(3);

        Aeq(4:6, obj.vars.p_2.i(:,j)) = eye(3);
        Aeq(4:6, obj.vars.p_2.i(:,j-1)) = -eye(3);
        Aeq(4:6, obj.vars.dp_2.i(:,j)) = -obj.dt*eye(3);

        Aeq(7:9, obj.vars.p_3.i(:,j)) = eye(3);
        Aeq(7:9, obj.vars.p_3.i(:,j-1)) = -eye(3);
        Aeq(7:9, obj.vars.dp_3.i(:,j)) = -obj.dt*eye(3);

        Aeq(10:12, obj.vars.p_4.i(:,j)) = eye(3);
        Aeq(10:12, obj.vars.p_4.i(:,j-1)) = -eye(3);
        Aeq(10:12, obj.vars.dp_4.i(:,j)) = -obj.dt*eye(3);

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(12, obj.nv);
        beq = zeros(12, 1);

        Aeq(1:3, obj.vars.dp_1.i(:,j)) = eye(3);
        Aeq(1:3, obj.vars.dp_1.i(:,j-1)) = -eye(3);
        Aeq(1:3, obj.vars.ddp_1.i(:,j)) = -obj.dt*eye(3);

        Aeq(4:6, obj.vars.dp_2.i(:,j)) = eye(3);
        Aeq(4:6, obj.vars.dp_2.i(:,j-1)) = -eye(3);
        Aeq(4:6, obj.vars.ddp_2.i(:,j)) = -obj.dt*eye(3);

        Aeq(7:9, obj.vars.dp_3.i(:,j)) = eye(3);
        Aeq(7:9, obj.vars.dp_3.i(:,j-1)) = -eye(3);
        Aeq(7:9, obj.vars.ddp_3.i(:,j)) = -obj.dt*eye(3);

        Aeq(10:12, obj.vars.dp_4.i(:,j)) = eye(3);
        Aeq(10:12, obj.vars.dp_4.i(:,j-1)) = -eye(3);
        Aeq(10:12, obj.vars.ddp_4.i(:,j)) = -obj.dt*eye(3);

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end
    end
    
    function obj = addDynamicConstraints(obj)
      % Adds the constraints that define the centroidal dynamics
      % of the robot's motion.
      % This includes the linear and angular momentum dynamics,
      % the contact forces contribution and the convex difference
      % of the angular momentum contribution of each contact force.

      % defines the centroidal dynamics of the motion
      for j = 1:obj.ts
        % defines the linear dynamics
        Aeq = sparse(3, obj.nv);
        beq = [0;0;-obj.mass*obj.g];

        Aeq(:, obj.vars.ddr.i(:,j)) = obj.mass*eye(3);
        Aeq(:, obj.vars.f_e_1.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.f_e_2.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.f_e_3.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.f_e_4.i(:,j)) = -eye(3);

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % defines the auxiliar variables
      for j = 1:obj.ts
        % defines the variable l
        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);

        Aeq(:, obj.vars.l_1.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.p_1.i(:,j)) = eye(3);
        Aeq(:, obj.vars.r.i(:,j)) = -eye(3);

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);
         
        Aeq(:, obj.vars.l_2.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.p_2.i(:,j)) = eye(3);
        Aeq(:, obj.vars.r.i(:,j)) = -eye(3);
         
        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);

        Aeq(:, obj.vars.l_3.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.p_3.i(:,j)) = eye(3);
        Aeq(:, obj.vars.r.i(:,j)) = -eye(3);
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);

        Aeq(:, obj.vars.l_4.i(:,j)) = -eye(3);
        Aeq(:, obj.vars.p_4.i(:,j)) = eye(3);
        Aeq(:, obj.vars.r.i(:,j)) = -eye(3);
         
        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % encodes input limits
      a_max = [2; 1; .6];
      v_max = [.5; .4; .5];

      for j = 2:obj.ts-1
        Ai = sparse(6, obj.nv);
        bi = [a_max; a_max];

        Ai(1:3, obj.vars.ddr.i(:,j)) = eye(3);
        Ai(4:6, obj.vars.ddr.i(:,j)) = -eye(3);
        obj = obj.addLinearConstraints(Ai,bi,[],[]);

        Ai = sparse(6, obj.nv);
        bi = [v_max; v_max];

        Ai(1:3, obj.vars.dr.i(:,j)) = eye(3);
        Ai(4:6, obj.vars.dr.i(:,j)) = -eye(3);
        obj = obj.addLinearConstraints(Ai,bi,[],[]);
      end
    end

    function obj = addRunningCost(obj)
      % Adds the running costs for the optimization problem

      for j = 1:obj.ts
        % minimizes the CoM acceleration
        Qi = sparse(obj.nv, obj.nv);
        Qi(obj.vars.ddr.i(:,j), obj.vars.ddr.i(:,j)) = diag([1 1 1]);
        obj = obj.addCost(obj.q_a*Qi, [], []);

        % minimizes the contact forces
        Qi = sparse(obj.nv, obj.nv);
        Qi(obj.vars.f_e_1.i(:,j), obj.vars.f_e_1.i(:,j)) = eye(3);
        Qi(obj.vars.f_e_2.i(:,j), obj.vars.f_e_2.i(:,j)) = eye(3);
        Qi(obj.vars.f_e_3.i(:,j), obj.vars.f_e_3.i(:,j)) = eye(3);
        Qi(obj.vars.f_e_4.i(:,j), obj.vars.f_e_4.i(:,j)) = eye(3);
        obj = obj.addCost(obj.q_f*Qi, [], []);
      end
    end

    function obj = addCVXDC(obj)
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
      obj = obj.addVariable('a_p_d_1', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_d_1', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_p_e_1', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_m_e_1', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_p_f_1', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_m_f_1', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('a_p_d_2', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_d_2', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_p_e_2', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_m_e_2', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_p_f_2', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_m_f_2', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('a_p_d_3', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_d_3', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_p_e_3', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_m_e_3', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_p_f_3', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_m_f_3', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('a_p_d_4', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_d_4', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_p_e_4', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('b_m_e_4', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_p_f_4', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('c_m_f_4', 'C', [2, obj.ts],-inf, inf);

      % defines the angular dynamics
      for j = 1:obj.ts
        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);
        
        Aeq(:, obj.vars.u_plus_1.i(:,j)) = obj.normalization*eye(3)/4;
        Aeq(:, obj.vars.u_min_1.i(:,j)) = -obj.normalization*eye(3)/4;
        
        Aeq(:, obj.vars.u_plus_2.i(:,j)) = obj.normalization*eye(3)/4;
        Aeq(:, obj.vars.u_min_2.i(:,j)) = -obj.normalization*eye(3)/4;
        
        Aeq(:, obj.vars.u_plus_3.i(:,j)) = obj.normalization*eye(3)/4;
        Aeq(:, obj.vars.u_min_3.i(:,j)) = -obj.normalization*eye(3)/4;
        
        Aeq(:, obj.vars.u_plus_4.i(:,j)) = obj.normalization*eye(3)/4;
        Aeq(:, obj.vars.u_min_4.i(:,j)) = -obj.normalization*eye(3)/4;
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % defines the descomposition elements
      for j = 1:obj.ts
        % for each leg
        for l = 1:4
          if l == 1
            apd_idx = obj.vars.a_p_d_1.i(:,j);
            amd_idx = obj.vars.a_m_d_1.i(:,j);

            bpe_idx = obj.vars.b_p_e_1.i(:,j);
            bme_idx = obj.vars.b_m_e_1.i(:,j);

            cpf_idx = obj.vars.c_p_f_1.i(:,j);
            cmf_idx = obj.vars.c_m_f_1.i(:,j);

            l_1_idx = obj.vars.l_1.i(1,j);
            l_2_idx = obj.vars.l_1.i(2,j);
            l_3_idx = obj.vars.l_1.i(3,j);

            f_1_idx = obj.vars.f_e_1.i(1,j);
            f_2_idx = obj.vars.f_e_1.i(2,j);
            f_3_idx = obj.vars.f_e_1.i(3,j);
          elseif l == 2
            apd_idx = obj.vars.a_p_d_2.i(:,j);
            amd_idx = obj.vars.a_m_d_2.i(:,j);

            bpe_idx = obj.vars.b_p_e_2.i(:,j);
            bme_idx = obj.vars.b_m_e_2.i(:,j);

            cpf_idx = obj.vars.c_p_f_2.i(:,j);
            cmf_idx = obj.vars.c_m_f_2.i(:,j);

            l_1_idx = obj.vars.l_2.i(1,j);
            l_2_idx = obj.vars.l_2.i(2,j);
            l_3_idx = obj.vars.l_2.i(3,j);

            f_1_idx = obj.vars.f_e_2.i(1,j);
            f_2_idx = obj.vars.f_e_2.i(2,j);
            f_3_idx = obj.vars.f_e_2.i(3,j);
          elseif l == 3
            apd_idx = obj.vars.a_p_d_3.i(:,j);
            amd_idx = obj.vars.a_m_d_3.i(:,j);

            bpe_idx = obj.vars.b_p_e_3.i(:,j);
            bme_idx = obj.vars.b_m_e_3.i(:,j);

            cpf_idx = obj.vars.c_p_f_3.i(:,j);
            cmf_idx = obj.vars.c_m_f_3.i(:,j);

            l_1_idx = obj.vars.l_3.i(1,j);
            l_2_idx = obj.vars.l_3.i(2,j);
            l_3_idx = obj.vars.l_3.i(3,j);

            f_1_idx = obj.vars.f_e_3.i(1,j);
            f_2_idx = obj.vars.f_e_3.i(2,j);
            f_3_idx = obj.vars.f_e_3.i(3,j);
          else
            apd_idx = obj.vars.a_p_d_4.i(:,j);
            amd_idx = obj.vars.a_m_d_4.i(:,j);

            bpe_idx = obj.vars.b_p_e_4.i(:,j);
            bme_idx = obj.vars.b_m_e_4.i(:,j);

            cpf_idx = obj.vars.c_p_f_4.i(:,j);
            cmf_idx = obj.vars.c_m_f_4.i(:,j);

            l_1_idx = obj.vars.l_4.i(1,j);
            l_2_idx = obj.vars.l_4.i(2,j);
            l_3_idx = obj.vars.l_4.i(3,j);

            f_1_idx = obj.vars.f_e_4.i(1,j);
            f_2_idx = obj.vars.f_e_4.i(2,j);
            f_3_idx = obj.vars.f_e_4.i(3,j);
          end            

          % defines the leg a and d
          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, apd_idx) = -eye(2);
          
          Aeq(1, l_3_idx) = -1/obj.l_norm;
          Aeq(1, f_2_idx) = 1/obj.f_norm;

          Aeq(2, l_2_idx) = 1/obj.l_norm;
          Aeq(2, f_3_idx) = 1/obj.f_norm;

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, amd_idx) = -eye(2);
          
          Aeq(1, l_3_idx) = -1/obj.l_norm;
          Aeq(1, f_2_idx) = -1/obj.f_norm;

          Aeq(2, l_2_idx) = 1/obj.l_norm;
          Aeq(2, f_3_idx) = -1/obj.f_norm;

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          % defines the leg b and e
          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, bpe_idx) = -eye(2);
          
          Aeq(1, l_3_idx) = 1/obj.l_norm;
          Aeq(1, f_1_idx) = 1/obj.f_norm;

          Aeq(2, l_1_idx) = -1/obj.l_norm;
          Aeq(2, f_3_idx) = 1/obj.f_norm;

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, bme_idx) = -eye(2);
          
          Aeq(1, l_3_idx) = 1/obj.l_norm;
          Aeq(1, f_1_idx) = -1/obj.f_norm;

          Aeq(2, l_1_idx) = -1/obj.l_norm;
          Aeq(2, f_3_idx) = -1/obj.f_norm;

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          % defines the leg c and f
          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, cpf_idx) = -eye(2);
          
          Aeq(1, l_2_idx) = -1/obj.l_norm;
          Aeq(1, f_1_idx) = 1/obj.f_norm;

          Aeq(2, l_1_idx) = 1/obj.l_norm;
          Aeq(2, f_2_idx) = 1/obj.f_norm;

          obj = obj.addLinearConstraints([],[],Aeq,beq);

          Aeq = sparse(2, obj.nv);
          beq = zeros(2, 1);

          Aeq(:, cmf_idx) = -eye(2);
          
          Aeq(1, l_2_idx) = -1/obj.l_norm;
          Aeq(1, f_1_idx) = -1/obj.f_norm;

          Aeq(2, l_1_idx) = 1/obj.l_norm;
          Aeq(2, f_2_idx) = -1/obj.f_norm;

          obj = obj.addLinearConstraints([],[],Aeq,beq);
        end
      end

      % defines the convex part of the cross product
      sides = 4;
      As = [1,1;-1,1;1,-1;-1,-1];
      bs = [1;1;1;1];

      for j = 1:obj.ts
        % defines the variable U_plus_1
        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_p_d_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_1.i(1,j)) = -bs;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_p_e_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_1.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_p_f_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_1.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_m_d_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_1.i(1,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_m_e_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_1.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_m_f_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_1.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);
        
        % defines the variable U_plus_2
        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_p_d_2.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_2.i(1,j)) = -bs;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_p_e_2.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_2.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_p_f_2.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_2.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_m_d_2.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_2.i(1,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_m_e_2.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_2.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_m_f_2.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_2.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);
        
        % defines the variable U_plus 3
        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_p_d_3.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_3.i(1,j)) = -bs;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_p_e_3.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_3.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_p_f_3.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_3.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_m_d_3.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_3.i(1,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_m_e_3.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_3.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_m_f_3.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_3.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);
        
        % defines the variable U_plus 4
        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_p_d_4.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_4.i(1,j)) = -bs;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_p_e_4.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_4.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_p_f_4.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_4.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_m_d_4.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_4.i(1,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_m_e_4.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_4.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_m_f_4.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_4.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);
      end

      for j = 1:obj.ts
        % defines the variable U_plus_1
        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_p_d_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_1.i(1,j)) = -bs;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_p_e_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_1.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_p_f_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_plus_1.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.a_m_d_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_1.i(1,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.b_m_e_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_1.i(2,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(sides, obj.nv);
        bi = zeros(sides, 1);

        Ai(1:sides, obj.vars.c_m_f_1.i(1:2,j)) = As;
        Ai(1:sides, obj.vars.u_min_1.i(3,j)) = -bs;
        
        obj = obj.addLinearConstraints(Ai, bi, [],[]);
      end

      % minimizes the centroidal angular momentum of the motion
      for j = 1:obj.ts
        Qi = sparse(obj.nv, obj.nv);
        Qi(obj.vars.u_min_1.i(:,j), obj.vars.u_min_1.i(:,j)) = eye(3);
        Qi(obj.vars.u_min_2.i(:,j), obj.vars.u_min_2.i(:,j)) = eye(3);
        Qi(obj.vars.u_min_3.i(:,j), obj.vars.u_min_3.i(:,j)) = eye(3);
        Qi(obj.vars.u_min_4.i(:,j), obj.vars.u_min_4.i(:,j)) = eye(3);
        Qi(obj.vars.u_plus_1.i(:,j), obj.vars.u_plus_1.i(:,j)) = eye(3);
        Qi(obj.vars.u_plus_2.i(:,j), obj.vars.u_plus_2.i(:,j)) = eye(3);
        Qi(obj.vars.u_plus_3.i(:,j), obj.vars.u_plus_3.i(:,j)) = eye(3);
        Qi(obj.vars.u_plus_4.i(:,j), obj.vars.u_plus_4.i(:,j)) = eye(3);

        obj = obj.addCost(obj.q_dko*Qi, [], []);
      end
    end

    function obj = addCVXCC(obj)
      % introduces the quadratic approximation of bilinear terms
      % convex difference dissagregates for contact torques
      obj = obj.addVariable('a_p_b_1', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_b_1', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('a_p_b_2', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_b_2', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('a_p_b_3', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_b_3', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('a_p_b_4', 'C', [2, obj.ts],-inf, inf);
      obj = obj.addVariable('a_m_b_4', 'C', [2, obj.ts],-inf, inf);

      obj = obj.addVariable('u_p_1', 'C', [2, obj.ts],0, inf);
      obj = obj.addVariable('u_p_2', 'C', [2, obj.ts],0, inf);
      obj = obj.addVariable('u_p_3', 'C', [2, obj.ts],0, inf);
      obj = obj.addVariable('u_p_4', 'C', [2, obj.ts],0, inf);

      obj = obj.addVariable('u_m_1', 'C', [2, obj.ts],0, inf);
      obj = obj.addVariable('u_m_2', 'C', [2, obj.ts],0, inf);
      obj = obj.addVariable('u_m_3', 'C', [2, obj.ts],0, inf);
      obj = obj.addVariable('u_m_4', 'C', [2, obj.ts],0, inf);

      obj = obj.addVariable('vp_1', 'C', [1, obj.ts],0, inf);
      obj = obj.addVariable('vp_2', 'C', [1, obj.ts],0, inf);
      obj = obj.addVariable('vp_3', 'C', [1, obj.ts],0, inf);
      obj = obj.addVariable('vp_4', 'C', [1, obj.ts],0, inf);

      % velocity module
      As = [1,1,1;1,1,-1;1,-1,1;1,-1,-1;-1,1,1;-1,1,-1;-1,-1,1;-1,-1,-1];;

      for j = 1:obj.ts
        Ai = sparse(8, obj.nv);
        bi = zeros(8, 1);
        
        Ai(:, obj.vars.dp_1.i(:,j)) = As;
        Ai(:, obj.vars.vp_1.i(1,j)) = -1;
        
        obj = obj.addLinearConstraints(Ai,bi,[],[]);

        Ai = sparse(8, obj.nv);
        bi = zeros(8, 1);
        
        Ai(:, obj.vars.dp_2.i(:,j)) = As;
        Ai(:, obj.vars.vp_2.i(1,j)) = -1;
        
        obj = obj.addLinearConstraints(Ai,bi,[],[]);

        Ai = sparse(8, obj.nv);
        bi = zeros(8, 1);
        
        Ai(:, obj.vars.dp_3.i(:,j)) = As;
        Ai(:, obj.vars.vp_3.i(1,j)) = -1;
        
        obj = obj.addLinearConstraints(Ai,bi,[],[]);

        Ai = sparse(8, obj.nv);
        bi = zeros(8, 1);
        
        Ai(:, obj.vars.dp_4.i(:,j)) = As;
        Ai(:, obj.vars.vp_4.i(1,j)) = -1;
        
        obj = obj.addLinearConstraints(Ai,bi,[],[]);
      end

      % defines the complementary dynamics
      for j = 1:obj.ts
        Aeq = sparse(3, obj.nv);
        beq = zeros(3, 1);
        
        Aeq(:, obj.vars.u_p_1.i(:,j)) = eye(2)/4;
        Aeq(:, obj.vars.u_m_1.i(:,j)) = -eye(2)/4;
        
        Aeq(:, obj.vars.u_p_2.i(:,j)) = eye(2)/4;
        Aeq(:, obj.vars.u_m_2.i(:,j)) = -eye(2)/4;
        
        Aeq(:, obj.vars.u_p_3.i(:,j)) = eye(2)/4;
        Aeq(:, obj.vars.u_m_3.i(:,j)) = -eye(2)/4;
        
        Aeq(:, obj.vars.u_p_4.i(:,j)) = eye(2)/4;
        Aeq(:, obj.vars.u_m_4.i(:,j)) = -eye(2)/4;
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % f + z
        Aeq = sparse(4, obj.nv);
        beq = zeros(4, 1);
        
        Aeq(1, obj.vars.a_p_b_1.i(1,j)) = -1;
        Aeq(1, obj.vars.f_e_1.i(3,j)) = 1;
        Aeq(1, obj.vars.p_1.i(3,j)) = 1;

        Aeq(2, obj.vars.a_p_b_2.i(1,j)) = -1;
        Aeq(2, obj.vars.f_e_2.i(3,j)) = 1;
        Aeq(2, obj.vars.p_2.i(3,j)) = 1;

        Aeq(3, obj.vars.a_p_b_3.i(1,j)) = -1;
        Aeq(3, obj.vars.f_e_3.i(3,j)) = 1;
        Aeq(3, obj.vars.p_3.i(3,j)) = 1;

        Aeq(4, obj.vars.a_p_b_4.i(1,j)) = -1;
        Aeq(4, obj.vars.f_e_4.i(3,j)) = 1;
        Aeq(4, obj.vars.p_4.i(3,j)) = 1;
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % f - z
        Aeq = sparse(4, obj.nv);
        beq = zeros(4, 1);
        
        Aeq(1, obj.vars.a_m_b_1.i(1,j)) = -1;
        Aeq(1, obj.vars.f_e_1.i(3,j)) = 1;
        Aeq(1, obj.vars.p_1.i(3,j)) = -1;

        Aeq(2, obj.vars.a_m_b_2.i(1,j)) = -1;
        Aeq(2, obj.vars.f_e_2.i(3,j)) = 1;
        Aeq(2, obj.vars.p_2.i(3,j)) = -1;

        Aeq(3, obj.vars.a_m_b_3.i(1,j)) = -1;
        Aeq(3, obj.vars.f_e_3.i(3,j)) = 1;
        Aeq(3, obj.vars.p_3.i(3,j)) = -1;

        Aeq(4, obj.vars.a_m_b_4.i(1,j)) = -1;
        Aeq(4, obj.vars.f_e_4.i(3,j)) = 1;
        Aeq(4, obj.vars.p_4.i(3,j)) = -1;
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % f + v_z
        Aeq = sparse(4, obj.nv);
        beq = zeros(4, 1);
        
        Aeq(1, obj.vars.a_p_b_1.i(2,j)) = -1;
        Aeq(1, obj.vars.f_e_1.i(3,j)) = 1;
        Aeq(1, obj.vars.vp_1.i(1,j)) = 1;

        Aeq(2, obj.vars.a_p_b_2.i(2,j)) = -1;
        Aeq(2, obj.vars.f_e_2.i(3,j)) = 1;
        Aeq(2, obj.vars.vp_2.i(1,j)) = 1;

        Aeq(3, obj.vars.a_p_b_3.i(2,j)) = -1;
        Aeq(3, obj.vars.f_e_3.i(3,j)) = 1;
        Aeq(3, obj.vars.vp_3.i(1,j)) = 1;

        Aeq(4, obj.vars.a_p_b_4.i(2,j)) = -1;
        Aeq(4, obj.vars.f_e_4.i(3,j)) = 1;
        Aeq(4, obj.vars.vp_4.i(1,j)) = 1;
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);

        % f - v_z
        Aeq = sparse(4, obj.nv);
        beq = zeros(4, 1);
        
        Aeq(1, obj.vars.a_m_b_1.i(2,j)) = -1;
        Aeq(1, obj.vars.f_e_1.i(3,j)) = 1;
        Aeq(1, obj.vars.vp_1.i(1,j)) = -1;

        Aeq(2, obj.vars.a_m_b_2.i(2,j)) = -1;
        Aeq(2, obj.vars.f_e_2.i(3,j)) = 1;
        Aeq(2, obj.vars.vp_2.i(1,j)) = -1;

        Aeq(3, obj.vars.a_m_b_3.i(2,j)) = -1;
        Aeq(3, obj.vars.f_e_3.i(3,j)) = 1;
        Aeq(3, obj.vars.vp_3.i(1,j)) = -1;

        Aeq(4, obj.vars.a_m_b_4.i(2,j)) = -1;
        Aeq(4, obj.vars.f_e_4.i(3,j)) = 1;
        Aeq(4, obj.vars.vp_4.i(1,j)) = -1;
        
        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      % defines the linear decomposition of f*z = 0
      for j = 1:obj.ts
        % defines the variable U_1
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_1.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_p_1.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_1.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_m_1.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        % defines the variable U_2
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_2.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_p_2.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_2.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_m_2.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        % defines the variable U_3
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_3.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_p_3.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_3.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_m_3.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        % defines the variable U_4
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_4.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_p_4.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_4.i(1,j)) = [1;-1];
        Ai(:, obj.vars.u_m_4.i(1,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);
      end

      % defines the linear decomposition of f*|dp| = 0
      for j = 1:obj.ts
        % defines the variable U_1
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_1.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_p_1.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_1.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_m_1.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        % defines the variable U_2
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_2.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_p_2.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_2.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_m_2.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        % defines the variable U_3
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_3.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_p_3.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_3.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_m_3.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        % defines the variable U_4
        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_p_b_4.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_p_4.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);

        Ai = sparse(2, obj.nv);
        bi = zeros(2, 1);

        Ai(:, obj.vars.a_m_b_4.i(2,j)) = [1;-1];
        Ai(:, obj.vars.u_m_4.i(2,j)) = -1;

        obj = obj.addLinearConstraints(Ai, bi, [],[]);
      end

      % minimizes the upper bounds to bilinear terms
      for j = 1:obj.ts
        Qi = sparse(obj.nv, obj.nv);
        Qi(obj.vars.u_m_1.i(:,j), obj.vars.u_m_1.i(:,j)) = eye(2);
        Qi(obj.vars.u_m_2.i(:,j), obj.vars.u_m_2.i(:,j)) = eye(2);
        Qi(obj.vars.u_m_3.i(:,j), obj.vars.u_m_3.i(:,j)) = eye(2);
        Qi(obj.vars.u_m_4.i(:,j), obj.vars.u_m_4.i(:,j)) = eye(2);
        Qi(obj.vars.u_p_1.i(:,j), obj.vars.u_p_1.i(:,j)) = eye(2);
        Qi(obj.vars.u_p_2.i(:,j), obj.vars.u_p_2.i(:,j)) = eye(2);
        Qi(obj.vars.u_p_3.i(:,j), obj.vars.u_p_3.i(:,j)) = eye(2);
        Qi(obj.vars.u_p_4.i(:,j), obj.vars.u_p_4.i(:,j)) = eye(2);

        Qi(obj.vars.v_1.i(:,j), obj.vars.v_1.i(:,j)) = eye(2);
        Qi(obj.vars.v_2.i(:,j), obj.vars.v_2.i(:,j)) = eye(1);
        Qi(obj.vars.v_3.i(:,j), obj.vars.v_3.i(:,j)) = eye(1);
        Qi(obj.vars.v_4.i(:,j), obj.vars.v_4.i(:,j)) = eye(1);

        obj = obj.addCost(obj.q_dko*Qi, [], []);
      end
    end

    function obj = addKinematicConstraints(obj)
      % Adds MI constraints that assign each end-effector movement to 
      % a transfer-phase, defined by the gait constraints      
      r_v = obj.r_v;

      % constrains the CoM position
      for j = 2:obj.ts
        Aeq = sparse(2, obj.nv);
        beq = zeros(2,1);

        Aeq(:, obj.vars.r_nominal.i(1:2,j)) = 4*eye(2);
        Aeq(:, obj.vars.p_1.i(1:2,j)) = -eye(2);
        Aeq(:, obj.vars.p_2.i(1:2,j)) = -eye(2);
        Aeq(:, obj.vars.p_3.i(1:2,j)) = -eye(2);
        Aeq(:, obj.vars.p_4.i(1:2,j)) = -eye(2);
        % 
        Ai = sparse(6, obj.nv);
        bi = zeros(6,1);
        
        Ai(1, obj.vars.r.i(1,j)) = 1;
        Ai(1, obj.vars.r_nominal.i(1,j)) = -1;
        bi(1,1) = r_v(1);

        Ai(2, obj.vars.r.i(1,j)) = -1;
        Ai(2, obj.vars.r_nominal.i(1,j)) = 1;
        bi(2,1) = r_v(1);

        Ai(3, obj.vars.r.i(2,j)) = 1;
        Ai(3, obj.vars.r_nominal.i(2,j)) = -1;
        bi(3,1) = r_v(2);

        Ai(4, obj.vars.r.i(2,j)) = -1;
        Ai(4, obj.vars.r_nominal.i(2,j)) = 1;
        bi(4,1) = r_v(2);

        Ai(5, obj.vars.r.i(3,j)) = 1;
        Ai(5, obj.vars.r_nominal.i(3,j)) = -1;
        bi(5,1) = r_v(3);

        Ai(6, obj.vars.r.i(3,j)) = -1;
        Ai(6, obj.vars.r_nominal.i(3,j)) = 1;
        bi(6,1) = r_v(3);

        obj = obj.addLinearConstraints(Ai,bi,Aeq,beq);
      end

      % constrains the end effector swing height      
      for j = 1:obj.Ns:obj.ts-obj.Ns+1
        for l = 1:obj.Ns-2*obj.Nt
          % constrains in z
          Ai = sparse(2, obj.nv);
          bi = [obj.swing_ht; 0];
          Ai(1, obj.vars.p_1.i(3,j+obj.Nt-1+l)) = 1;
          Ai(1, obj.vars.p_1.i(3,j)) = -(obj.Ns-obj.Nt-l-1)/(obj.Ns-obj.Nt-1);
          Ai(1, obj.vars.p_1.i(3,j+obj.Ns-1)) = -l/(obj.Ns-obj.Nt-1);
          obj = obj.addLinearConstraints(Ai,bi,[],[]);

          Ai = sparse(2, obj.nv);
          bi = [obj.swing_ht; 0];
          Ai(1, obj.vars.p_2.i(3,j+obj.Nt-1+l)) = 1;
          Ai(1, obj.vars.p_2.i(3,j)) = -(obj.Ns-obj.Nt-l-1)/(obj.Ns-obj.Nt-1);
          Ai(1, obj.vars.p_2.i(3,j+obj.Ns-1)) = -l/(obj.Ns-obj.Nt-1);
          obj = obj.addLinearConstraints(Ai,bi,[],[]);

          Ai = sparse(2, obj.nv);
          bi = [obj.swing_ht; 0];
          Ai(1, obj.vars.p_3.i(3,j+obj.Nt-1+l)) = 1;
          Ai(1, obj.vars.p_3.i(3,j)) = -(obj.Ns-obj.Nt-l-1)/(obj.Ns-obj.Nt-1);
          Ai(1, obj.vars.p_3.i(3,j+obj.Ns-1)) = -l/(obj.Ns-obj.Nt-1);
          obj = obj.addLinearConstraints(Ai,bi,[],[]);

          Ai = sparse(2, obj.nv);
          bi = [obj.swing_ht; 0];
          Ai(1, obj.vars.p_4.i(3,j+obj.Nt-1+l)) = 1;
          Ai(1, obj.vars.p_4.i(3,j)) = -(obj.Ns-obj.Nt-l-1)/(obj.Ns-obj.Nt-1);
          Ai(1, obj.vars.p_4.i(3,j+obj.Ns-1)) = -l/(obj.Ns-obj.Nt-1);
          obj = obj.addLinearConstraints(Ai,bi,[],[]);
        end
      end 
    end

    function obj = addFlatFrictionConesConstraints(obj)
      % Constrains the friction force to lie within
      % a contact cone in flat ground

      % adds a reward to the cone robustness
      c = sparse(obj.nv, 1);
      c(obj.vars.alpha.i(:,:),1) = 1;
      obj = obj.addCost([], -obj.q_cws*c, []);

      % computes the cone at flat ground
      display('using flat terrain friction cones');
      theta = linspace(0,2*pi,obj.num_edges+1);
      theta = theta(1:end-1);
      edges_0 = [obj.mu_ground*cos(theta);obj.mu_ground*sin(theta);ones(1,obj.num_edges)];

      % constrains the forces to a cone
      for j = 1:obj.ts
        Aeq = sparse(12,obj.nv);
        beq = zeros(12,1);

        Aeq(1:3,obj.vars.f_e_1.i(:,j)) = eye(3);
        for e = 1:obj.num_edges
          Aeq(1:3,obj.vars.lambda_e_1.i(e,j)) = -edges_0(:,e);
        end
        Aeq(3,obj.vars.alpha.i(1,j)) = -1;

        Aeq(4:6,obj.vars.f_e_2.i(:,j)) = eye(3);
        for e = 1:obj.num_edges
          Aeq(4:6,obj.vars.lambda_e_2.i(e,j)) = -edges_0(:,e);
        end
        Aeq(6,obj.vars.alpha.i(2,j)) = -1;

        Aeq(7:9,obj.vars.f_e_3.i(:,j)) = eye(3);
        for e = 1:obj.num_edges
          Aeq(7:9,obj.vars.lambda_e_3.i(e,j)) = -edges_0(:,e);
        end
        Aeq(9,obj.vars.alpha.i(3,j)) = -1;

        Aeq(10:12,obj.vars.f_e_4.i(:,j)) = eye(3);
        for e = 1:obj.num_edges
          Aeq(10:12,obj.vars.lambda_e_4.i(e,j)) = -edges_0(:,e);
        end
        Aeq(12,obj.vars.alpha.i(4,j)) = -1;

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end
    end
    % end of methods
  end
end