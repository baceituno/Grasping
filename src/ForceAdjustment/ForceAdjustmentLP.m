classdef ForceAdjustmentLP < Quad_MixedIntegerConvexProgram
% Developed by Bernardo Aceituno-C (Mechatronics Group, USB C Laboratory)
%            and Hongkai Dai (Toyota Research Institute)
  properties
    n_contacts
    safe_regions

    tau_max = 1;
    mu_object = 1;
    num_edges = 4;
  end

  methods
    function obj = ForceAdjustmentLP(G, normals)
      % Constructs the optimization problem and declares the variables for each contact
      % @param n_contacts: number of fingers of the gripper
      assert(nargin > 0);

      % sets up the optimization
      obj = obj@Quad_MixedIntegerConvexProgram();

      obj.n_contacts = size(normals,2);
      % Contact forces
      obj = obj.addVariable('f_e', 'C', [3, obj.n_contacts],-inf, inf);

      % contact surface normal and force cones
      obj = obj.addVariable('alpha', 'C', [1, obj.n_contacts],0, inf);
      obj = obj.addVariable('epsilon', 'C', [1, 1],0.01, inf);
      obj = obj.addVariable('lambda_e', 'C', [obj.num_edges, obj.n_contacts],0, inf);

      % adds the cone edges
      region_edges = cell(obj.n_contacts,1);
      theta = linspace(0,2*pi,obj.num_edges+1);
      theta = theta(1:end-1);
      edges_0 = [obj.mu_object*cos(theta);obj.mu_object*sin(theta);ones(1,obj.num_edges)]; 

      % adds the grasp matrix
      Aeq = sparse(6,obj.nv);
      beq = zeros(6,1);
      Aeq(:,obj.vars.f_e.i(:,:)) = G;
      obj = obj.addLinearConstraints([], [], Aeq, beq);

      % computes all the cones of the regions
      for i = 1:obj.n_contacts
        R_fc = rotateVectorToAlign([0;0;1],normals(:,i));
        region_edges{i} = R_fc*edges_0;
      end

      % for each contact
      for i = 1:obj.n_contacts
        % constrains that the force stays in the cone of the region 
        % if the contact is assigned to it
        Aeq = sparse(3,obj.nv);
        beq = zeros(3,1);

        % H_{r,i} => f_{i,j} in FC_e
        Aeq(1:3,obj.vars.f_e.i(:,i)) = eye(3);
        Aeq(1:3,obj.vars.alpha.i(1,i)) = -normals(:,i);

        % for each edge adds a positive weight
        for e = 1:obj.num_edges
          Aeq(1:3,obj.vars.lambda_e.i(e,i)) = -region_edges{i}(1:3,e);
        end

        obj = obj.addLinearConstraints([], [], Aeq, beq);
      end      

      for j = 1:obj.n_contacts
        Ai = sparse(2, obj.nv);
        bi = [obj.tau_max;0];
        
        Ai(1, obj.vars.f_e.i(:,j)) = normals(:,j);
        Ai(2, obj.vars.epsilon.i(1,1)) = 1;
        Ai(2, obj.vars.alpha.i(1,j)) = -1;
        obj = obj.addLinearConstraints(Ai,bi,[],[]);
      end

      c = sparse(obj.nv, 1);
      c(obj.vars.epsilon.i(:,:),1) = 1;
      obj = obj.addCost([], -c, []);
    end
    % end of methods
  end
end