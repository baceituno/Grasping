%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function planner = PlanGraspFromPolygon(safe_regions, n_contacts, options)
% planGraspFromPolygon: finds a set of contact locations and reaction forces
% that grasp the object with a relaxed angular equilibirium condition
% @params safe_regions a list of planar polytopic regions into which the
%                      contacts must be placed. Each
%                      contact must fall within at least one of the
%                      defined safe_regions. safe_regions should be a list
%                      of objects or structures which have fields 'A',
%                      'b', 'point', and 'normal' which define a region in
%                      v = [x,y,z,yaw] as
%                      A*v <= b AND normal'*v(1:3) == normal'*point
% @param options a struct of options:
% @option quad_approx: use a quadratic approximation of bilinear term.
%                      if set to false, the planner uses the linear approx.
% @option palm_pos: position of the palm of them robot
% @option use_kin:  constrain the fingers with the kinematics of the hand.

  % defines the default options
  if nargin < 3; options = struct(); end
  if nargin < 2; n_contacts = 3; end

  % checks the unused option
  if ~isfield(options, 'quad_approx'); options.quad_approx = false; end
  if ~isfield(options, 'palm_pos'); options.palm_pos = [-0.03,0,0]'; end
  if ~isfield(options, 'use_kin'); options.use_kin = false; end
  if ~isfield(options, 'logvars'); options.logvars = true; end

  assert(n_contacts > 2)

  % defines the optimization problem
  planner = MixedIntegerGraspPlanningProblem(safe_regions, n_contacts);

  % sets up the costs
  planner.q_cws = 1e-1;
  planner.q_u = 1;

  % add regional constraints
  if ~isfield(safe_regions(1), 'isV')
    for j = 1:length(safe_regions)
      sizecheck(safe_regions(j).A, [NaN, 3]);
      sizecheck(safe_regions(j).b, [size(safe_regions(j).A, 1), 1]);
      sizecheck(safe_regions(j).point, [3,1]);
      sizecheck(safe_regions(j).normal, [3,1]);
      safe_regions(j).normal = safe_regions(j).normal/norm(safe_regions(j).normal);
      for k = 1:size(safe_regions(j).A, 1)
        n = norm(safe_regions(j).A(k,:));
        safe_regions(j).A(k,:) = safe_regions(j).A(k,:)/n;
        safe_regions(j).b(k) = safe_regions(j).b(k)/n;
      end
    end
    planner = planner.addConvexRegions();
  else
    for j = 1:length(safe_regions)
      sizecheck(safe_regions(j).normal, [3,1]);
      safe_regions(j).normal = safe_regions(j).normal/norm(safe_regions(j).normal);
    end
    planner = planner.addVregions(options.logvars);
  end
  
  % add force closure constraints
  planner = planner.addForceClosureConstraints();
  if options.quad_approx
    planner = planner.addQuadConvexDecompositionofBilinearTerms();
  else
    planner = planner.addLinConvexDecompositionofBilinearTerms();
  end
    
  planner = planner.addFrictionConesConstraints();

  if options.use_kin && n_contacts == 3
    planner = planner.addKinematicConstraints(options.palm_pos);
  else
    planner = planner.constrainRegions();
  end

  % solves the problem
  t0 = tic();
  planner = planner.solveMosek();
  toc(t0)

end