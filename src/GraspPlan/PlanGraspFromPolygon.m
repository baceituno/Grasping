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
% @option Not yet used

  % defines the default options
  if nargin < 3; options = struct(); end
  if nargin < 2; n_contacts = 3; end

  % checks the unused option
  if ~isfield(options, 'quad_approx'); options.quad_approx = false; end
    if ~isfield(options, 'lin_sides'); options.lin_sides = 4; end

  assert(n_contacts > 2)

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

  % defines the optimization problem
  planner = MixedIntegerGraspPlanningProblem(safe_regions, n_contacts);

  % sets up the costs
  planner.q_cws = 1e-2;
  planner.q_u = 1e1;

  % add constraints
  planner = planner.addConvexRegions();
  planner = planner.addForceClosureConstraints();
  if options.quad_approx
    planner = planner.addQuadConvexDecompositionofBilinearTerms();
  else
    planner = planner.addLinConvexDecompositionofBilinearTerms(options.lin_sides);
  end
    
  planner = planner.addFrictionConesConstraints();
  % planner = planner.addKinematicConstraints();

  % solves the problem
  planner = planner.solve();
  
end