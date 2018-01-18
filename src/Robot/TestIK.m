%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualizes a grasp

checkDependency('lcmgl');
hand = RigidBodyManipulator([getDrakePath,'/examples/IRB140/urdf/irb_140_robotiq_simple_ati.urdf'],struct('floating',false));

% gets the initial configuration
hand_tip1 = hand.findLinkId('finger_1_link_3');
hand_tip2 = hand.findLinkId('finger_2_link_3');
hand_tip3 = hand.findLinkId('finger_middle_link_3');
palm = hand.findLinkId('palm');

% compites the robot
hand = compile(hand);

q0 = getZeroConfiguration(hand);
kinsol0 = hand.doKinematics(q0,0*q0,struct('use_mex',false));
palm_pos0 = hand.forwardKin(kinsol0,palm,[0;0;0]);

v = hand.constructVisualizer();

% transforms the grasp
p(:,1) = p(:,1) + palm_pos0;
p(:,2) = p(:,2) + palm_pos0;
p(:,3) = p(:,3) + palm_pos0;

% draws the box
verts = bsxfun(@times,palm_pos0+[0.23;0;0],ones(1,8))+repmat(box_size/2,1,8).*[1 1 1 1 -1 -1 -1 -1;1 1 -1 -1 1 1 -1 -1;1 -1 1 -1 1 -1 1 -1];
lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton,'box');
lcmgl.glColor3f(0,0,1);
lcmgl.polyhedron(verts(1,:),verts(2,:),verts(3,:));
lcmgl.switchBuffers();

for i = 1:3
	lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton,'p_n');
	lcmgl.glColor3f(1,0,0);
	lcmgl.box(p(:,i),[0.05,0.05,0.05]);
	lcmgl.switchBuffers();
end

% does the IK
ik_args = {};

ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
					hand,hand_tip1,[0;0;0],p(:,1),p(:,1))}];

ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
					hand,hand_tip2,[0;0;0],p(:,2),p(:,2))}];

ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
					hand,hand_tip3,[0;0;0],p(:,3),p(:,3))}];

ik_prob = InverseKinematics(hand, q0, ik_args{:});
cost = Point(hand.getStateFrame,1);
cost = 0*double(cost);
ik_prob = ik_prob.setQ(diag(cost(1:15)));
q = ik_prob.solve(q0);

v.draw(0,q)