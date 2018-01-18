classdef Arm < RigidBodyManipulator
	properties
		% Body indexes
		f1
		f2
		f3
		palm
	end

	methods
		function obj = Arm(options)
			% Defines a RigidBodyManipulator for an arm with IRB140 hand

			% sets up the default options
			if nargin < 1; options = struct('floating', false); end
			
			% defines the manipulator
			obj = RigidBodyManipulator([getDrakePath,'/examples/IRB140/urdf/irb_140_robotiq_simple_ati.urdf'],options);
			obj.UpdateJointlimits();

			% compiles the manipulator
			obj = compile(obj);
		end

		function obj = UpdateJointlimits(obj)
			% Updates the kinematic joint limits for the hand 

			% gets the default limits
			[hand_jl,hand_ju] = obj.getJointLimits();

			% updates the limits
			hand_ju(hand.getBody(hand_tip1).position_num) = pi/12;
			hand_ju(hand.getBody(hand_tip2).position_num) = pi/12;
			hand_ju(hand.getBody(hand_tip3).position_num) = pi/12;
			hand_ju(hand.getBody(obj.findJointId('finger_1_joint_1')).position_num) = pi/4;
			hand_ju(hand.getBody(obj.findJointId('finger_2_joint_1')).position_num) = pi/4;
			hand_ju(hand.getBody(obj.findJointId('finger_middle_joint_1')).position_num) = pi/4;
			hand_ju(hand.getBody(obj.findJointId('finger_1_joint_0')).position_num) = pi/6;
			hand_ju(hand.getBody(obj.findJointId('finger_2_joint_0')).position_num) = pi/6;
			hand_ju(hand.getBody(obj.findJointId('finger_middle_joint_0')).position_num) = pi/6;

			% updates the joint limits
			obj = obj.setJointLimits(hand_jl,hand_ju);
		end

		function obj = GetBodyIdx(obj)
			% returns the body indexes for the three fingers and the palm

			% geths the link Id for the fingers
			obj.f1 = obj.findLinkId('finger_1_link_3');
			obj.f2 = obj.findLinkId('finger_2_link_3');
			obj.f3 = obj.findLinkId('finger_middle_link_3');

			% gets the link Id for the palm
			obj.palm = obj.findLinkId('palm');
		end

		function q = SolveIK(obj,p)
			% dos the invers kinematics for the arm
			% so that hand is located at positions p

			% defines the constraints	
			ik_args = {};
			
			ik_args = [ik_args, {constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
								obj,obj.f1,[0;0;0],p(:,1),p(:,1))}];

			ik_args = [ik_args, {constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
								obj,obj.f2,[0;0;0],p(:,2),p(:,2))}];

			ik_args = [ik_args, {constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
								obj,obj.f3,[0;0;0],p(:,3),p(:,3))}];

			% gets the nominal config
			q0 = q0 = getZeroConfiguration(obj);

			% sets up the ik NLP
			ik_prob = InverseKinematics(obj, q0, ik_args{:});
			cost = Point(obj.getStateFrame,1);
			cost = 0*double(cost);
			ik_prob = ik_prob.setQ(diag(cost(1:15)));

			% solves the NLP for the configuration q
			q = ik_prob.solve(q0);
		end
	end
end