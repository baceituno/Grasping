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
			obj = obj@RigidBodyManipulator([getDrakePath,'/examples/IRB140/urdf/irb_140_robotiq_simple_ati.urdf'],options);
			
			% gets the link ID for the finger tips
			obj.f1 = obj.findLinkId('finger_1_link_3');
			obj.f2 = obj.findLinkId('finger_2_link_3');
			obj.f3 = obj.findLinkId('finger_middle_link_3');

			% gets the link Id for the palm
			obj.palm = obj.findLinkId('palm');

			obj = obj.UpdateJointlimits();

			% compiles the manipulator
			obj = compile(obj);
		end

		function ref = getPalmPos(obj)
			q0 = getZeroConfiguration(obj);
			kinsol0 = obj.doKinematics(q0,0*q0,struct('use_mex',false));
			palm_pos0 = obj.forwardKin(kinsol0,obj.palm,[0;0;0]);
			ref = palm_pos0+[0.23;0;-0.1]''
		end

		function obj = UpdateJointlimits(obj)
			% Updates the kinematic joint limits for the hand 

			% gets the default limits
			[hand_jl,hand_ju] = obj.getJointLimits();

			% updates the limits
			hand_ju(obj.getBody(obj.f1).position_num) = pi/12;
			hand_ju(obj.getBody(obj.f2).position_num) = pi/12;
			hand_ju(obj.getBody(obj.f3).position_num) = pi/12;
			hand_ju(obj.getBody(obj.findJointId('finger_1_joint_1')).position_num) = pi/4;
			hand_ju(obj.getBody(obj.findJointId('finger_2_joint_1')).position_num) = pi/4;
			hand_ju(obj.getBody(obj.findJointId('finger_middle_joint_1')).position_num) = pi/4;
			hand_ju(obj.getBody(obj.findJointId('finger_1_joint_0')).position_num) = pi/6;
			hand_ju(obj.getBody(obj.findJointId('finger_2_joint_0')).position_num) = pi/6;
			hand_ju(obj.getBody(obj.findJointId('finger_middle_joint_0')).position_num) = pi/6;

			% updates the joint limits
			obj = obj.setJointLimits(hand_jl,hand_ju);
		end

		function q = SolveIK(obj,p)
			% dos the invers kinematics for the arm
			% so that hand is located at positions p

			tol = 0.01*ones(3,1);

			% defines the constraints	
			ik_args = {};
			
			ik_args = [ik_args, {constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
								obj,obj.f1,[0;0;0],p(:,1) - tol,p(:,1) + tol)}];

			ik_args = [ik_args, {constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
								obj,obj.f2,[0;0;0],p(:,2) - tol,p(:,2) + tol)}];

			ik_args = [ik_args, {constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,...
								obj,obj.f3,[0;0;0],p(:,3) - tol,p(:,3) + tol)}];

			% gets the nominal config
			q0 = getZeroConfiguration(obj);

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