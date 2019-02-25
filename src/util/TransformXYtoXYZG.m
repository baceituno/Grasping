function [x,y,g,qx,q_y,q_z,q_w] = TransformXYtoXYGQ(p1,p2)
	% Transforms finger positions of one hand in the plane (X,Y)^2 
	% to Yumi coordinates (X,Y,Z,G,Theta)

	% position of the hand
	x = (p1(1)+p2(1))/0.2;
	y = (p1(2)+p2(2))/0.2;

	% gripper separation
	g = 5*sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);

	% computes the angle of the hand
	theta = atan((p2(2)-p1(2))/(p2(1)-p1(1)));

	% defines the type of hand mode to use
	roll = pi/2;
	pitch = 0;
	yaw = theta;

	% from rpy to quat
	cy = cos(yaw*0.5);
    sy = sin(yaw*0.5);
    cr = cos(roll*0.5);
    sr = sin(roll*0.5);
    cp = cos(pitch*0.5);
    sp = sin(pitch*0.5);

    q_w = cy*cr*cp + sy*sr*sp
    q_x = cy*sr*cp - sy*cr*sp
    q_y = cy*cr*sp + sy*sr*cp
    q_z = sy*cr*cp - cy*sr*sp
end