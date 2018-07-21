function mat = PlanarRotMat(ang)
% return the 2D rotation matrix for the angle given 
mat = [cos(ang) -sin(ang);...
       sin(ang) cos(ang)];
end