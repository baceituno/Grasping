close all; clear all; clc;

p = PlanarShape('Poly2');
p_ = [-5,5;0,0];

%%%%%%%%%%%%
% Concepts %
%%%%%%%%%%%%

% C-SPACE
figure(11)

q = [linspace(0,-3,15),linspace(-3,0,35);...
	 linspace(-0,-3,25),linspace(-3,5,25);...
	  linspace(0,pi/4,25),linspace(pi/4,-pi/2,25)];

f1 = draw_cspace(p,q)

% C-Obstacle
figure(33)

q = [linspace(-5,-5,50);...
	 linspace(-2,2,30),linspace(2,0,20);...
	  linspace(0,pi/6,50)];

f2 = draw_cage(p,p_,q)

% Free-Space
figure(44)

q = [linspace(0,-3,15),linspace(-3,0,35);...
	 linspace(0,2,25),linspace(2,0,25);...
	  linspace(-pi/10,pi/10,50)];

f3 = draw_cage(p,p_,q)

% ESCAPING
figure(55)

q = [linspace(0,-3,15),linspace(-3,0,35);...
	 linspace(0,3,25),linspace(3,5,25);...
	  linspace(0,pi/4,50)];

f4 = draw_cage(p,p_,q)

%%%%%%%%%%%%%%%%%%
% Model elements %
%%%%%%%%%%%%%%%%%%

p_ = [-2,2;0,0];

figure(55)

q = [linspace(-0.8,0.8,20),linspace(0.8,-0.8,20),linspace(-0.8,0,10);...
	 linspace(-1,1,10),linspace(1,-1,10),linspace(-1,-1,10),linspace(-1,1,10),linspace(1,0,10);...
	  linspace(0,0,50)];

f5 = draw_cage(p,p_,q)

figure(110)

q = [linspace(-0.8,0.8,20),linspace(0.8,-0.8,20),linspace(-0.8,0,10);...
	 linspace(-1,1,10),linspace(1,-1,10),linspace(-1,-1,10),linspace(-1,1,10),linspace(1,0,10);...
	  linspace(-pi/10,pi/10,50)];

f6 = draw_cage(p,p_,q)

figure(220)

q = [linspace(-1,1,10),linspace(1,0,5),zeros(1,30),linspace(0,1,10),linspace(1,0,5);...
	 zeros(1,60);...
	 -pi/4*ones(1,15),linspace(-pi/4,pi/4,30),pi/4*ones(1,15)];

f7 = draw_cage(p,p_,q)
