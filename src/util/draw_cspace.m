function frames = draw_cspace(p,q)
	for t = 1:size(q,2)
		clf
		subplot(1,2,1)
		v = [];
		for i = 1:length(p.polygons)
			th = q(3,t);
			rotmat = [cos(th),-sin(th);sin(th),cos(th)];
			v = q(1:2,t) + rotmat*p.polygons{i}.v;
			hullIndices = convhull(v(1,:),v(2,:));
			h = patch(v(1,hullIndices),v(2,hullIndices),[217,218,219]/255);
			set(h,'facealpha',.1);
			axis square
			hold on
		end

		xlim([-10,10]);
		ylim([-10,10]);
		set(gcf, 'color', 'white');
		grid off;
		set(gca,'Visible','off');
		% title('Workspace')

		% draws the C-space
		subplot(1,2,2)

		r = 0.3;
		
		[x,y,z] = sphere;
		r = 0.5;
		surf(x*r+q(1,t),y*r+q(2,t),z*r+q(3,t));

		xlim([-10,10]);
		ylim([-10,10]);
		zlim([-10,10]);
		set(gcf, 'color', 'white');
		grid off;
		% set(gca,'Visible','off');
		% title('C-Space')
		pause(0.1)

		frames(t) = getframe(gcf);
		drawnow
	end
end