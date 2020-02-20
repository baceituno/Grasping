function frames = draw_cage(p,p_,q)
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
		r = 0.3;
		th = 0:pi/50:2*pi;
		for j = 1:size(p_,2)
			x = p_(1,j);
			y = p_(2,j);
			xunit = r * cos(th) + x;
			yunit = r * sin(th) + y;
			h = fill(xunit, yunit,'k');
			set(h,'facealpha',201/255);
			hold on
		end

		xlim([-10,10]);
		ylim([-10,10]);
		set(gcf, 'color', 'white');
		grid off;
		set(gca,'Visible','off');
		title('Workspace')

		% draws the C-space slice
		subplot(1,2,2)
		for j = 1:size(p_,2)
			for i = 1:length(p.polygons)
				th = q(3,t);
				rotmat = [cos(th),-sin(th);sin(th),cos(th)];
				v = -rotmat*p.polygons{i}.v + p_(:,j);
				h = fill(v(1,:),v(2,:),[209,203,193]/255);
				axis square
				set(h,'facealpha',201/255);
				hold on
			end
		end

		r = 0.3;
		th = 0:pi/50:2*pi;
		xunit = r * cos(th) + q(1,t);
		yunit = r * sin(th) + q(2,t);
		h = fill(xunit, yunit,'r');
		set(h,'facealpha',.9);
		xlim([-10,10]);
		ylim([-10,10]);
		set(gcf, 'color', 'white');
		grid off;
		set(gca,'Visible','off');
		title('Free-Space')
		pause(0.1)

		frames(t) = getframe(gcf);
		drawnow
	end
end