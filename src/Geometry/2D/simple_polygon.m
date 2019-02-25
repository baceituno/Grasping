function [shape] = simple_polygon(numSides,irregularity,spikeyness)

    if nargin < 2; irregularity = 0.1; end;
    if nargin < 3; spikeyness = 0.1; end;

    % The following algortih Returns a list of vertices for a random polygon, in CCW order.
    % Adapted from the solution: https://stackoverflow.com/questions/8997099/algorithm-to-generate-random-2d-polygon

    ctrX = 0;
    ctrY = 0; 
    aveRadius = 1; 
    numVerts = numSides;

    irregularity = clip( irregularity, 0,1 ) * 2*pi/ numVerts;
    spikeyness = clip( spikeyness, 0,1 ) * aveRadius;

    % generate n angle steps
    angleSteps = [];
    lower = (2*pi / numVerts) - irregularity;
    upper = (2*pi / numVerts) + irregularity;
    sum = 0;
    for i =1:numVerts
        tmp = unifrnd(lower, upper);
        angleSteps(i) = tmp;
        sum = sum + tmp;
    end

    % normalize the steps so that point 0 and point n+1 are the same
    k = sum / (2*pi);
    for i =1:numVerts
        angleSteps(i) = angleSteps(i) / k;
    end

    % now generate the points
    points = [];
    angle = unifrnd(0, 2*pi);
    for i =1:numVerts
        r_i = clip( normrnd(aveRadius, spikeyness), 0, 2*aveRadius);
        x = ctrX + r_i* cos(angle);
        y = ctrY + r_i* sin(angle);
        points(i,:)= [(x),(y)];
        angle = angle + angleSteps(i);
    end

    % centers the polygon
    x = points(:,1)-mean(points(:,1));
    y = points(:,2)-mean(points(:,2));

    % x = [1,0,-1,0]';
    % y = [0,1,0,-1]';

    figure(1984)
    fill(x,y,'b')

    % describes the line segments
    nv = length(x);
    pre_lines = {};

    for i = 1:nv
        pre_lines{i} = struct('v1', [], 'v2', [], 'opp', [], 'non_cod', [], 'angle', [], 'isCV', false);
        pre_lines{i}.v1 = [x(i);y(i)];
        idx_2 = i+1;
        if i == nv
            idx_2 = 1;
        end
        pre_lines{i}.v2 = [x(idx_2);y(idx_2)];

        % computes the angle of the segment, given that the vertexes are counted counterclockwise
        dy = abs(pre_lines{i}.v2(2)-pre_lines{i}.v1(2));
        dx = abs(pre_lines{i}.v2(1)-pre_lines{i}.v1(1));

        if pre_lines{i}.v2(2) > pre_lines{i}.v1(2)
            if pre_lines{i}.v2(1) > pre_lines{i}.v1(1)
                pre_lines{i}.angle = -(pi/2-atan(dy/dx));
            elseif pre_lines{i}.v2(1) == pre_lines{i}.v1(1)
                pre_lines{i}.angle = 0;
            else
                pre_lines{i}.angle = (pi/2-atan(dy/dx));
            end
        elseif pre_lines{i}.v2(2) == pre_lines{i}.v1(2)
            if pre_lines{i}.v2(1) > pre_lines{i}.v1(1)
                pre_lines{i}.angle = -pi/2;
            else
                pre_lines{i}.angle = pi/2;
            end
        else
            if pre_lines{i}.v2(1) > pre_lines{i}.v1(1)
                pre_lines{i}.angle = -pi/2 - atan(dy/dx);
            else
                pre_lines{i}.angle = pi/2 + atan(dy/dx);
            end
        end
    end

    k = 1;
    lines = {};

    % loads vertices to check for convexity
    px = x;
    py = y;
    signPoly = sign(1);

    % accounts for concave vertices
    for i = 1:nv
        % assings each line segment
        lines{k} = pre_lines{i};
        k = k + 1;

        % indexes
        idx_0 = i;
        idx_1 = i + 1;
        idx_2 = i + 2;

        if i == nv-1
            idx_2 = 1;
        end

        if i == nv
            idx_1 = 1;
            idx_2 = 2;
        end

        v1 = [px(idx_1) - px(idx_0), py(idx_1) - py(idx_0)];
        v2 = [px(idx_2) - px(idx_1), py(idx_2) - py(idx_1)]; 
        curr_signPoly = sign(det([v1; v2]));

        % check that the signs match
        if ~isequal(curr_signPoly, signPoly)
            lines{k} = struct('v1', [], 'v2', [], 'opp', [], 'non_cod', [], 'angle', 0, 'isCV', true);
            lines{k}.v1 = [pre_lines{idx_0}.v2(1);pre_lines{idx_0}.v2(2)];
            lines{k}.v2 = [pre_lines{idx_0}.v2(1);pre_lines{idx_0}.v2(2)];
            k = k + 1;
        end
    end



    nl = length(lines);

    for i = 1:nl
        % facet normal
        if lines{i}.isCV
            if i == 1
                deg = lines{i+1}.angle + lines{end}.angle;
                normal = [cos(deg/2);sin(deg/2)];
                normal = normal/norm(normal);
            elseif i == nl
                deg = lines{i-1}.angle + lines{1}.angle;
                normal = [cos(deg/2);sin(deg/2)];
                normal = normal/norm(normal);                       
            else
                deg = lines{i-1}.angle + lines{i+1}.angle;
                normal = [cos(deg/2);sin(deg/2)];
                normal = normal/norm(normal);
            end
        else
            normal = [cos(lines{i}.angle);sin(lines{i}.angle)];
            normal = normal/norm(normal);
        end 

        lines{i}.normal = normal;

        % opposition of the faces
        opposite = [];

        if lines{i}.isCV
            if i == 1
                opposite = [i+2:nl-1];
            elseif i == nl
                opposite = [2:i-2];
            else
                opposite = [1:i-2,i+2:nl];
            end
        else
            for j = 1:nl
                if lines{j}.isCV && j ~= i
                    opposite = [opposite, j];
                end
                if abs(lines{i}.angle - lines{j}.angle) == pi
                    opposite = [opposite, j];
                end
            end
        end

        lines{i}.opp = opposite;

        % parallelity of the faces
        paral = [];

        if lines{i}.isCV
            if i == 1
                kdx1 = nl;
                kdx2 = i+1;
            elseif i == nl
                kdx1 = i-i;
                kdx2 = 1;
            else
                kdx1 = i-1;
                kdx2 = i+1;
            end
            
            for j = 1:nl
                if lines{j}.isCV
                    idx1 = j - 1;
                    idx2 = j + 1;

                    if idx1 < 1; idx1 = nl; end;
                    if idx2 > nl; idx2 = 1; end;

                    dif_angs1 = abs(lines{kdx1}.angle - lines{idx1}.angle);
                    dif_angs2 = abs(lines{kdx1}.angle - lines{idx2}.angle);

                    dif_angs3 = abs(lines{kdx2}.angle - lines{idx1}.angle);
                    dif_angs4 = abs(lines{kdx2}.angle - lines{idx2}.angle);

                    if (dif_angs1 == 0) || (dif_angs1 == pi) || (dif_angs2 == 0) || (dif_angs2 == pi) || (dif_angs3 == 0) || (dif_angs3 == pi) || (dif_angs4 == 0) || (dif_angs4 == pi)
                        paral = [paral, j];
                    end
                else
                    if kdx1 < 1; kdx1 = 1; end;
                    if kdx2 > nl; kdx2 = nl; end;
                    dif_angs1 = abs(lines{kdx1}.angle - lines{j}.angle);
                    dif_angs2 = abs(lines{kdx2}.angle - lines{j}.angle);
                    if (dif_angs1 == 0) || (dif_angs1 == pi) || (dif_angs2 == 0) || (dif_angs2 == pi)
                        paral = [paral, j];
                    end
                end
            end

        else
            for j = 1:nl
                if lines{j}.isCV
                    idx1 = j - 1;
                    idx2 = j + 1;

                    if idx1 < 1
                        idx1 = nl;
                    end

                    if idx2 > nl
                        idx2 = 1;
                    end

                    dif_angs1 = abs(lines{i}.angle - lines{idx1}.angle);
                    dif_angs2 = abs(lines{i}.angle - lines{idx2}.angle);

                    if (dif_angs1 == 0) || (dif_angs1 == pi) || (dif_angs2 == 0) || (dif_angs2 == pi)
                        paral = [paral, j];
                    end
                else
                    dif_angs = abs(lines{i}.angle - lines{j}.angle);
                    if (dif_angs == 0) || (dif_angs == pi)
                        paral = [paral, j];
                    end
                end
            end
        end

        lines{i}.parallel_fac = paral;

        % co-directionaly of the faces
        noncod = [];

        if lines{i}.isCV
            if i == 1
                noncod = [i+2:nl-1];
            elseif i == nl
                noncod = [2:i-2];
            else
                noncod = [1:i-2,i+2:nl];
            end
        else
            for j = 1:nl
                if lines{j}.isCV && j ~= i
                    if j == 1
                        if i ~= nl && i ~= 2
                            noncod = [noncod, j];
                        end
                    elseif j == nl
                        if i ~= nl-1 && i ~= 1
                            noncod = [noncod, j];
                        end
                    else
                        if i ~= j-1 && i ~= j+1
                            noncod = [noncod, j];
                        end
                    end
                end
                if ~lines{j}.isCV
                    dif_angs = abs(lines{i}.angle - lines{j}.angle);
                    if (dif_angs > 0) && (dif_angs ~= pi);
                        noncod = [noncod, j];
                    end
                end
            end
        end

        lines{i}.non_cod = noncod;

        % co-equality of the faces
        noneq = [];

        if lines{i}.isCV
            if i == 1
                noneq = [i+2:nl-1];
            elseif i == nl
                noneq = [2:i-2];
            else
                noneq = [1:i-2,i+2:nl];
            end
        else
            for j = 1:nl
                if lines{j}.isCV && j ~= i
                    if j == 1
                        if i ~= nl && i ~= 2
                            noneq = [noneq, j];
                        end
                    elseif j == nl
                        if i ~= nl-1 && i ~= 1
                            noneq = [noneq, j];
                        end
                    else
                        if i ~= j-1 && i ~= j+1
                            noneq = [noneq, j];
                        end
                    end
                end
                if ~lines{j}.isCV
                    dif_angs = abs(lines{i}.angle - lines{j}.angle);
                    if (dif_angs > 0);
                        noneq = [noneq, j];
                    end
                end
            end
        end

        lines{i}.non_eq = noneq;
    end

    % computes the set of convex regions that cover the object
    regions = {};
    k = 1;
    for i = 1:nl
        idx_1 = i-1;
        idx_2 = i+1;

        if i == 1; idx_1 = nl; end;
        if i == nl; idx_2 = 1; end;

        if lines{i}.isCV
            regions{k} = struct('A', [], 'b', []);

            dx1 = lines{idx_1}.v2(1) - lines{idx_1}.v1(1);
            dy1 = lines{idx_1}.v2(2) - lines{idx_1}.v1(2); 
            
            dx2 = lines{idx_2}.v2(1) - lines{idx_2}.v1(1);
            dy2 = lines{idx_2}.v2(2) - lines{idx_2}.v1(2); 

            a1 = dy1;
            b1 = -dx1;
            c1 = a1*lines{idx_1}.v2(1) + b1*lines{idx_1}.v2(2);

            a2 = dy2;
            b2 = -dx2;
            c2 = a2*lines{idx_2}.v2(1) + b1*lines{idx_2}.v2(2);

            regions{k}.A = -[a1,b1;a2,b2];
            regions{k}.b = -[c1;c2];

            j = i+1;
            for l = 1:nl-1
                if j > nl
                    j = j - nl;
                end
                if j == 1; idx_1 = nl; end;
                if j == nl; idx_2 = 1; end;
                if lines{j}.isCV
                    dx1 = lines{idx_1}.v2(1) - lines{idx_1}.v1(1);
                    dy1 = lines{idx_1}.v2(2) - lines{idx_1}.v1(2); 
                    
                    dx2 = lines{idx_2}.v2(1) - lines{idx_2}.v1(1);
                    dy2 = lines{idx_2}.v2(2) - lines{idx_2}.v1(2); 

                    a1 = dy1;
                    b1 = -dx1;
                    c1 = a1*lines{idx_1}.v2(1) + b1*lines{idx_1}.v2(2);

                    a2 = dy2;
                    b2 = -dx2;
                    c2 = a2*lines{idx_2}.v2(1) + b1*lines{idx_2}.v2(2);

                    regions{k}.A = [regions{k}.A;-a1,-b1;-a2,-b2];
                    regions{k}.b = [regions{k}.b;-c1;-c2];
                else
                    break;
                end
                j = j + 1;
            end

            k = k + 1;
        elseif ~lines{idx_1}.isCV && ~lines{idx_2}.isCV
            regions{k} = struct('A', [], 'b', []);

            dx = lines{i}.v2(1) - lines{i}.v1(1);
            dy = lines{i}.v2(2) - lines{i}.v1(2); 

            a = dy;
            b = -dx;
            c = a*lines{i}.v2(1) + b*lines{i}.v2(2);

            regions{k}.A = -[a,b];
            regions{k}.b = -c;

            k = k + 1;
        end
    end

    % checks if polygon is convex
    px = x;
    py = y;
    isConvex = false;
    yet = false;

    numPoints = numSides;
    if numPoints < 4
        isConvex = true;
        yet = true;
    end

    if yet == false
        % % can determine if the polygon is convex based on the direction the angles
        % % turn.  If all angles are the same direction, it is convex.
        v1 = [px(1) - px(end), py(1) - py(end)];
        v2 = [px(2) - px(1), py(2) - py(1)];
        signPoly = sign(det([v1; v2]));

        % % check subsequent vertices
        for k = 2:numPoints-1
            v1 = v2;
            v2 = [px(k+1) - px(k), py(k+1) - py(k)]; 
            curr_signPoly = sign(det([v1; v2]));
            % check that the signs match
            if not (isequal(curr_signPoly, signPoly))
                isConvex = false;
                yet = true;
                break;
            end
        end
    end

    if yet == false
        % % check the last vectors
        v1 = v2;
        v2 = [px(1) - px(end), py(1) - py(end)];
        curr_signPoly = sign(det([v1; v2]));
        if not (isequal(curr_signPoly, signPoly))
            isConvex = false;
        else
            isConvex = true;
        end
    end

    isConvex

    polygons = {};
    iv = 0;
    G = [];

    if isConvex
        % segments the object in convex polygons
        v = [x'; y'];

        polygons{1}.nv = nv;
        polygons{1}.iv = 0;
        polygons{1}.v = v;
        polygons{1}.center = [mean(v(1,:)),mean(v(2,:))]';
        M = 1;
        numvert = nv;
        G = 1;
    else
        dt = delaunayTriangulation(x,y);   % delaunay triangulation 
        points = dt.Points;    % points 
        tri = dt.ConnectivityList;    % nodal connectivity
        x = points(:,1);  
        pX = x(tri);
        y = points(:,2);
        pY = y(tri);

        X = [];
        Y = [];

        for i = 1:length(pX)
            p1 = [pX(i,1);pY(i,1)];
            p2 = [pX(i,2);pY(i,2)];
            p3 = [pX(i,3);pY(i,3)];

            p = (p1+p2+p3)/3;

            if inpolygon(p(1),p(2),px,py)
                X = [X;pX(i,:)];
                Y = [Y;pY(i,:)];
            end
        end

        % segments the object into triangles
        M = size(X,1);
        numvert = M*3;

        for i = 1:M
            polygons{i} = struct('v', [], 'nv', 3, 'iv', iv, 'center', []);
            iv = iv + 3;
            v = [];
            for j = 1:3
                vert = [X(i,j);Y(i,j)];
                v = [v, vert];
            end
            polygons{i}.v = v;
            polygons{i}.center = [mean(v(1,:)),mean(v(2,:))]';
        end

        G = zeros(M,M);
        shareEdge = false;

        % checks for connections
        for i = 1:M
            for j = 1:M
                if i ~= j
                    shareEdge = false;

                    for k = 1:polygons{i}.nv
                        for l = 1:polygons{j}.nv
                            idx1 = k + 1;
                            idx2 = l + 1;

                            if k == polygons{i}.nv; idx1 = 1; end;
                            if l == polygons{j}.nv; idx2 = 1; end;

                            dx1 = polygons{i}.v(1,k) - polygons{j}.v(1,l);
                            dy1 = polygons{i}.v(2,k) - polygons{j}.v(2,l);

                            dx2 = polygons{i}.v(1,idx1) - polygons{j}.v(1,idx2);
                            dy2 = polygons{i}.v(2,idx1) - polygons{j}.v(2,idx2);    

                            if dx1 == 0 && dx2 == 0 && dy1 == 0 && dy2 == 0
                                shareEdge = true;
                                break;
                            end

                            dx1 = polygons{i}.v(1,k) - polygons{j}.v(1,idx2);
                            dy1 = polygons{i}.v(2,k) - polygons{j}.v(2,idx2);

                            dx2 = polygons{i}.v(1,idx1) - polygons{j}.v(1,l);
                            dy2 = polygons{i}.v(2,idx1) - polygons{j}.v(2,l);    

                            if dx1 == 0 && dx2 == 0 && dy1 == 0 && dy2 == 0
                                shareEdge = true;
                                break;
                            end
                        end
                        if shareEdge
                            break;
                        end
                    end

                    shareEdgeVertex = false;
                    for k = 1:polygons{i}.nv
                        for l = 1:polygons{j}.nv

                            dx = polygons{i}.v(1,k) - polygons{j}.v(1,l);
                            dy = polygons{i}.v(2,k) - polygons{j}.v(2,l);

                            if dx == 0 && dy == 0
                                shareEdgeVertex = true;
                                break;
                            end
                        end
                        if shareEdgeVertex
                            break;
                        end
                    end

                    if shareEdge || shareEdgeVertex
                        G(i,j) = 1;
                    end
                end
            end
        end
    end

    % forms shape structure
    shape = struct('polygons',{polygons},'regions',{regions},'lines',{lines},'nv',numvert,'G',G,'isConvex',isConvex);
end