function [vw, phi, JN, JT] = polyGeometry(vb_om, vb_stat, c, th)
%POLYGEOMETRY generates necessary distances and jacobians for
%multi-polygonal system. All polygons are assumed to be convex.
% Arguments are cell arrays indexed by body number.
% This function assumse q = [c{1}; th{1}; c{2}; th{2}; ...]
%   vb_om: body-coordinate vertices of mobile bodies in counter-clockwise order
%   vb_stat: body-coordinate vertices of static bodies in counter-clockwise order
%   c: world-coordinate centers of bodies
%   th: angles between body coordinates and world coordinates
    R = inline('[cos(t), -sin(t); sin(t), cos(t)]');
    epsilon = 1e-8;
    n_om = numel(vb_om);
    n_stat = numel(vb_stat);
    nb = n_om + n_stat;
    nq = 3*n_om;
    
    vb = [vb_om, vb_stat];
    th_stat = cell(1,n_stat);
    c_stat = cell(1,n_stat);
    [th_stat{:}] = deal(0);
    [c_stat{:}] = deal([0;0]);
    th = [th, th_stat];
    c = [c, c_stat];
    
    %calculate world-frame vertices edge normals, and boundary inequality
    vw = cell(1,nb);
    
    n_vert = cell(1,nb);
    
    n_hat = vb;
    t_hat = vb;
    
    b = cell(1,nb);
    for i=1:nb
        [~, n_vert{i}] = size(vb{i});
        vw{i} = R(th{i})*vb{i} + repmat(c{i},[1 n_vert{i}]);
        for j=1:n_vert{i}
            v1 = vw{i}(:,j);
            v2 = vw{i}(:,mod(j,n_vert{i}) + 1);
            t_hat{i}(:,j) = (v2-v1)/norm(v2-v1);
            n_hat{i}(:,j) = -R(pi/2)*t_hat{i}(:,j);
        end
        b{i} = diag(n_hat{i}'*vw{i});
    end
    
    
    
    if (nargout > 1)
        % calculate phi
        n_phi = 0;
        for i=1:n_om
            for j=(i+1):nb
                n_phi = n_phi + n_vert{i} + n_vert{j};
            end
        end
        ncon = 1;
        v_pen = cell(nb,nb);
        v_pen_idx = cell(nb,nb);
        phi = zeros(n_phi,1);
        if (nargout > 2)
           JN = zeros(n_phi, nq);
           JT = zeros(2*n_phi, nq);
        end
        for i=1:n_om
            for j=(i+1):nb
                    % model each contact as penetration depth of a point
                    % on one polygon into another polygon
                    v_pen{i,j} = max(n_hat{j}'*vw{i} - repmat(b{j},[1 n_vert{i}]));
                    
                    v_pen{j,i} = max(n_hat{i}'*vw{j} - repmat(b{i},[1 n_vert{j}]));

                    phi(ncon:(ncon+n_vert{i}+n_vert{j}-1)) = [v_pen{i,j} v_pen{j,i}];
                    
                    % construct normal and tangential jacobians according
                    % to Zhou et al. 2018:
                    % https://doi.org/10.1177/0278364918755536
                    if (nargout > 2)
                       lociter = 0;
                       for id=1:n_vert{i}
                           j_close = v2bQP(vw{i}(:,id),n_hat{j}',b{j});
                           j_depth = n_hat{j}'*j_close - b{j};
                           edge = find(j_depth == max(j_depth));
                           vjac = velocityJacobian(c{j}, j_close, c{i}, vw{i}(:,id));
                           if (j <= n_om)
                            JN(ncon+lociter,qidx(j)) = n_hat{j}(:,edge(1))'*vjac(:,1:3);
                            JT(2*(ncon+lociter)-1,qidx(j)) = t_hat{j}(:,edge(1))'*vjac(:,1:3);
                            JT(2*(ncon+lociter),qidx(j)) = -t_hat{j}(:,edge(1))'*vjac(:,1:3);
                            
                           end
                           if (i <= n_om)
                            JN(ncon+lociter,qidx(i)) = n_hat{j}(:,edge(1))'*vjac(:,4:end);
                            JT(2*(ncon+lociter)-1,qidx(i)) = t_hat{j}(:,edge(1))'*vjac(:,4:end);
                            JT(2*(ncon+lociter),qidx(i)) = -t_hat{j}(:,edge(1))'*vjac(:,4:end);
                           end
                           lociter = lociter+1;
                       end
                       for id=1:n_vert{j}
                           i_close = v2bQP(vw{j}(:,id),n_hat{i}',b{i});
                           i_depth = n_hat{i}'*i_close - b{i};
                           edge = find(i_depth == max(i_depth));
                           vjac = velocityJacobian(c{i}, i_close, c{j}, vw{j}(:,id));
                           if (i <= n_om)
                            JN(ncon+lociter,qidx(i)) = n_hat{i}(:,edge(1))'*vjac(:,1:3);
                            JT(2*(ncon+lociter)-1,qidx(i)) = t_hat{i}(:,edge(1))'*vjac(:,1:3);
                            JT(2*(ncon+lociter),qidx(i)) = -t_hat{i}(:,edge(1))'*vjac(:,1:3);
                           end
                           if (j <= n_om)
                            JN(ncon+lociter,qidx(j)) = n_hat{i}(:,edge(1))'*vjac(:,4:end);
                            JT(2*(ncon+lociter)-1,qidx(j)) = t_hat{i}(:,edge(1))'*vjac(:,4:end);
                            JT(2*(ncon+lociter),qidx(j)) = -t_hat{i}(:,edge(1))'*vjac(:,4:end);
                           end
                           lociter = lociter+1;
                       end
                    end
                ncon = ncon + n_vert{i} + n_vert{j};
            end
        end
    end
    
    
    
    function [p] = v2bQP(v,A,b)
        % Finds closest point in polygon A * x \leq b by solving the QP
        % 
        % min_p1   norm(p1 - v)^2
        %   s.t.   A*p1 \leq b
        f = -2*v;
        H = 2*eye(2);
        options =  optimoptions('quadprog','Display','off');
        p = quadprog(H,f,A,b,[],[],[],[],[],options);
    end
    

    function out = qidx(i)
       out = (1:3) + 3*(i-1);
    end


    function out = velocityJacobian(c1,p1,c2,p2)
        p1x = p1(1);
        p1y = p1(2);
        c1x = c1(1);
        c1y = c1(2);
        p2x = p2(1);
        p2y = p2(2);
        c2x = c2(1);
        c2y = c2(2);
        out1 = [ -1,  0, p1y - c1y;
                     0, -1, c1x - p1x];
        out2 = [1, 0, c2y - p2y;
                   0, 1, p2x - c2x];
        out = [out1 out2];
        
    end
end

