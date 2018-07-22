function polygonalPegInHole(save_video)
if (nargin < 1)
   save_video = false; 
end
R = inline('[cos(t), -sin(t); sin(t), cos(t)]');
relpose = @(dqO,dqM) [(dqO(1:2) + R(dqO(3))*dqM(1:2))' (dqO(3)+dqM(3)) (dqO(1:2) + R(dqO(3))*dqM(4:5))' (dqO(3)+dqM(6))]';
tolerance = .01;
O_b = [1, -3; 1, 3; -1, 3; -1, -3]';
M1_b = [1, -1; 1, 1; -1, 1]';
M2_b = [1, -1; 1, 1; -1, 1]';
S1_b = [-1-tolerance, -10; -1-tolerance, 0; -19-tolerance, 0; -19-tolerance, -10]';
S2_b = diag([-1, 1]) * S1_b;

b_om = {O_b, M1_b, M2_b};
b_stat = {S1_b, S2_b};
bod_c = {[0 0]',[0 0]',[0 0]'};
bod_th = {0,0,0};
[~, phi, JN, JT] = polyGeometry(b_om, b_stat, bod_c, bod_th);


%% colors
r = [1 0 0 1];
g = [0 1 0 1];
b = [0 0 1 1];
y = [.9 .9 0 1];
bl = [0 0 0 1];
wh = [1 1 1 1];
lgr = [0.7 0.7 0.7 1];
n_i = 0;
bod_cols = {b,r,r,y,y};


simfun = @finiteFBTS;
geoFun = @(qOc, qMc) PPIHGeoFun(qOc, qMc, b_om, b_stat);
params = @() PPIHParams(.01);

u = @sideSlide;
qO0 = [-2 7 pi/3]';
qM0 = relpose(qO0,[-1.5, 5, 3*pi/4, 1.5, 5, -pi/4]');
h = 0.05;
nsteps = ceil(12/h);
videoskip = 1;
videospeed = 1;
trajectoryskip = 9;
fnum = @(num) @() figure(num);

[t_sim, qO_sim, qM_sim, ctime_sim] = simulateQuasiStatics(simfun, geoFun, u, params, qO0, qM0, nsteps, h);
if (save_video)
    generateVideo(t_sim, qO_sim, qM_sim, b_om, b_stat, fnum(3), bod_cols, videoskip, videospeed, 'peginhole.avi');
end
tsplit = 3.5;
figure(4);
clf;
plotTrajectory(t_sim(1:round(tsplit/h)-1), qO_sim(:,1:round(tsplit/h)), qM_sim(:,1:round(tsplit/h)), b_om, b_stat, @() subplot(1,2,1), bod_cols, trajectoryskip, @(c,t) .4*(c/t)^4);
plotTrajectory(t_sim(round(tsplit/h):end), qO_sim(:,round(tsplit/h):end), qM_sim(:,round(tsplit/h):end), b_om, b_stat, @() subplot(1,2,2), bod_cols, trajectoryskip, @(c,t) .4*(c/t)^4);
    
    function [c, th] = qToCth(q)
        nq = numel(q);
        nb = nq/3;
        c = cell(1,nb);
        th = cell(1,nb);
        for i=1:nb
            c{i} = q([3*i-2,3*i-1]);
            th{i} = q(3*i);
        end 
    end

    function [phi, JNO, JNM, JTO, JTM]  = PPIHGeoFun(qO, qM, b_om, b_stat)
        q = [qO; qM];
        [c, th] = qToCth(q);
        if (nargout < 2)
            [~, phi] = polyGeometry(b_om, b_stat, c, th);
        else
            [~, phi, JN, JT] = polyGeometry(b_om, b_stat, c, th);
            JNO = JN(:,1:3);
            JNM = JN(:,4:end);
            JTO = JT(:,1:3);
            JTM = JT(:,4:end);
        end
    end

    function [mu, A, B, fbscale] = PPIHParams(fbscale)
       mu = eye(58);
       A = eye(3);
       B = eye(6);
    end

    function plotPoly(poly, color)
         poly1 = patch(poly(1,:),poly(2,:),color(1:3));
        set(poly1,'facealpha',color(4));
    end

    function plotPose(b_om, b_stat, qO, qM, cols, fselect, clear)
        q = [qO; qM];
        fselect();
        if (clear)
            hold off
            clf;
        end
        [c, th] = qToCth(q);
        vb = polyGeometry(b_om, b_stat, c, th);
        for i=1:numel(vb)
           plotPoly(vb{i},cols{i});
           hold on
        end
        axis equal;
        plotStyleTrajectory();
    end

    function generateVideo(t, O, M, b_om, b_stat, fselect, cols, skip, speed, filename)
        v = VideoWriter(filename);
        open(v);
        t = [t(1:skip:end), t(1)];
        O = O(:,1:skip:end);
        M = M(:,1:skip:end);
        for i=1:(numel(t)-1)
            tic
            
            plotPose(b_om, b_stat, O(:,i),M(:,i), cols, fselect, 1);
            plotStyleVideo();
            pause((t(i+1)-t(i))/speed - toc);
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v);
    end

    function plotTrajectory(t, O, M, b_om, b_stat, fselect, cols, skip, alpha)
        
        t = [t(1:skip:end), t(1)];
        O = O(:,1:skip:end);
        M = M(:,1:skip:end);
        for i=1:(numel(t)-1)
            cols_j = cols;
            for j=1:numel(b_om)
                cols_j{j}(4) = alpha(i,numel(t)); 
            end
            plotPose(b_om, b_stat, O(:,i),M(:,i), cols_j, fselect, 0);
        end
    end

    function u = sideSlide(qcur, tcur)
        push_pos = [-1-cos(17*pi/16), -2, 17*pi/16];
        qO = qcur(1:3);
        qM = qcur(4:end);
        all_id = ones(1,6);
        th_id = [0 0 1 0 0 1];
        u = zeros(6,1);
        relpose = @(dqO,dqM) [(dqO(1:2) + R(dqO(3))*dqM(1:2))' (dqO(3)+dqM(3)) (dqO(1:2) + R(dqO(3))*dqM(4:5))' (dqO(3)+dqM(6))]';
        if (tcur <= 1)
            K = diag(all_id);
            qM_des = relpose(qO,[-1.5, 5, 3*pi/4, 1.5, 5, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 1.5)
            K = diag(all_id);
            qM_des = relpose(qO,[-1.5, 0, 3*pi/4, 1.5, 0, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 1.75)
            K = diag(all_id);
            qM_des = relpose(qO,[0, 0, 3*pi/4, 0, 0, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 2.5)
            K = diag(all_id);
            qM_des = relpose([-2 7 0]',[0, 0, 3*pi/4, 0, 0, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 3.5)
            K = diag([1 1 3 1 1 3]);
            qM_des = relpose([-2 2.5 0]',[0, 0, 3*pi/4, 0, 0, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 4.5)
            K = diag([1 1 1 1 1 1]);
            qM_des = relpose([0 3.5 0]',[0, 0, 3*pi/4, 0, 0, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 7.5)
            K = diag(all_id);
            qM_des = relpose([0 3 0]',[0, -.5, 3*pi/4, 0, -.5, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 8)
            K = diag(all_id);
            qM_des = relpose(qO,[-1.5, 0, 3*pi/4, 1.5, 0, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 8.5)
            K = diag(all_id);
            qM_des = relpose(qO,[-1.5, 0, 3*pi/4, 1.5, 5, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 9)
            K = diag(all_id);
            qM_des = relpose(qO,[-1.5, 0, 3*pi/4, .5, 5, -pi/4]');
            u = 5*K*(qM_des - qM);
        elseif (tcur <= 100)
            K = diag(all_id);
            qM_des = relpose([0 -3 0]',[-1.5, 0, 3*pi/4, .5, 3+sqrt(2), -pi/4]');
            u = 5*K*(qM_des - qM);
            u(1:3) = [0 0 0]';
        end
        
    end

    function plotStyleTrajectory
        xlim([-10 10]);
        ylim([-5 15]);
        xlabel('$x_W$','Interpreter','latex');
        ylabel('$y_W$','Interpreter','latex');
        set(gca,'TickLabelInterpreter','latex')
    end

    function plotStyleVideo
        xlim([-10 10]);
        ylim([-5 15]);
        axis equal
        axis off
        set(gcf,'color','w');
    end
end

