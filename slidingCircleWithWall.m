function slidingCircleWithWall(save_video)
% SLIDINGCIRCLEWITHWALL generates figures and video for a a rolling
% against a wall example.
if (nargin < 1)
   save_video = false; 
end
syms('qOx','qOy', 'qOth', 'qMx1', 'qMx2', 'qMy1', 'qMy2', 'real');
syms('vOx', 'vOy', 'vOth', 'vMx1', 'vMx2', 'vMy1', 'vMy2', 'real');
qO = [qOx qOy qOth]';
qM = [qMx1 qMy1 qMx2 qMy2]';
vO = [vOx vOy vOth]';
vM = [vMx1 vMy1 vMx2 vMy2]';

radius = 1;

S = inline('[0 -t; t 0]');
renorm = @(t) sqrt(t(1)^2 + t(2)^2);

qM1c = [qMx1 qMy1]';
qM2c = [qMx2 qMy2]';
qOc = [qOx qOy]';
vM1c = [vMx1 vMy1]';
vM2c = [vMx2 vMy2]';
vOc = [vOx vOy]';

phi_1 = renorm(qM1c-qOc) - radius;
phi_2 = renorm(qM2c-qOc) - radius;
phi_3 = qOy - radius;
phi_4 = qMy1;
phi_5 = qMy2;

nhat_1 = (qM1c - qOc)/renorm(qM1c-qOc);
nhat_2 = (qM2c - qOc)/renorm(qM2c-qOc);
nhat_3 = sym([0; -1]);
nhat_4 = sym([0; 1]);
nhat_5 = sym([0; 1]);

vc_1 = vM1c - (vOc + S(vOth)*radius*nhat_1);
vc_2 = vM2c - (vOc + S(vOth)*radius*nhat_2);
vc_3 = -(vOc + S(vOth)*radius*nhat_3);
vc_4 = vM1c;
vc_5 = vM2c;

phi = [phi_1 phi_2 phi_3 phi_4 phi_5]';
nhat = [nhat_1 nhat_2 nhat_3 nhat_4 nhat_5];
vc = [vc_1 vc_2 vc_3 vc_4 vc_5];
suffix = '_tfw';

generatePM(qO, qM, vO, vM, phi, vc, nhat, suffix);


%% colors
r = [1 0 0 1];
g = [0 1 0 1];
b = [0 0 1 1];
y = [.9 .9 0 1];
bl = [0 0 0 1];
wh = [1 1 1 1];
lgr = [0.7 0.7 0.7 1];


%% squeeze into wall
qO0w = [-2, 2, 0]';
fstart = [-4; 4];
del = [0.25; 0.25];
qM0w = [fstart + del; fstart - del];
uw = @(zz,t) [1, -1, 1, -1]';


videoskip = 2;
videoname = 'wallroll.avi';

h = 0.025;
nsteps = ceil(10/h);
geoFun = @(q0, qM) scWallGeo(q0, qM, suffix);

simfun = @finiteFBTS;

paramFun = @() scWallParams(0);
figure(1);
nfin = 14;
nplot = 14;
t_n = cell(1,nfin);
qO_n = cell(1,nfin);
qM_n = cell(1,nfin);
cols_n = cell(1,nfin);
thickness_n = cell(1,nfin);
[t_n{1}, qO_n{1}, qM_n{1}] = simulateQuasiStatics(simfun, geoFun, uw, paramFun, qO0w, qM0w, nsteps, h);
cols_n{1} = {bl,bl,bl};
thickness_n{1} = 5;
plotTrajectory(qO_n{1}, qM_n{1}, @() figure(3), cols_n{1}, thickness_n{1}, true);
for n=2:nfin
    c_n= 2^(-n+4);
    paramFun = @() scWallParams(c_n);
    [t_n{n}, qO_n{n}, qM_n{n}] = simulateQuasiStatics(simfun, geoFun, uw, paramFun, qO0w, qM0w, nsteps, h);
    figure(1);
    hue = (pi/2)*(n/nplot);
    curcol = y*cos(hue) + b*sin(hue);
    curcol(4) = 1;
    cols_n{n} = {curcol, curcol, curcol};
    thickness_n{n} = 3;
    if (n <= nplot)
        plotTrajectory(qO_n{n}, qM_n{n}, @() figure(3), cols_n{n}, thickness_n{n}, true);
    end
end
if (save_video)
    generateVideo(t_n(1:nplot), qO_n(1:nplot), qM_n(1:nplot), @() figure(4), cols_n(1:nplot), thickness_n(1:nplot), videoskip, videoname);
end

end
    
    
    function [mu, A, B, fbscale] = scWallParams(fbscale)
        A = eye(3);
        B = eye(4);
        mu = diag([1 1 1 .5 .5]);
    end
    
    function [phi, JNO, JNM, JTO, JTM] = scWallGeo(qO, qM, suffix)
        q = [qO; qM];
        phifun = str2func(['phi',suffix]);
        phi = phifun(q);
        if (nargout > 1)
            JNOfun = str2func(['JNO',suffix]);
            JNO = JNOfun(q);
            JNMfun = str2func(['JNM',suffix]);
            JNM = JNMfun(q);
            JTOfun = str2func(['JTO',suffix]);
            JTO = JTOfun(q);
            JTMfun = str2func(['JTM',suffix]);
            JTM = JTMfun(q);
        end
    end
    
    function generateVideo(t, O, M, fselect, cols, thickness, frameskip, filename)
        v = VideoWriter(filename);
        open(v);
        for i=1:numel(t)
            t{i} = [t{i}(1:frameskip:end), t{i}(1)];
            O{i} = O{i}(:,1:frameskip:end);
            M{i} = M{i}(:,1:frameskip:end);
        end
        for i=1:(numel(t{1})-1)
            tic
            for n=1:numel(t)
                plotTrajectory(O{n}(:,1:i), M{n}(:,1:i), fselect, cols{n}, thickness{n}, true);
            end
            plotStyleReel();
            frame = getframe(gcf);
            writeVideo(v,frame);
            hold off
            clf
            
        end
        close(v);
    end

    function plotPose(O, M, fselect, cols, thickness, plotFingers)
        if (nargin < 6)
           plotFingers = false; 
        end
        fselect();
        
        r_vec = 1;
        th_vec = 0:pi/50:2*pi;
        xunit = r_vec * cos(th_vec) + O(1);
        yunit = r_vec * sin(th_vec) + O(2);
        plot(xunit, yunit,'LineWidth',thickness,'Color',cols{1});
        axis equal;
        hold on;
        xth = [O(1), O(1) + cos(O(3))];
        yth = [O(2), O(2) + sin(O(3))];
        plot(xth, yth,'LineWidth',thickness,'Color',cols{2});
        if (plotFingers)
            plot(M(1),M(2),'o','MarkerEdgeColor',cols{3}(1:3),'MarkerFaceColor',cols{3}(1:3));
            plot(M(3),M(4),'o','MarkerEdgeColor',cols{3}(1:3),'MarkerFaceColor',cols{3}(1:3));
        end
        patch([15 15 (-15) (-15)],[(-1) 0 0 (-1)],[.5 .5 .5]);
    end
    
    
    function plotTrajectory(O, M, fselect, cols, thickness, plotFingers)
        if (nargin < 6)
           plotFingers = false; 
        end
        fselect();
        hold on;
        plot(O(1,:),O(2,:),':','LineWidth',thickness,'Color',cols{1});
        hold on
        plotPose(O(:,end), M(:,end), fselect, cols, thickness, plotFingers);
        plotStyleTrajectory();
    end
    
    function plotStyleTrajectory
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$x_W$','Interpreter','latex');
        ylabel('$y_W$','Interpreter','latex');
        axis equal
        xl = 7;
        yl = 7;
        xlim(xl*[-1,1]);
        ylim(yl*[-1,1] + 4);
    end
    
    function plotStyleReel
        axis equal
        xl = 9;
        yl = 9;
        xlim(xl*[-1,1]);
        ylim(yl*[-1,1] + 4);
        axis off
        set(gcf,'color','w');
    end
