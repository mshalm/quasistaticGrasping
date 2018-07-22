function slidingCircle(save_video)
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

R = inline('[cos(t), -sin(t); sin(t), cos(t)]');
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

nhat_1 = (qM1c - qOc)/renorm(qM1c-qOc);
nhat_2 = (qM2c - qOc)/renorm(qM2c-qOc);

vc_1 = vM1c - (vOc + S(vOth)*radius*nhat_1);
vc_2 = vM2c - (vOc + S(vOth)*radius*nhat_2);

phi = [phi_1 phi_2]';
nhat = [nhat_1 nhat_2];
vc = [vc_1 vc_2];
suffix = '_tfb';

generatePM(qO, qM, vO, vM, phi, vc, nhat, suffix);


%% colors
r = [1 0 0 1];
g = [0 1 0 1];
b = [0 0 1 1];
y = [.9 .9 0 1];
bl = [0 0 0 1];
wh = [1 1 1 1];
lgr = [0.7 0.7 0.7 1];
n_i = 12;



%% even two-finger push
qO0{1} = [0, 3, 0]';
qM0{1} = [0.5, 5, -0.5, 5]';
u{1} = @(zz,t) [0, -1, 0, -1]';

%% lopsided two-finger push
qO0{2} = [0, 5, 0]';
qM0{2} = [0.5, 8, -0.5, 7]';
u{2} = @(zz,t) [0, -2, 0, -1]';

%% semi-circular push
center = [0 0]';
rad = 5;
grip = pi/6;
freq = pi/10;
rb = rad + sin(grip);
rs = rad - sin(grip);
qO0{3} = [center(1) - rad, center(2), 0]';
qM0{3} = [qO0{3}(1:2)+ R(grip)*[0;1]; qO0{3}(1:2) + R(-grip)*[0;1]];
u{3} = @(zz,t) [rb*freq*sin(freq*t),-rb*freq*cos(freq*t),rs*freq*sin(freq*t),-rs*freq*cos(freq*t)]';


videoskip = 2;
videonames = {'evenpush.avi','unevenpush.avi','circlepush.avi'};

h = 0.025;
nsteps = ceil(10/h);
geoFun = @(q0, qM) scGeo(q0, qM, suffix);

simfun = @finiteFBTS;

for i=1:3
    paramFun = @() scParams(0);
    figure(1);
    nfin = 22;
    nplot = 14;
    c_vector = [];
    e_vector = [];
    t_n = cell(1,nfin);
    qO_n = cell(1,nfin);
    qM_n = cell(1,nfin);
    cols_n = cell(1,nfin);
    thickness_n = cell(1,nfin);
    [t_n{1}, qO_n{1}, qM_n{1}] = simulateQuasiStatics(simfun, geoFun, u{i}, paramFun, qO0{i}, qM0{i}, nsteps, h);
    cols_n{1} = {bl,bl,bl};
    thickness_n{1} = 5;
    plotTrajectory(qO_n{1}, qM_n{1}, @() subplot(2,3,i), cols_n{1}, thickness_n{1}, true);
    for n=2:nfin
        c_n= 2^(-n+4);
        c_vector = [c_vector c_n];
        paramFun = @() scParams(c_n);
        [t_n{n}, qO_n{n}, qM_n{n}] = simulateQuasiStatics(simfun, geoFun, u{i}, paramFun, qO0{i}, qM0{i}, nsteps, h);
        e_vector = [e_vector norm(qO_n{n}(:,end)-qO_n{1}(:,end))];
        figure(1);
        hue = (pi/2)*(n/nplot);
        curcol = y*cos(hue) + b*sin(hue);
        curcol(4) = 1;
        cols_n{n} = {curcol, curcol, curcol};
        thickness_n{n} = 3;
        if (n <= nplot)
            plotTrajectory(qO_n{n}, qM_n{n}, @() subplot(2,3,i), cols_n{n}, thickness_n{n}, true);
        end
        subplot(2,3,i+3);
        
        
    end
    plot(-log2(c_vector),log2(e_vector));
    plotStyleConvergence();
    if (save_video)
        generateVideo(t_n(1:nplot), qO_n(1:nplot), qM_n(1:nplot), @() figure(2), cols_n(1:nplot), thickness_n(1:nplot), videoskip, videonames{i});
    end
end


end
    
    function [mu, A, B, fbscale] = scParams(fbscale)
        A = eye(3);
        B = eye(4);
        mu = eye(2);
    end
    
    function [phi, JNO, JNM, JTO, JTM] = scGeo(qO, qM, suffix)
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
            plotStyleVideo();
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
        ylim(yl*[-1,1]);
    end
    
    function plotStyleVideo
        axis equal
        xl = 9;
        yl = 9;
        xlim(xl*[-1,1]);
        ylim(yl*[-1,1]);
        axis off
        set(gcf,'color','w');
    end
    
    function plotStyleConvergence
        set(gca,'TickLabelInterpreter','latex')
        grid on
        xlabel('$-\log (c_i)$','Interpreter','latex');
        ylabel('$\log (e_i)$','Interpreter','latex');
        axis tight;
        axis square;
    end