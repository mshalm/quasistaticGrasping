function [t,qO,qM,ctime] = simulateQuasiStatics(simfun, geoFun, u, ...
    params, qO0, qM0, nsteps, h)
epsilon = 1e-6;
nqO = numel(qO0);
nqM = numel(qM0);
[mu, A, B, fbscale] = params();
qO = zeros(nqO,nsteps+1);
qM = zeros(nqM,nsteps+1);
t = zeros(1,nsteps+1);
ctime = zeros(1,nsteps+1);
qO(:,1) = qO0;
qM(:,1) = qM0;
t(1) = 0;
ctime(1) = 0;
sim_broke = false;
phi_negative = false;
for i=1:nsteps
    tc = t(i);
    qOc = qO(:,i);
    qMc = qM(:,i);
    qc = [qOc; qMc];
    [phistart] = geoFun(qOc, qMc);
    [~,phisort] = sort(phistart);
    num_tot = numel(phistart);
    num_con = nnz(phistart < 1e-3);
    begin = ceil(log2(max([1 num_con])/num_tot));
    [phic,NOc,NMc,LOc,LMc] = geoFun(qOc, qMc);
    uc = u(qc,tc);
    
    for imprep = 1:1
        if (sim_broke)
            qOp = qOc;
            qMp = qMc;
           break 
        end
        for sel=begin:0
            num_cur_con = ceil(2^(sel)*num_tot);
            ri = phisort(1:num_cur_con);
            rid = zeros(1,2*num_cur_con);
            rid(1:2:(2*num_cur_con-1)) = 2*ri - 1;
            rid(2:2:end) = 2*ri;
            NOct = NOc(ri,:);
            NMct = NMc(ri,:);
            LOct = LOc(rid,:);
            LMct = LMc(rid,:);
            phict = phic(ri);
            muct = mu(ri,ri);
            try
                [qOp, qMp, z, ctemp] = simfun(qOc, qMc, NOct, NMct, ... 
                    LOct, LMct, uc, A, fbscale*B, muct, phict, h);
            catch except
                if (~sim_broke)
                    warning(['PATH failed to solve time-step at ' ...
                        't = %.4f. Manipulator and object jamming is ' ...
                        'assumed for the remainder of the trajectory.'],...
                        t(i));
                end
                sim_broke = true;
                break
            end
            ctime(i) = ctime(i) + ctemp;
            phip = geoFun(qOp, qMp);
            ppe = phip(phisort((num_cur_con + 1) : end));
            if (nnz(ppe < -epsilon) == 0)
                break;
            end
        end

        if (sim_broke)
           break 
        end
        phip = geoFun(qOp, qMp);
        if (nnz(phip < -epsilon) == 0)
                break; 
        end
        [~,NOc,NMc,LOc,LMc] = geoFun(qOp, qMp);
    end
    phip = geoFun(qOp, qMp);
    if (nnz(phip < -epsilon) > 0)
        if (~phi_negative)
        warning(['phi below negative tolerance (%.4e) beginning ' ...
             'at t = %.4f with value %.4e'], -epsilon, t(i), min(phip));
        end
        phi_negative = true;
    end
       
    qO(:,i+1) = qOp;
    qM(:,i+1) = qMp;
    t(i+1) = tc + h;
end 



end
