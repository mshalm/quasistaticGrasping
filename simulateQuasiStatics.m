function [t,qO,qM,ctime] = simulateQuasiStatics(simfun, geoFun, u, params, qO0, qM0, nsteps, h)
epsilon = 1e-8;
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
    sim_broke = false;
    for imprep = 1:1
        
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
                [qOp, qMp, z, ctemp] = simfun(qOc, qMc, NOct, NMct, LOct, LMct, uc, A, fbscale*B, muct, phict, h);
            catch except
                disp('SIMULATION FAILED!')
                t(i)
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
        %{
        try
            [qOp, qMp, z] = simfun(qOc, qMc, NOc, NMc, LOc, LMc, uc, A, fbscale*B, mu, phic, h);
        catch except
            disp('SIMULATION FAILED!')
            t(i)
            break
        end
        %}
        if (sim_broke)
           break 
        end
        phip = geoFun(qOp, qMp);
        if (nnz(phip < -epsilon) == 0)
                break; 
        end
        %qOc = qOp;
        %qMc = qMp;
        %qc = [qOc; qMc];
        [~,NOc,NMc,LOc,LMc] = geoFun(qOp, qMp);
    end
    phip = geoFun(qOp, qMp);
    if (nnz(phip < -epsilon) > 0)
        t(i)
        min(phip)
        disp('negative phi!');
        %break;
    end
       
    qO(:,i+1) = qOp;
    qM(:,i+1) = qMp;
    t(i+1) = tc + h;
end 



end
