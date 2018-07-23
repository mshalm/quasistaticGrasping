function [qOp, qMp, z, ctime] = finiteFBTS(qO, qM, NO, NM, LO, LM, u, ...
    A, B, mu, phi, h)
% FINITEFBTS solves finite-feedback timestepping LCP, and computes
% new positions after a time-step of h seconds.
JO = [NO; LO];
JM = [NM; LM];
[m, ~] = size(NO);
o = 1;
E = eye(m);
E = repmat(E,2*o,1);
E = reshape(E, [], size(E,1))';
Z = zeros(m);
ZE = E*0;
z = zeros(m,1);
ze = zeros(2*o*m,1);
M = [JO*A*JO', [Z; E];
     [mu, -E'], Z];
M = M + [JM*B*JM', [Z; ZE];
         [Z, ZE'], Z];
q = [h*JM*u; z] + [phi; ze; z];
tic;
z = pathlcp(M,q);
ctime = toc;
n = z(1:m);
l = z((m + 1):(m * (2 * o + 1)));
f = [n; l];
qOp = qO + A*JO'*f;
qMp = qM + B*JM'*f + h*u;
end

