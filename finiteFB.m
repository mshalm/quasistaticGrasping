function [qOp, qMp] = finiteFB(qO, qM, NO, NM, LO, LM, u, A, B, mu, phi, h)
JO = [NO; LO];
JM = [NM; LM];
[m, ~] = size(N0);
o = 1;
E = eye(m);
E = repmat(E,2*o,1);
E = reshape(E, [], size(E,1))';
Z = zeros(m);
ZE = E*0;
z = zeros(m,1);
M = [JO*A*JO', [Z; E];
     [mu, -E'], Z];
M = M + [JM*B*JM', [Z; ZE];
         [Z, ZE'], Z];
q = [JM*u; z];

z = pathlcp(M,q);
n = z(1:m);
l = z((m + 1):(m * (2 * o + 1)));
f = [n; l];
qOp = qO + h*A*JO'*f;
qMp = qM + h*(u+B*JM'*f);
end

