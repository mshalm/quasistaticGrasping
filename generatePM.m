function [suffix] = generatePM(qO, qM, vO, vM, phi, vc, n_hat, suffix)
% GENERATEPM generates .m files for all the necessary geometric functions
% for a given system. inputs:
%   qO: symbolic vector of object generalized coordinates
%   qM: symbolic vector of manipulator generalized coordinates
%   vO: symbolic vector of object generalized velocities vO = qO_dot
%   vM: symbolic vector of manipulator generalized velocites vM = qM_dot
%   phi: symbolic vector for gap functions in terms of qO and qM
%   vc: symbolic matrix containing world-frame contact velocity i in it's
%       i-th column
%   n_hat: symbolic matrix containing world-frame contact normal i in it's
%       i-th column
%   suffix: character array; filename suffix for generated files.

Rz = inline('[cos(t), -sin(t); sin(t), cos(t)]');
nO = numel(qO);
nM = numel(qM);
m = numel(phi);
q = [qO; qM];
v = [vO; vM];
x = [q; v];
n_hat = sym(n_hat);

% contact geometry
JN = sym(zeros(m, nO + nM));
t_hat = Rz(sym(pi)/2)*n_hat;
Ldim = int32(2*m);
JT = sym(zeros(Ldim, numel(v)));

%for all contacts
for i=1:m
  JN(i,:) = n_hat(:,i).'*jacobian(vc(:,i),v);
  JT((2*i-1):(2*i),:) = [t_hat(:,i), -t_hat(:,i)].'*jacobian(vc(:,i),v);
end

% simplify expressions
JN = simplify(JN);
JNO = JN(:, 1:nO);
JNM = JN(:, (nO+1):end);

JT = simplify(JT);
JTO = JT(:, 1:nO);
JTM = JT(:, (nO+1):end);

phi = simplify(phi);

% write .m files
matlabFunction(JNO,'File',['JNO', suffix],'Vars',{q});
matlabFunction(JNM,'File',['JNM', suffix],'Vars',{q});
matlabFunction(JTO,'File',['JTO', suffix],'Vars',{q});
matlabFunction(JTM,'File',['JTM', suffix],'Vars',{q});
matlabFunction(phi,'File',['phi', suffix],'Vars',{q});
end