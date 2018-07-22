function [suffix] = generatePM(qO, qM, vO, vM, phi, vc, n_hat, suffix)

Rz = inline('[cos(t), -sin(t); sin(t), cos(t)]');
nO = numel(qO);
nM = numel(qM);
m = numel(phi);
q = [qO; qM];
v = [vO; vM];
x = [q; v];
n_hat = sym(n_hat);

% contact geometry
JN = sym(zeros(m, nO + nM))%;jacobian(phi,q);
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

%{
writeMatrixFunction(A, [], '', ['A', suffix], );
writeMatrixFunction(B, [], '', ['B', suffix], 0);
writeMatrixFunction(NO, q, 'q', ['NO', suffix], 0);
writeMatrixFunction(NM, q, 'q', ['NM', suffix], 0);
writeMatrixFunction(LO, q, 'q', ['LO', suffix], 0);
writeMatrixFunction(LM, q, 'q', ['LM', suffix], 0);
writeMatrixFunction(phi, q, 'q', ['phi', suffix], 0);
writeMatrixFunction(mu, [], '', ['mu', suffix], 0);
%}
end