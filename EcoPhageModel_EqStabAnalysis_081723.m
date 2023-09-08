%% Phage Eco Model Equilibria/Stability for Alseth et al. 
% Sarah Sundius, 08/17/23

clear; clc;


%% Equilibria + stability analysis for 2 species gLV with phage

syms P A V b11 b12 b21 b22 c rP rA

% Define model
dPdt = rP*P - (b11*P^2 + b12*P*A + c*P*V);
dAdt = rA*A - (b21*P*A + b22*A^2);

eqnP = dPdt == 0;
eqnA = dAdt == 0;

eqns = [eqnP, eqnA];
vars = [P, A];

% Solve for equilibria
[Peq, Aeq] = vpasolve(eqns,vars);
eq1 = [Peq(1), Aeq(1)];
eq2 = [Peq(2), Aeq(2)];
eq3 = [Peq(3), Aeq(3)];
eq4 = [Peq(4), Aeq(4)];

% Calculate Jacobian
model = [dPdt; dAdt];
J = jacobian(model,vars);

% Evaluate Jacobian at each equilibrium
J1 = subs(J,vars,eq1);
J2 = subs(J,vars,eq2);
J3 = subs(J,vars,eq3);
J4 = subs(J,vars,eq4);

% Find eigenvectors/values of the Jacobian
[V1,D1] = eig(J1);
[V2,D2] = eig(J2);
[V3,D3] = eig(J3);
[V4,D4] = eig(J4);

% Calculate tr/det for stability conditions
cond1_1 = simplify(J1(1,1)+J1(2,2));          % < 0
cond2_1 = simplify(J1(1,1)*J1(2,2) - J1(1,2)*J1(2,1));  % > 0

cond1_2 = simplify(J2(1,1)+J2(2,2));          % < 0
cond2_2 = simplify(J2(1,1)*J2(2,2) - J2(1,2)*J2(2,1));  % > 0

cond1_3 = simplify(J3(1,1)+J3(2,2));          % < 0
cond2_3 = simplify(J3(1,1)*J3(2,2) - J3(1,2)*J3(2,1));  % > 0

cond1_4 = simplify(J4(1,1)+J4(2,2));          % < 0
cond2_4 = simplify(J4(1,1)*J4(2,2) - J4(1,2)*J4(2,1));  % > 0

