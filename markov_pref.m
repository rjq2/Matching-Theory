function [time_medium,number_stable_matchings,M,P,B] = markov_pref(pref_1,pref_2,all_matchings_M)
M = zeros(34,34);


for i = 1:length(M(:,1))
    for j =1:length(M(:,1))
        M(i,j) = markov_trasition(pref_1,pref_2,all_matchings_M(i,:),all_matchings_M(j,:));
    end
end


filas_con_uno = any(diag(M) == 1, 2);

% Obtiene las diagonales con valor '1' y elimina esas filas de la matriz original
filas_con_uno_valores = M(filas_con_uno, :);
M_2 = M(~filas_con_uno,:);
Q= M(~filas_con_uno,~filas_con_uno);
R = zeros(length(Q),length(M));
for i = 1:length(M)
    if any(M(i,i) == 1)
R(:,i) = [M_2(:,i)];
    end
end

columnas_sin_ceros = ~all(R == 0);

% Elimina las columnas de la diagonal con valor 1 
R = R(:, columnas_sin_ceros);

zero = zeros(length(R(1,:)),length(Q));
identity = eye(length(R(1,:)));

% Forma Canónica

P = [Q,R;zero,identity];

% Matriz N del ppt
N = inv((eye(length(Q))-Q));

% Tiempos esperados de absorcion 

t = N*ones(length(N),1);
time_medium = t(1);
% Probabilidad de Absorción

B = N*R;

number_stable_matchings = length(B(1,:));
% Al parecer la matriz no es ergodica