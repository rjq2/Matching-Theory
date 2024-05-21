clear all 
clc
% Define los números disponibles (matchings)
numbers = [0, 1, 2, 3];
% almacenamos las permutaciones
permutations = [];
% Genera todas las permutaciones de longitud 3 en las que el 0 puede salir
% hasta 3 veces
for i = 1:length(numbers)
    for j = 1:length(numbers)
        for k = 1:length(numbers)
            % Crea la permutación
            permut = [numbers(i), numbers(j), numbers(k)];
            % Verifica las restricciones
            if sum(permut == 0) <= 3 && sum(permut == 1) <= 1 && sum(permut == 2) <= 1 && sum(permut == 3) <= 1
                permutations = [permutations; permut];
            end
        end
    end
end


all_matchings_M = permutations;

pref_f = [3,1,3;1,3,2;2,2,1];
pref_w = [2,3,1;3,2,2;1,1,3];

% another example

pref_f = [3,2,1;1,3,2;2,1,3];
pref_w = [3,2,2;2,3,1;1,1,3];

[time_medium,number_stable_matchings,M,P,B]=markov_pref(pref_f,pref_w,all_matchings_M);

