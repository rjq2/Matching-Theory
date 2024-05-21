function [m_stable,max_m_stable,max_f_stable,pref_m] = stable_matchings4(pref_1,pref_2,all_matchings_M)
counts = 1;
count = 1;
fragment_t = [];
fragment_2 = zeros(2,2,2);

for s = 1:length(all_matchings_M(:,1))
    MatchingM = all_matchings_M(s,:);
    MatchingF = zeros(1,4);
    for j = 1:4
    if ~any(MatchingM == j)
MatchingF(j) = 0 ;
    else 
    [~,column] = find(MatchingM == j);
    MatchingF(j) = column;
    end
    end
for i = 1:4
if MatchingM(i) == 0 || (MatchingM(i) == 2 && ...
        pref_1(1,i) > pref_1(2,i)) || ...
        (MatchingM(i) == 3 && pref_1(1,i) > pref_1(3,i))|| ...
        (MatchingM(i) == 4 && pref_1(1,i) > pref_1(4,i))
b_1(i) =  1;
else 
    b_1(i) =  0;
end
end

%En el siguiente vector vemos que hombres quieren bloquear con la mujer 2 el
% matching dado
for i = 1:4
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(2,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 3 && pref_1(2,i) > pref_1(3,i))|| ...
        (MatchingM(i) == 4 && pref_1(2,i) > pref_1(4,i))
b_2(i) =  1;
else 
    b_2(i) =  0;
end
end

%En el siguiente vector vemos que hombres quieren bloquear con la mujer 3 el
% matching dado
for i = 1:4
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(3,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 2 && pref_1(3,i) > pref_1(2,i))|| ...
        (MatchingM(i) == 4 && pref_1(3,i) > pref_1(4,i))
b_3(i) =  1;
else 
    b_3(i) =  0;
end
end

%En el siguiente vector vemos que hombres quieren bloquear con la mujer 4 el
% matching dado
for i = 1:4
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(4,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 2 && pref_1(4,i) > pref_1(2,i))|| ...
        (MatchingM(i) == 3 && pref_1(4,i) > pref_1(3,i))
b_4(i) =  1;
else 
    b_4(i) =  0;
end
end




% Juntamos todos los vectores para crear una matriz de posibles bloqueos
b(:,:,1)= [b_1',b_2',b_3',b_4'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ahora vemos los posibles bloqueos de las mujeres en el matching dado
for i = 1:4
if MatchingF(i) == 0 || (MatchingF(i) == 2 && ...
        pref_2(1,i) > pref_2(2,i)) || ...
        (MatchingF(i) == 3 && pref_2(1,i) > pref_2(3,i))|| ...
        (MatchingF(i) == 4 && pref_2(1,i) > pref_2(4,i))
    
b_1(i) =  1;
else 
    b_1(i) =  0;
end
end

for i = 1:4
if MatchingF(i) == 0 || (MatchingF(i) == 1 && ...
        pref_2(2,i) > pref_2(1,i)) || ...
        (MatchingF(i) == 3 && pref_2(2,i) > pref_2(3,i))|| ...
        (MatchingF(i) == 4 && pref_2(2,i) > pref_2(4,i))
b_2(i) =  1;
else 
    b_2(i) =  0;
end
end

for i = 1:4
if MatchingF(i) == 0 || (MatchingF(i) == 1 && ...
        pref_2(3,i) > pref_2(1,i)) || ...
        (MatchingF(i) == 2 && pref_2(3,i) > pref_2(2,i))|| ...
        (MatchingF(i) == 4 && pref_2(3,i) > pref_2(4,i))
    
b_3(i) =  1;
else 
    b_3(i) =  0;
end
end

for i = 1:4
if MatchingF(i) == 0 || (MatchingF(i) == 1 && ...
        pref_2(4,i) > pref_2(1,i)) || ...
        (MatchingF(i) == 2 && pref_2(4,i) > pref_2(2,i))|| ...
        (MatchingF(i) == 3 && pref_2(4,i) > pref_2(3,i))
    
b_4(i) =  1;
else 
    b_4(i) =  0;
end
end



b(:,:,2)= [b_1',b_2',b_3',b_4'];


% Intersectamos las matrices de bloqueo, llegado a si a una valor 1 si el
% par de bloqueo es efectivo , si esta matriz esta llena de ceros hemos
% llegado a un matching estable
blocking_pair = b(:,:,1) & b(:,:,2)';


if ~any(blocking_pair==1,"all")
m_stable(count,:) = MatchingM;
count = count + 1;
end

counts = counts + 1;
end

    for k = 1:length(m_stable(:,1))
    for i = 1:4
    j=m_stable(k,i);
    pref_m(k,i) = pref_1(j,i);
    end
    end

if length(pref_m(:,1)==1)
m_optimal = pref_m;
else
m_optimal = max(pref_m);
end

for j = 1:4
    h=m_optimal(j) ;
    max_m_stable(j)= find(pref_1(:,j)== h);
end

if length(pref_m(:,1)==1)
f_optimal = pref_m;
else
f_optimal = min(pref_m);
end

for j = 1:4
    h =f_optimal(j);
    max_f_stable(j)= find(pref_1(:,j)== h);
end 
    end
