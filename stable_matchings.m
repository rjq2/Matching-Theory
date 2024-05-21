function [m_stable,max_m_stable,m_r,f_ex,fragment_2] = stable_matchings(pref_1,pref_2,all_matchings_M)
counts = 1;
count = 1;
fragment_t = [];
fragment_2 = zeros(2,2,2);

for s = 1:length(all_matchings_M(:,1))
    MatchingM = all_matchings_M(s,:);
    MatchingF = zeros(1,3);
    for j = 1:3
    if ~any(MatchingM == j)
MatchingF(j) = 0 ;
    else 
    [~,column] = find(MatchingM == j);
    MatchingF(j) = column;
    end
    end
    
for i = 1:3
if MatchingM(i) == 0 || (MatchingM(i) == 2 && ...
        pref_1(1,i) > pref_1(2,i)) || ...
        (MatchingM(i) == 3 && pref_1(1,i) > pref_1(3,i))
b_1(i) =  1;
else 
    b_1(i) =  0;
end
end

%En el siguiente vector vemos que hombres quieren bloquear con la mujer 2 el
% matching dado
for i = 1:3
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(2,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 3 && pref_1(2,i) > pref_1(3,i))
b_2(i) =  1;
else 
    b_2(i) =  0;
end
end

%En el siguiente vector vemos que hombres quieren bloquear con la mujer 3 el
% matching dado
for i = 1:3
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(3,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 2 && pref_1(3,i) > pref_1(2,i))
b_3(i) =  1;
else 
    b_3(i) =  0;
end
end

% Juntamos todos los vectores para crear una matriz de posibles bloqueos
b(:,:,1)= [b_1',b_2',b_3'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ahora vemos los posibles bloqueos de las mujeres en el matching dado
for i = 1:3
if MatchingF(i) == 0 || (MatchingF(i) == 2 && ...
        pref_2(1,i) > pref_2(2,i)) || ...
        (MatchingF(i) == 3 && pref_2(1,i) > pref_2(3,i))
b_1(i) =  1;
else 
    b_1(i) =  0;
end
end

for i = 1:3
if MatchingF(i) == 0 || (MatchingF(i) == 1 && ...
        pref_2(2,i) > pref_2(1,i)) || ...
        (MatchingF(i) == 3 && pref_2(2,i) > pref_2(3,i))
b_2(i) =  1;
else 
    b_2(i) =  0;
end
end

for i = 1:3
if MatchingF(i) == 0 || (MatchingF(i) == 1 && ...
        pref_2(3,i) > pref_2(1,i)) || ...
        (MatchingF(i) == 2 && pref_2(3,i) > pref_2(2,i))
b_3(i) =  1;
else 
    b_3(i) =  0;
end
end
b(:,:,2)= [b_1',b_2',b_3'];


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
    for i = 1:3
    j=m_stable(k,i);
    pref_m(k,i) = pref_1(j,i);
    end
    end

if length(pref_m(:,1)==1)
m_optimal = pref_m;
else
m_optimal = max(pref_m);
end

for j = 1:3
    h=m_optimal(j) ;
    max_m_stable(j)= find(pref_1(:,j)== h);
end

if length(pref_m(:,1)==1)
f_optimal = pref_m;
else
f_optimal = min(pref_m);
end



for j = 1:3
    h =f_optimal(j);
    max_f_stable(j)= find(pref_1(:,j)== h);
end

%%%% Fragmentos

    for g = 1:length(m_stable(:,1))
    for j = 1:3
    if ~any(m_stable(g,:) == j)
f_stable(g,j) = 0 ;
    else 
    [~,column] = find(m_stable(g,:) == j);
    f_stable(g,j) = column;
    end
    end

        for k = 1:length(f_stable(:,1))
    for i = 1:3
    j=f_stable(k,i);
    pref_f(k,i) = pref_2(j,i);
    end
    end
f_ex=zeros(length(f_stable(:,1)),3,2);
m_r(:,:,1) = m_stable; 
m_r(:,:,2) = pref_m;
f_ex(:,:,1) = f_stable;
f_ex(:,:,2) = pref_f;

count = 1;
b = m_r(:,:,2)==3;
a = f_ex(:,:,2)==3;
for i = 1:3
    for j = 1:length(f_stable(:,1))
        if b(j,i) == 1
            h = m_r(j,i,1);
            if a(j,h) == 1
                fragment_t(count,:) = [i,h];
                count = count+1;
            end
        end
    end
end

 full_set = [1, 2, 3];

% ahora buenos casos
% 1 y 2

for j = 1:length(f_stable(:,1))
m_n = setdiff(full_set,[m_r(j,1,1), m_r(j,2,1)]);
if (m_r(j,1,2)>pref_1(m_n,1))&(m_r(j,2,2)>pref_1(m_n,2))
    m_n = setdiff(full_set,[f_ex(j,m_r(j,1,1),1), m_r(j,m_r(j,2,1),1)]);
   if (f_ex(j,m_r(j,1,1),2)>pref_2(m_n,m_r(j,1,1)))&(f_ex(j,m_r(j,2,1),2)>pref_2(m_n,m_r(j,2,1)))
       fragment_2(:,:,j) = [1,2;m_r(j,1,1),m_r(j,2,1)];
   end
end
end

% 1 y 3 

for j = 1:length(f_stable(:,1))
m_n = setdiff(full_set,[m_r(j,1,1), m_r(j,3,1)]);    
if (m_r(j,1,2)>pref_1(m_n,1))&(m_r(j,3,2)>pref_1(m_n,3))
    m_n = setdiff(full_set,[f_ex(j,m_r(j,1,1),1), m_r(j,m_r(j,3,1),1)]);
   if (f_ex(j,m_r(j,1,1),2)>pref_2(m_n,m_r(j,1,1)))&(f_ex(j,m_r(j,3,1),2)>pref_2(m_n,m_r(j,3,1)))
       fragment_2(:,:,j) = [1,3;m_r(j,1,1),m_r(j,3,1)];
   end
end
end

% 2 y 3

for j = 1:length(f_stable(:,1))
 m_n = setdiff(full_set,[m_r(j,2,1), m_r(j,3,1)]);   
if (m_r(j,2,2)>pref_1(m_n,2))&(m_r(j,3,2)>pref_1(m_n,3))
    m_n = setdiff(full_set,[f_ex(j,m_r(j,2,1),1), m_r(j,m_r(j,3,1),1)]);
   if (f_ex(j,m_r(j,2,1),2)>pref_2(m_n,m_r(j,2,1)))&(f_ex(j,m_r(j,3,1),2)>pref_2(m_n,m_r(j,3,1)))
       fragment_2(:,:,j) = [1,3;m_r(j,1,1),m_r(j,3,1)];
   end
end
end



    



    end


















