function [p] = markov_trasition(pref_1,pref_2,MatchingM,actual_match_M)

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

for i = 1:3
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(2,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 3 && pref_1(2,i) > pref_1(3,i))
b_2(i) =  1;
else 
    b_2(i) =  0;
end
end

for i = 1:3
if MatchingM(i) == 0 || (MatchingM(i) == 1 && ...
        pref_1(3,i) > pref_1(1,i)) || ...
        (MatchingM(i) == 2 && pref_1(3,i) > pref_1(2,i))
b_3(i) =  1;
else 
    b_3(i) =  0;
end
end

b(:,:,1)= [b_1',b_2',b_3'];

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

blocking_pair = b(:,:,1) & b(:,:,2)';

if ~any(blocking_pair(:) == 1)
 if MatchingM == actual_match_M 
     p = 1;
 else
p = 0;
 end
else
index= 1;
for  i = 1:3
for  j = 1:3
    if blocking_pair(i,j) == 1
        new_match_M = MatchingM;
        [row, column] = find(MatchingM == j);
        new_match_M(column)= 0;
        new_match_M(i)= j;
        new_match_M_cell(index,:) = new_match_M;
        index = index + 1;
    end
end
end

new_match_F_cell = new_match_M_cell;

for i = 1:index-1
    for j = 1:3
    if ~any(new_match_M_cell(i,:) == j)
new_match_F_cell(i,j) = 0 ;
    else 
    [~,column] = find(new_match_M_cell(i,:) == j);
    new_match_F_cell(i,j) = column;
    end
    end
end

%%% Vemos si coincide un match 
coincidence = ismember(new_match_M_cell, actual_match_M, 'rows');
 n_coincidence = sum(coincidence);

p =  n_coincidence/(length(new_match_M_cell(:,1)));

end
end
