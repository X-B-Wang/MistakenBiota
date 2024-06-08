% Deal with the overlap fossils
% Among the overlapping fossils, the one with smallest #Number is used in Delaunay
% Using a overlap[i]=i array to record the overlapping
% And the connections will be updated by the headmost fossil
% for each overlapping pair fossils (i, j)
% the overlaps of fossil(i) is: 1) found, 2) non-existent
% thus let overlap(j) = overlap(i) in order to direct fossil(j) to its headmost overlapping fossil
overlap = 1:totFossil;
for i = 1:totFossil
    ov = find(X(i+1:end)==X(i) & Y(i+1:end)==Y(i));
    if ov
        overlap(i+ov) = overlap(i);
    end
end
overlapFossils = find(overlap~=1:totFossil);

% Delaunay-Voronoi
DT = delaunay(X, Y);

% Count the connections between delaunay-pairs in triangles
connections = zeros(totFossil);
endDT = length(DT);
for i = 1:endDT
    a = sort(DT(i, :));
    % connections(i, j): i<j
    connections(a(1), a(2)) = 1;
    connections(a(1), a(3)) = 1;
    connections(a(2), a(3)) = 1;
end

% Using overlapping array to update the ignored fossils
% By their headmost ancestors
% Has been tested for accuracy
for i = overlapFossils
    % overlap(i) <= i
    % sum(connections(:, i)) + sum(connections(i, :))
    j = overlap(i);
    connections(1:i-1, i) = [connections(1:j-1, j); 1; connections(j, j+1:i-1)'];
    connections(i, i+1:totFossil) = connections(j, i+1:totFossil);
end

% Sparse,m<=3*n-6(+overlaps)
connections = sparse(connections);
% idx = find(connections==1);
% [idx_1, idx_2] = find(connections==1);

% Convert the fossil neighbours to type neighbours
rec = zeros(20);
s = 0;
for i = 1:20
    for j = i:20
        totC = sum(connections(fossilType==i, fossilType==j), 'all');
        s = s + totC;
        rec(i, j) = totC;
        rec(j, i) = totC;
    end
end
