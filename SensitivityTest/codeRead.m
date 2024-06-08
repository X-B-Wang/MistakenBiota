% Mistaken Biota at Mistaken Point, Newfoundland
% (Clapham et al., 2003)
clear all; close all;
addpath(genpath('../coding'));

colorMap = [[49 130 189]; [107 174 214]; [158 202 225]; [198 219 239]; ...
            [230 85 13]; [253 141 60]; [253 174 107]; [253 208 162]; ...
            [49 163 84]; [116 196 118]; [161 217 155]; [199 233 192]; ...
            [132 60 57]; [173 73 74]; [214 97 107]; [231 150 156]; ...
            [123 65 115]; [165 81 148]; [206 109 189]; [222 158 214]] / 255;
colorSilver = [189 195 199] / 255;
        
% All fossil types, =20
fossilNames = ["Aspidella", "Bradgatia", "Charnia", "Charnia II", "Charniodiscus", ...
               "Discs", "Discs+Stems", "Feather dusters", "Fronds", "Hiemalora", ...
               "Holdfasts", "Ivesia", "Lobate Discs", "Networks", "Ostrich feather", ...
               "Pectinifrons", "Fractofusus", "Spoon Fronds", "Triangles", "Unknown"];
surfaces = ["BCSurface", "DSurface", "ESurface", "GSurface", "LMPSurface", "PCSurface"];
BCSurface = [3; 6; 7; 12; 17; 20];
DSurface = [1; 2; 3; 6; 8; 9; 12; 14; 16; 17; 20];
ESurface = [1; 2; 3; 5; 6; 7; 8; 10; 12; 13; 14; 17; 18; 19; 20];
GSurface = [2; 3; 5; 8; 9; 11; 12; 20];
LMPSurface = [3; 4; 5; 6; 7; 8; 9; 10; 12; 15; 17; 20];
PCSurface = [3; 12; 19; 20];
surfaceMapObj = containers.Map({'BCSurface', 'DSurface', 'ESurface', 'GSurface', 'LMPSurface', 'PCSurface'}, {BCSurface, DSurface, ESurface, GSurface, LMPSurface, PCSurface});

kase = 5;
curSurface = surfaceMapObj(surfaces(kase));
N = length(curSurface);
% Input dataset, initialization
X = []; Y = [];
fossilType = [];
% Input dataset, initialization
fixedFossil = 3;
for i = 1:N
    [num, txt, raw] = xlsread('data/' + surfaces(kase) + '.xlsx', fossilNames(curSurface(i)));
    X = [X; num(:, 1)]; Y = [Y; num(:, 2)];
    fossilType = [fossilType; curSurface(i) + zeros(length(num(:, 1)), 1)];
end
totFossil = length(fossilType);


% Using convex for random distributions
% Make all fossil locates inside the convex hull
% Convex Hull: (1) == (end)
% Initialize a convex hull, and set a record of it
convHull = convhull(X, Y);
convX = X(convHull); convY = Y(convHull);
convN = length(convHull);
% Vectors of each triangle in triangle split
X_ = convX - convX(1);
Y_ = convY - convY(1);
area = zeros(convN, 1);
% Triangle i is composed of points(1, i, i+1), where 2<=i<=convN-2
% (1) == (convN), thus the last triangle is (1, convN-2, convN-1)
area(2:convN-2) = (X_(2:convN-2) .* Y_(3:convN-1) - X_(3:convN-1) .* Y_(2:convN-2)) / 2;
% Using the area of each triangle as the probability of being select
prob = area' / sum(area);

codeDelaunay
actual = rec;

fixingRatio = [0:0.05:0.9 0.91:0.01:1]; % Ratio of spindles that keep their original positions
preference = 0 * fixingRatio;
set(gcf, 'unit', 'centimeters', 'position', [0 0 10 2*length(fixingRatio)]);
for rkase = 1:length(fixingRatio)
    rkase
    randomKase = 1000;
    totRec = zeros(randomKase, N, N);
    for k = 1:randomKase

    %     % Select a triangle by the area-based probabilities
    %     res = randsrc(totFossil, 1, [1:convN; prob]);
    %     % Set a random point inside the parallelogram (doubled triangle)
    %     tt_1 = rand(totFossil, 1); tt_2 = rand(totFossil, 1);
    %     % If the point locates outside the triangle, make it home
    %     rev = (tt_1 + tt_2) > 1;
    %     tt_1 = abs(rev - tt_1); tt_2 = abs(rev - tt_2);
    %     X = X_(res) .* tt_1 + X_(res + 1) .* tt_2 + convX(1);
    %     Y = Y_(res) .* tt_1 + Y_(res + 1) .* tt_2 + convY(1);
    %     % Steady spindles with random others
    %     X(fossilType==fixedFossil) = spindleX; Y(fossilType==fixedFossil) = spindleY;
    % 
    %     codeDelaunay

        % Steady fossils and steady spindles, only others' type would change
        % Shuffle fossil type with given locations
        shuff = 1:totFossil;
        shuffMark = (fossilType~=fixedFossil);
        shuffMark(fossilType==fixedFossil) = (randperm(sum(fossilType==fixedFossil)) > sum(fossilType==fixedFossil) * fixingRatio(rkase));
        shuffArray = shuff(shuffMark);
        shuff(shuffMark) = shuffArray(randperm(sum(shuffMark)));
        newFossil = fossilType(shuff);

        rec = zeros(20);
        for i = 1:20
            for j = 1:20
                rec(i, j) = sum(connections(newFossil==i, newFossil==j), 'all');
            end
        end
        rec = rec + rec';
        for i = 1:20
            rec(i, i) = bitshift(rec(i, i), -1);
        end
        
        totRec(k, :, :) = rec(curSurface, curSurface);
    end
    
    subplot(length(fixingRatio), 1, rkase);
    curData = totRec(:, curSurface==fixedFossil, curSurface==fixedFossil);
    curData = curData(:);
    h = histogram(curData, 'lineWidth', 0.25, 'faceColor', '#4DBEEE', 'faceAlpha', 1); hold on;
    set(gca, 'yGrid', 'on');
    h.Normalization = 'probability';
    while h.NumBins > 20
        fewerbins(h);
    end
    yy = ylim;
    % Gaussian distribution
    [muHat1, sigmaHat1] = normfit(curData);
    tmp = [];
    for bin = 1:h.NumBins
        tmp = [tmp bin+zeros(1, round(randomKase*h.Values(bin)))];
    end
    [muHat, sigmaHat] = normfit(tmp);
    G_x = muHat - 3 * sigmaHat + 6 * sigmaHat * (0:0.01:1);
    G_y = normpdf(G_x, muHat, sigmaHat);
    G_x = (G_x - muHat) / sigmaHat * sigmaHat1 + muHat1;
    plot(G_x, G_y, 'color', '#D95319', 'lineWidth', 1); hold on;
    confidence = [0.99 0.95 0.90];
    confColor = ["#FF0000"; "#7E2F8E"; "#77AC30"; "#EDB120"];
    realC = actual(fixedFossil, fixedFossil);
    preference(rkase) = (realC - muHat1) / sigmaHat1;
    if abs(realC - muHat1) <= 0.0001
        preference(rkase) = 0;
    end
    xlim([1500 1850]);
    set(gca, 'tickdir', 'out', 'yTick', []);
    if rkase ~= length(fixingRatio)
        set(gca, 'xTick', []);
    end
    yl = ylim;
    plot([realC realC], yl, 'color', '#FF00FF'); hold on;
    ylabel(num2str(fixingRatio(rkase)*100, '%.0f'));
end
saveas(gcf, surfaces(kase) + '_' + fossilNames(fixedFossil), 'svg');
close all;

plot(fixingRatio, preference); hold on;
saveas(gcf, surfaces(kase) + '_' + fossilNames(fixedFossil) + '_preference', 'svg');
close all;