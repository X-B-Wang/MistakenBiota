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
               "Holdfasts", "Ivesheadiomorphs", "Lobate Discs", "Networks", "Ostrich feather", ...
               "Pectinifrons", "Fractofusus", "Spoon Fronds", "Thectardis", "Unknown"];
surfaces = ["BCSurface", "DSurface", "ESurface", "GSurface", "LMPSurface", "PCSurface"];
BCSurface = [3; 6; 7; 12; 17; 20];
DSurface = [1; 2; 3; 6; 8; 9; 12; 14; 16; 17; 20];
ESurface = [1; 2; 3; 5; 6; 7; 8; 10; 12; 13; 14; 17; 18; 19; 20];
GSurface = [2; 3; 5; 8; 9; 11; 12; 20];
LMPSurface = [3; 4; 5; 6; 7; 8; 9; 10; 12; 15; 17; 20];
PCSurface = [3; 12; 19; 20];
surfaceMapObj = containers.Map({'BCSurface', 'DSurface', 'ESurface', 'GSurface', 'LMPSurface', 'PCSurface'}, {BCSurface, DSurface, ESurface, GSurface, LMPSurface, PCSurface});

kase = 3;
maxDist = 50;
selected = 13;
curSurface = surfaceMapObj(surfaces(kase));
N = length(curSurface);
% Input dataset, initialization
X = []; Y = [];
fossilType = [];
% Input dataset, initialization
for i = 1:N
    [num, txt, raw] = xlsread('data/' + surfaces(kase) + '.xlsx', fossilNames(curSurface(i)));
    X = [X; num(:, 1)]; Y = [Y; num(:, 2)];
    fossilType = [fossilType; curSurface(i) + zeros(length(num(:, 1)), 1)];
end
totFossil = length(fossilType);
% X = X / 10; Y = Y / 10;
% maxDist = maxDist / 10;

% Initialize a convex hull
convHull = convhull(X, Y);
convN = length(convHull);
totArea = polyarea(X(convHull), Y(convHull));
% The lines of the convex hull could be expressed as Ax+By+C=0
% conv_A_B=sqrt(A^2+B^2)
convA = zeros(convN - 1, 1);
convB = zeros(convN - 1, 1);
convC = zeros(convN - 1, 1);
conv_A_B = zeros(convN - 1, 1);
for i = 1:convN - 1
    convA(i) = Y(convHull(i)) - Y(convHull(i + 1));
    convB(i) = X(convHull(i + 1)) - X(convHull(i));
    convC(i) = X(convHull(i)) * Y(convHull(i + 1)) - X(convHull(i + 1)) * Y(convHull(i));
    conv_A_B(i) = sqrt(convA(i) * convA(i) + convB(i) * convB(i));
end


codeDelaunay
selectedArray = (fossilType==selected);
selectedC = connections(selectedArray, selectedArray);
selectedC = selectedC + selectedC';
sumSelected = sum(selectedC);

selectedX = X(selectedArray);
selectedY = Y(selectedArray);

xMin = floor(min(selectedX)); xMax = ceil(max(selectedX));
yMin = floor(min(selectedY)); yMax = ceil(max(selectedY));



peak = zeros(xMax - xMin, yMax - yMin);
peakC = zeros(xMax - xMin, yMax - yMin);
peakM = zeros(xMax - xMin, yMax - yMin);
mk = zeros(xMax - xMin, yMax - yMin);
for i = 1:xMax-xMin
    for j = 1:yMax-yMin
        for k = 1:convN - 1
            x1 = X(convHull(k)); y1 = Y(convHull(k));
            x2 = X(convHull(k+1)); y2 = Y(convHull(k+1));
            x3 = i + xMin; y3 = j + yMin;
            val = x1 * y2 + x2 * y3 + x3 * y1 - y2 * x3 - y1 * x2 - y3 * x1;
            mk(i, j) = mk(i, j) | (val < 0);
        end
        
        for k = 1:length(selectedX)
            xx = selectedX(k) - xMin - i;
            yy = selectedY(k) - yMin - j;
            d = sqrt(xx * xx + yy * yy);
            if d < maxDist
                peak(i, j) = peak(i, j) + 1;
                peakC(i, j) = peakC(i, j) + sumSelected(k);
                if sumSelected(k) > peakM(i, j)
                    peakM(i, j) = sumSelected(k);
                end
            end
        end
    end
end

set(gcf, 'unit', 'centimeters', 'position', [0 0 20 20]);
contourf(peak' / pi / maxDist / maxDist * 10 * 10); colorbar;
exportgraphics(gca, surfaces(kase) + '_' + fossilNames(selected) + '_heatmap' + '.png', 'Resolution', 600);
close all;

set(gcf, 'unit', 'centimeters', 'position', [0 0 20 20]);
contourf(peakC' / pi / maxDist / maxDist * 10 * 10); colorbar;
exportgraphics(gca, surfaces(kase) + '_' + fossilNames(selected) + '_connections' + '.png', 'Resolution', 600);
close all;

set(gcf, 'unit', 'centimeters', 'position', [0 0 20 20]);
contourf(peakM'); colorbar;
exportgraphics(gca, surfaces(kase) + '_' + fossilNames(selected) + '_max_connections' + '.png', 'Resolution', 600);
close all;

set(gcf, 'unit', 'centimeters', 'position', [0 0 20 20]);
averC = peakC ./ peak;
averC(isnan(averC)) = 0;
averC(peak<=1) = 0;
contourf(averC'); colorbar;
exportgraphics(gca, surfaces(kase) + '_' + fossilNames(selected) + '_average_connections' + '.png', 'Resolution', 600);
close all;

set(gcf, 'unit', 'centimeters', 'position', [0 0 20 20]);
colormap('jet'); cmap = colormap;
maxC = max(sumSelected);
for i = 0:maxC
    scatter(selectedX(sumSelected==i)-xMin, selectedY(sumSelected==i)-yMin, 'markerFaceColor', cmap(round(i/maxC*255)+1, :), 'markerEdgeColor', 'none'); hold on;
end
xlim([0 xMax-xMin]); ylim([0 yMax-yMin]);
set(gca, 'xTick', 0:100:xMax-xMin, 'yTick', 0:100:yMax-yMin);
box on; axis equal;
colorbar('Ticks', (0:maxC)/maxC, 'tickLabels', split(num2str(0:maxC)));
exportgraphics(gca, surfaces(kase) + '_' + fossilNames(selected) + '.png', 'Resolution', 600);
close all;
% contourf(peak(1:200, 900:1000)); colorbar;

% h = voronoi(X, Y, DT); hold on;
% idx = sumSelected == 8
% scatter(selectedX(sumSelected==8), selectedY(sumSelected==8), 'markerfacecolor', 'r', 'markeredgecolor', 'none'); hold on;
% for i = 1:1497
%     if selectedC(i, 336) == 1
%         scatter(selectedX(i), selectedY(i), 'markerfacecolor', 'none', 'markeredgecolor', 'k'); hold on;
%     end
% end
