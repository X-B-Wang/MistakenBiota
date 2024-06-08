% Mistaken Biota at Mistaken Point, Newfoundland
% (Clapham et al., 2003)
% clear all; close all;

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

for kase = [2]
    surfaces(kase)
    curSurface = surfaceMapObj(surfaces(kase));
    N = length(curSurface);
    X = []; Y = [];
    fossilType = [];
    % Input dataset, initialization
    for i = 1:N
        [num, txt, raw] = xlsread('data/' + surfaces(kase) + '.xlsx', fossilNames(curSurface(i)));
        scatter(num(:, 1), num(:, 2), 20, 'MarkerFaceColor', colorMap(curSurface(i), :), 'MarkerEdgeColor', 'none'); hold on;
        X = [X; num(:, 1)]; Y = [Y; num(:, 2)];
        fossilType = [fossilType; curSurface(i) + 0 * num(:, 1)];
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
    % Draw picture of fossil distributions, which is Voronoi
    set(gcf, 'unit', 'centimeters', 'position', [2 2 15 15]);
    h = voronoi(X, Y, DT);
    h(1).MarkerFaceColor = 'none';
    h(1).Marker = 'none';
    h(2).Color = colorSilver;
    legend(fossilNames(curSurface), 'Location', 'eastOutside');
    axis equal; box on;
    saveas(gcf, surfaces(kase), 'svg');
    close all;
    % Surface informations
    surfacialFossils = categorical(fossilType, curSurface, cellstr(fossilNames(curSurface)));
    set(gcf, 'unit', 'centimeters', 'position', [2 2 2 * length(curSurface) 10]);
    h = histogram(surfacialFossils', 'barWidth', 0.5);
    h.Normalization = 'probability';
    saveas(gcf, surfaces(kase) + '_fossils', 'svg');
    close all;
    % Output to Excel
    output = string(N + 1);
    output(1, 1) = N;
    output(1, 2:N+1) = fossilNames(curSurface);
    output(2:N+1, 1) = fossilNames(curSurface);
    output(2:N+1, 2:N+1) = string(rec(curSurface, curSurface));
    actual = rec;
    xlswrite('Results.xlsx', output, surfaces(kase));
    
    % Random distribution
    randomKase = 10000;
    totRec = zeros(randomKase, N, N);
    for k = 1:randomKase
        % Select a triangle by the area-based probabilities
        res = randsrc(totFossil, 1, [1:convN; prob]);
        % Set a random point inside the parallelogram (doubled triangle)
        tt_1 = rand(totFossil, 1); tt_2 = rand(totFossil, 1);
        % If the point locates outside the triangle, make it home
        rev = (tt_1 + tt_2) > 1;
        tt_1 = abs(rev - tt_1); tt_2 = abs(rev - tt_2);
        X = X_(res) .* tt_1 + X_(res + 1) .* tt_2 + convX(1);
        Y = Y_(res) .* tt_1 + Y_(res + 1) .* tt_2 + convY(1);
        
        codeDelaunay
        totRec(k, :, :) = rec(curSurface, curSurface);
    end
    outputMean = string(N + 1);
    outputMean(1, 1) = N;
    outputMean(1, 2:N+1) = fossilNames(curSurface);
    outputMean(2:N+1, 1) = fossilNames(curSurface);
    outputMean(2:N+1, 2:N+1) = string(rec(curSurface, curSurface));
    for i = 1:N
        for j = 1:N
            realC = actual(curSurface(i), curSurface(j));
            curData = totRec(:, i, j);
            curData = curData(:);
            [muHat1, sigmaHat1] = normfit(curData);
            res = (realC - muHat1) / sigmaHat1;
            outputMean(i+1, j+1) = string(res);
        end
    end
    xlswrite('RandomMean.xlsx', outputMean, surfaces(kase));
    
    % Check if the random procedure fits the Gaussian distribution
    % D surface: [2, 9, 10] = Bradgatia, Pectinifrons, Fractofusus
    % E surface: [4, 7, 12] = Charniodiscus, Feather dusters, Fractofusus
    % LMP surface: [1, 2, 10] = Charnia 'A', Charnia 'B', Ostrich feather
    % E surface: exhibit = [2, 3, 4, 7, 9, 10, 12, 14];
    % D surface: exhibit = [2, 3, 9, 10];
    % LMP surface: exhibit = [1, 2, 3, 6, 10, 11];
    if kase == 3
        exhibit = [2, 3, 4, 7, 12, 9, 10, 14];
    elseif kase == 2
        exhibit = [2, 3, 10, 9];
    elseif kase == 5
        exhibit = [1, 2, 3, 6, 11, 10];
    else
        exhibit = 1:N;
    end
    nn = length(exhibit);
    set(gcf, 'unit', 'centimeters', 'position', [0 0 nn*5 nn*3.2]);
    for i = 1:nn
        for j = 1:nn
            if i < j
                continue;
            end
            subplot(nn, nn, j + (nn - i) * nn);
            curData = totRec(:, exhibit(i), exhibit(j));
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
            if j == 1
                if curSurface(exhibit(i)) == 3
                    s = ylabel("\rmCharniid I", 'color', [0 0 0]);
                elseif curSurface(exhibit(i)) == 4
                    ylabel("\rmCharniid II", 'color', [0 0 0]);
                elseif curSurface(exhibit(i)) == 12 || curSurface(exhibit(i)) == 13 || curSurface(exhibit(i)) == 15 || curSurface(exhibit(i)) == 8
                    ylabel('\rm' + fossilNames(curSurface(exhibit(i))), 'color', [0 0 0]);
                else
                    ylabel('\it' + fossilNames(curSurface(exhibit(i))), 'color', [0 0 0]);
                end
            end
            if i == nn
                if curSurface(exhibit(j)) == 3
                    t = title("\rmCharniid I");
                elseif curSurface(exhibit(j)) == 4
                    title("\rmCharniid II");
                elseif curSurface(exhibit(j)) == 12 || curSurface(exhibit(j)) == 13 || curSurface(exhibit(j)) == 15 || curSurface(exhibit(j)) == 8
                    title('\rm' + fossilNames(curSurface(exhibit(j))));
                else
                    title('\it' + fossilNames(curSurface(exhibit(j))));
                end
                %set(gca, 'XAxisLocation', 'top');
            end
            confidence = [0.99 0.95 0.90 0.68];
            confColor = ["#FF0000"; "#7E2F8E"; "#0072BD"; "#77AC30"; "#EDB120"];
            realC = actual(curSurface(exhibit(i)), curSurface(exhibit(j)));
            xx = [min([rmoutliers(curData, 'mean')' G_x realC]), max([rmoutliers(curData, 'mean')' G_x realC])];
            xx(1) = 1.1 * xx(1) - 0.1 * xx(2);
            xx(2) = 12 / 11 * xx(2) - xx(1) / 11;
            if xx(2) == xx(1)
                xx(2) = xx(1) + 1;
            end
            if xx(1) < 0
                xx(1) = 0;
            end
            xlim(xx);
            for conf = 1:length(confidence)
                xLeft = quantile(curData, (1 - confidence(conf)) / 2);
                xRight = quantile(curData, (1 + confidence(conf)) / 2);
                plot([xLeft, xRight], [yy(2) yy(2)] * (1 + conf / 10), 'color', confColor(conf + 1), 'lineWidth', 4); hold on;
                if realC >= 0
                    if realC < xLeft || xRight < realC
                        scatter(realC, yy(2) * 1.5, 'v', 'markerFaceColor', confColor(conf), 'markerEdgeColor', 'k', 'sizeData', 75, 'lineWidth', 0.75); hold on;
                        realC = -1;
                    end
                end
            end
            %plot([quantile(curData, 0.5), quantile(curData, 0.5)], yy(2) * [1.05 1.6], 'Color', "#808A87", 'LineWidth', 1.5);
            if realC >= 0
                scatter(realC, yy(2) * 1.5, 'o', 'markerFaceColor', confColor(end), 'markerEdgeColor', 'k', 'sizeData', 35, 'lineWidth', 0.75); hold on;
            end
            set(gca, 'tickdir', 'out');
            ylim(yy * 1.7);
        end
    end
    exportgraphics(gcf, surfaces(kase) + '.pdf', 'Resolution', 600);
    %print(gcf, surfaces(kase), '-dpng', '-r0');
    %saveas(gcf, surfaces(kase) + '_Gaussian', 'svg');
    close all;
end