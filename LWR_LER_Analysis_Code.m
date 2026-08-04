%% LER/LWR Analysis
%
% This script:
%   1. Loads and calibrates an SEM image.
%   2. Detects edges using the Canny method.
%   3. Allows manual selection and approval of line regions.
%   4. Measures the local line width at each valid y-coordinate.
%   5. Fits a straight nominal line to each detected edge.
%   6. Calculates local edge deviations from the fitted nominal edge.
%   7. Calculates 3-sigma LER and LWR.
%   8. Creates pooled variables and a table for external plotting.

clear;
close all;
clc;

%% Adjustable User Image Settings

imageFile = '10000pC_vertical_20nm.tif';

% Image adjustment settings
brightness = 0;      % Typical range: -100 to +100
contrastFactor = 1;  % Typical range: 0.5 to 3.0

% Image rotation correction in degrees
rotationAngleDegrees = -1;

% Number of pixels removed from the bottom of the SEM image
bottomCropPixels = 80;

% Assumed edge-location uncertainty
pixelErrorPixels = 1;

% Use 3-sigma reporting for LER and LWR
sigmaMultiplier = 3;

% Minimum number of valid y-axis cuts required for fitting an edge
minimumValidCuts = 3;

% Automatically write the edge-deviation table to a CSV file
exportEdgeDeviationCSV = true;


%% Load SEM image

imageData = imread(imageFile);

% Convert to grayscale if necessary
if ndims(imageData) == 3 && size(imageData, 3) == 3
    grayImage = rgb2gray(imageData);
else
    grayImage = imageData;
end


%% Apply brightness and contrast correction

grayImage = double(grayImage);

grayImage = grayImage * contrastFactor + brightness;

grayImage(grayImage > 255) = 255;
grayImage(grayImage < 0) = 0;

grayImage = uint8(grayImage);


%% Display full SEM image for scale-bar measurement

figure;
imshow(grayImage);
title('Full SEM Image Including Scale Bar');

disp('Click the two endpoints of the scale bar.');

[xScale, yScale] = ginput(2);

if length(xScale) ~= 2 || length(yScale) ~= 2
    error('Two scale-bar endpoints were not selected.');
end

scaleBarLengthPixels = sqrt( ...
    (xScale(2) - xScale(1))^2 + ...
    (yScale(2) - yScale(1))^2);

if scaleBarLengthPixels <= 0
    error('The selected scale-bar length must be greater than zero.');
end

knownScaleBarLengthNM = input( ...
    'Enter the known scale-bar length in nanometers: ');

if isempty(knownScaleBarLengthNM) || ...
        ~isnumeric(knownScaleBarLengthNM) || ...
        knownScaleBarLengthNM <= 0
    error('The scale-bar length must be a positive numeric value.');
end


%% Calculate pixel-to-nanometer conversion

pixelToNMConversionFactor = ...
    knownScaleBarLengthNM / scaleBarLengthPixels;

pixelSizeNM = pixelToNMConversionFactor;

fprintf('Pixel size: %.4g nm/pixel\n', pixelSizeNM);

pixelErrorNM = pixelErrorPixels * pixelToNMConversionFactor;

fprintf( ...
    'Pixel error: %.4g nm, assuming %.3g pixel uncertainty\n', ...
    pixelErrorNM, ...
    pixelErrorPixels);


%% Rotate and crop the image

correctedGrayImage = imrotate( ...
    grayImage, ...
    rotationAngleDegrees, ...
    'bilinear', ...
    'crop');

if bottomCropPixels < 0 || ...
        bottomCropPixels >= size(correctedGrayImage, 1)
    error('bottomCropPixels is outside the allowed image range.');
end

if bottomCropPixels == 0
    croppedGrayImage = correctedGrayImage;
else
    croppedGrayImage = ...
        correctedGrayImage(1:end-bottomCropPixels, :);
end


%% Display corrected and cropped image

figure;
imshow(croppedGrayImage);
title('Corrected and Cropped SEM Image');

text( ...
    10, ...
    10, ...
    sprintf('Pixel size: %.4g nm', pixelSizeNM), ...
    'Color', 'red', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', 'white', ...
    'VerticalAlignment', 'top');


%% Apply Canny edge detection

try
    edges = edge(croppedGrayImage, 'canny');
catch ME
    error('Canny edge detection failed: %s', ME.message);
end

if isempty(edges) || ~any(edges(:))
    error(['No edges were detected. Check the image contrast, ' ...
        'brightness, cropping, or edge-detection settings.']);
end


%% Ask user for total lines to analyze

numBoxes = input('Enter the number of boxes to select: ');

if isempty(numBoxes) || ...
        ~isnumeric(numBoxes) || ...
        numBoxes < 1 || ...
        mod(numBoxes, 1) ~= 0
    error('The number of boxes must be a positive integer.');
end


%% Display image for box selection

figure;
imshow(croppedGrayImage);
hold on;
title('Detected Edges: Select Line Regions');


%% Initialize output variables

% One LWR value per accepted line
lineLWR = nan(numBoxes, 1);

% Left and right LER values stored separately
lineLERLeft = nan(numBoxes, 1);
lineLERRight = nan(numBoxes, 1);

% One combined vector containing two LER values per accepted line
lineLER = [];

% Local width vectors
lineWidths = cell(numBoxes, 1);

% Raw edge-position vectors
edgePositionsLeft = cell(numBoxes, 1);
edgePositionsRight = cell(numBoxes, 1);
edgeYPositions = cell(numBoxes, 1);

% Fitted nominal edge vectors
edgeFitsLeft = cell(numBoxes, 1);
edgeFitsRight = cell(numBoxes, 1);

% Local edge-deviation vectors
edgeDeviationsLeft = cell(numBoxes, 1);
edgeDeviationsRight = cell(numBoxes, 1);

% Fit coefficients: [slope, intercept]
edgeFitCoefficientsLeft = cell(numBoxes, 1);
edgeFitCoefficientsRight = cell(numBoxes, 1);

% Box coordinates
selectedBoundingBoxes = nan(numBoxes, 4);

% Track whether each box was accepted
acceptedRegion = false(numBoxes, 1);


%% Select and analyze each line region

for i = 1:numBoxes

    satisfied = false;

    while ~satisfied

        fprintf('\nSelecting box %d of %d.\n', i, numBoxes);
        disp('Click the top-left and bottom-right corners of the box.');

        [xBox, yBox] = ginput(2);

        if length(xBox) ~= 2 || length(yBox) ~= 2
            warning('Two box corners were not selected. Please try again.');
            continue;
        end

        bbox = [ ...
            min(xBox), ...
            min(yBox), ...
            abs(xBox(2) - xBox(1)), ...
            abs(yBox(2) - yBox(1))];

        if bbox(3) < 1 || bbox(4) < 1
            warning('The selected box is too small. Please reselect it.');
            continue;
        end

        hBox = rectangle( ...
            'Position', bbox, ...
            'EdgeColor', 'g', ...
            'LineWidth', 1);

        % Crop edge image to selected region
        roi = imcrop(edges, bbox);

        if isempty(roi)
            warning('The selected ROI is empty.');
            delete(hBox);
            continue;
        end

        % Find all detected edge pixels in the ROI
        [roiY, roiX] = find(roi);

        if length(roiY) <= 1

            warning('Not enough edge points were detected in this region.');

            userResponse = input( ...
                'Do you want to reselect this box? (y/n): ', ...
                's');

            if strcmpi(userResponse, 'y') || ...
                    strcmpi(userResponse, 'yes')

                delete(hBox);
                disp('Reselecting box.');

            else
                satisfied = true;
                delete(hBox);
                fprintf('Skipping region %d.\n', i);
            end

            continue;
        end


        %% Sort edge points by y-coordinate

        [sortedY, sortIndex] = sort(roiY);
        sortedX = roiX(sortIndex);

        uniqueY = unique(sortedY);


        %% Initialize temporary vectors for this ROI

        widthsPixels = [];
        leftEdgeXPixels = [];
        rightEdgeXPixels = [];
        validYPixels = [];

        % Plotting coordinates in full-image coordinates
        leftEdgePlotX = [];
        leftEdgePlotY = [];

        rightEdgePlotX = [];
        rightEdgePlotY = [];

        widthLineX = [];
        widthLineY = [];


        %% Find outer left and right edges at each y-coordinate

        for j = 1:length(uniqueY)

            currentY = uniqueY(j);

            xPositions = sortedX(sortedY == currentY);

            % Sort x explicitly in case find/sort behavior changes
            xPositions = sort(xPositions);

            if length(xPositions) >= 2

                currentLeftX = xPositions(1);
                currentRightX = xPositions(end);

                currentWidth = currentRightX - currentLeftX;

                widthsPixels = [
                    widthsPixels;
                    currentWidth
                ];

                leftEdgeXPixels = [
                    leftEdgeXPixels;
                    currentLeftX
                ];

                rightEdgeXPixels = [
                    rightEdgeXPixels;
                    currentRightX
                ];

                validYPixels = [
                    validYPixels;
                    currentY
                ];

                % Coordinates used only to display accepted edge points
                leftEdgePlotX = [
                    leftEdgePlotX;
                    bbox(1) + currentLeftX
                ];

                leftEdgePlotY = [
                    leftEdgePlotY;
                    bbox(2) + currentY
                ];

                rightEdgePlotX = [
                    rightEdgePlotX;
                    bbox(1) + currentRightX
                ];

                rightEdgePlotY = [
                    rightEdgePlotY;
                    bbox(2) + currentY
                ];

                widthLineX = [
                    widthLineX;
                    bbox(1) + currentLeftX;
                    bbox(1) + currentRightX;
                    NaN
                ];

                widthLineY = [
                    widthLineY;
                    bbox(2) + currentY;
                    bbox(2) + currentY;
                    NaN
                ];
            end
        end


        %% Check number of valid cuts

        if length(validYPixels) < minimumValidCuts

            warning( ...
                ['Only %d valid y-axis cuts were found. At least %d ' ...
                'are required.'], ...
                length(validYPixels), ...
                minimumValidCuts);

            userResponse = input( ...
                'Do you want to reselect this box? (y/n): ', ...
                's');

            if strcmpi(userResponse, 'y') || ...
                    strcmpi(userResponse, 'yes')

                delete(hBox);
                disp('Reselecting box.');

            else
                satisfied = true;
                delete(hBox);
                fprintf('Skipping region %d.\n', i);
            end

            continue;
        end


        %% Display detected edge positions and width cuts

        hLeft = plot( ...
            leftEdgePlotX, ...
            leftEdgePlotY, ...
            'r.', ...
            'MarkerSize', 8);

        hRight = plot( ...
            rightEdgePlotX, ...
            rightEdgePlotY, ...
            'b.', ...
            'MarkerSize', 8);

        hWidths = plot( ...
            widthLineX, ...
            widthLineY, ...
            'y-', ...
            'LineWidth', 1);

        tempPlots = [hLeft, hRight, hWidths];

        drawnow;

        disp('Review the edge-detection results in the figure.');

        userResponse = input( ...
            'Are you satisfied with this box selection? (y/n): ', ...
            's');


        %% Accept or reject selected region

        if strcmpi(userResponse, 'y') || ...
                strcmpi(userResponse, 'yes')

            satisfied = true;
            acceptedRegion(i) = true;
            selectedBoundingBoxes(i, :) = bbox;


            %% Convert coordinates and widths to nanometers

            widthsNM = ...
                widthsPixels * pixelToNMConversionFactor;

            leftEdgeNM = ...
                leftEdgeXPixels * pixelToNMConversionFactor;

            rightEdgeNM = ...
                rightEdgeXPixels * pixelToNMConversionFactor;

            edgeYNM = ...
                validYPixels * pixelToNMConversionFactor;


            %% Save raw width and edge-position data

            lineWidths{i} = widthsNM;

            edgePositionsLeft{i} = leftEdgeNM;
            edgePositionsRight{i} = rightEdgeNM;
            edgeYPositions{i} = edgeYNM;


            %% Calculate LWR

            currentLWR = ...
                sigmaMultiplier * std(widthsNM);

            lineLWR(i) = currentLWR;


            %% Fit nominal straight edges

            % polyfit returns:
            %   coefficient(1) = slope
            %   coefficient(2) = intercept
            leftFitCoefficients = ...
                polyfit(edgeYNM, leftEdgeNM, 1);

            rightFitCoefficients = ...
                polyfit(edgeYNM, rightEdgeNM, 1);

            leftEdgeFitNM = ...
                polyval(leftFitCoefficients, edgeYNM);

            rightEdgeFitNM = ...
                polyval(rightFitCoefficients, edgeYNM);


            %% Calculate local edge deviations

            leftEdgeDeviationNM = ...
                leftEdgeNM - leftEdgeFitNM;

            rightEdgeDeviationNM = ...
                rightEdgeNM - rightEdgeFitNM;


            %% Save fit and residual data

            edgeFitsLeft{i} = leftEdgeFitNM;
            edgeFitsRight{i} = rightEdgeFitNM;

            edgeFitCoefficientsLeft{i} = leftFitCoefficients;
            edgeFitCoefficientsRight{i} = rightFitCoefficients;

            edgeDeviationsLeft{i} = leftEdgeDeviationNM;
            edgeDeviationsRight{i} = rightEdgeDeviationNM;


            %% Calculate detrended LER

            currentLERLeft = ...
                sigmaMultiplier * std(leftEdgeDeviationNM);

            currentLERRight = ...
                sigmaMultiplier * std(rightEdgeDeviationNM);

            lineLERLeft(i) = currentLERLeft;
            lineLERRight(i) = currentLERRight;


            %% Print regional results

            fprintf( ...
                'Region %d: Mean line width = %.4g nm\n', ...
                i, ...
                mean(widthsNM));

            fprintf( ...
                'Region %d: LWR (%g-sigma) = %.4g nm\n', ...
                i, ...
                sigmaMultiplier, ...
                currentLWR);

            fprintf( ...
                'Region %d: Left-edge LER (%g-sigma) = %.4g nm\n', ...
                i, ...
                sigmaMultiplier, ...
                currentLERLeft);

            fprintf( ...
                'Region %d: Right-edge LER (%g-sigma) = %.4g nm\n', ...
                i, ...
                sigmaMultiplier, ...
                currentLERRight);

            fprintf( ...
                'Region %d: Number of valid y-axis cuts = %d\n', ...
                i, ...
                length(widthsNM));

        else

            % Delete rejected box and temporary edge plots
            delete(tempPlots);
            delete(hBox);

            disp('Reselecting box.');
        end
    end
end


%% Remove skipped regions from scalar result vectors

validRegionIndices = find(acceptedRegion);

validLineLWR = lineLWR(acceptedRegion);
validLineLERLeft = lineLERLeft(acceptedRegion);
validLineLERRight = lineLERRight(acceptedRegion);

% One column containing every valid edge-specific LER value
lineLER = [
    validLineLERLeft;
    validLineLERRight
];


%% Check whether any regions were accepted

if isempty(validRegionIndices)
    hold off;
    error('No valid regions were accepted. No results were calculated.');
end


%% Calculate average width for each accepted line

avgWidths = nan(numBoxes, 1);

for i = validRegionIndices'
    avgWidths(i) = mean(lineWidths{i});
end

disp(' ');
disp('Average width of each accepted line:');

for i = validRegionIndices'
    fprintf( ...
        'Region %d: Mean width = %.4g nm\n', ...
        i, ...
        avgWidths(i));
end


%% Pool local width measurements

validWidthCells = lineWidths(acceptedRegion);

allWidths = vertcat(validWidthCells{:});

avgWidthAllLines = mean(allWidths);
stdWidthAllLines = std(allWidths);
numberOfWidthMeasurements = length(allWidths);

fprintf('\nPooled local-width results:\n');
fprintf('Average local width: %.4g nm\n', avgWidthAllLines);
fprintf('Standard deviation of local widths: %.4g nm\n', ...
    stdWidthAllLines);
fprintf('Total number of local-width measurements: %d\n', ...
    numberOfWidthMeasurements);


%% Pool local edge deviations

validLeftDeviationCells = ...
    edgeDeviationsLeft(acceptedRegion);

validRightDeviationCells = ...
    edgeDeviationsRight(acceptedRegion);

allLeftEdgeDeviationsNM = ...
    vertcat(validLeftDeviationCells{:});

allRightEdgeDeviationsNM = ...
    vertcat(validRightDeviationCells{:});

allEdgeDeviationsNM = [
    allLeftEdgeDeviationsNM;
    allRightEdgeDeviationsNM
];


%% Calculate descriptive statistics for pooled deviations

meanLeftEdgeDeviationNM = mean(allLeftEdgeDeviationsNM);
stdLeftEdgeDeviationNM = std(allLeftEdgeDeviationsNM);

meanRightEdgeDeviationNM = mean(allRightEdgeDeviationsNM);
stdRightEdgeDeviationNM = std(allRightEdgeDeviationsNM);

meanAllEdgeDeviationsNM = mean(allEdgeDeviationsNM);
stdAllEdgeDeviationsNM = std(allEdgeDeviationsNM);

fprintf('\nLocal edge-deviation data:\n');

fprintf( ...
    'Number of left-edge deviation measurements: %d\n', ...
    length(allLeftEdgeDeviationsNM));

fprintf( ...
    'Number of right-edge deviation measurements: %d\n', ...
    length(allRightEdgeDeviationsNM));

fprintf( ...
    'Total pooled edge-deviation measurements: %d\n', ...
    length(allEdgeDeviationsNM));

fprintf( ...
    'Mean pooled edge deviation: %.4g nm\n', ...
    meanAllEdgeDeviationsNM);

fprintf( ...
    'Standard deviation of pooled edge deviations: %.4g nm\n', ...
    stdAllEdgeDeviationsNM);


%% Calculate summary LWR results

meanLWR = mean(validLineLWR);
stdLWR = std(validLineLWR);

fprintf('\nLWR summary across accepted regions:\n');

fprintf( ...
    'Mean %g-sigma LWR: %.4g nm\n', ...
    sigmaMultiplier, ...
    meanLWR);

fprintf( ...
    'Standard deviation of line-specific LWR values: %.4g nm\n', ...
    stdLWR);


%% Calculate summary LER results

meanLER = mean(lineLER);
stdLER = std(lineLER);

fprintf('\nLER summary across accepted edges:\n');

fprintf( ...
    'Mean %g-sigma LER: %.4g nm\n', ...
    sigmaMultiplier, ...
    meanLER);

fprintf( ...
    'Standard deviation of edge-specific LER values: %.4g nm\n', ...
    stdLER);


%% Create line-level summary table

lineSummaryTable = table( ...
    validRegionIndices, ...
    avgWidths(acceptedRegion), ...
    validLineLWR, ...
    validLineLERLeft, ...
    validLineLERRight, ...
    'VariableNames', { ...
        'Region', ...
        'MeanWidthNM', ...
        'LWR_3SigmaNM', ...
        'LeftLER_3SigmaNM', ...
        'RightLER_3SigmaNM' ...
    });


%% Create long-format edge-deviation table

regionColumn = [];
edgeSideColumn = strings(0, 1);
yPositionColumnNM = [];
rawEdgePositionColumnNM = [];
fittedEdgePositionColumnNM = [];
edgeDeviationColumnNM = [];

for i = validRegionIndices'

    numberOfLeftPoints = ...
        length(edgeDeviationsLeft{i});

    regionColumn = [
        regionColumn;
        repmat(i, numberOfLeftPoints, 1)
    ];

    edgeSideColumn = [
        edgeSideColumn;
        repmat("Left", numberOfLeftPoints, 1)
    ];

    yPositionColumnNM = [
        yPositionColumnNM;
        edgeYPositions{i}
    ];

    rawEdgePositionColumnNM = [
        rawEdgePositionColumnNM;
        edgePositionsLeft{i}
    ];

    fittedEdgePositionColumnNM = [
        fittedEdgePositionColumnNM;
        edgeFitsLeft{i}
    ];

    edgeDeviationColumnNM = [
        edgeDeviationColumnNM;
        edgeDeviationsLeft{i}
    ];


    numberOfRightPoints = ...
        length(edgeDeviationsRight{i});

    regionColumn = [
        regionColumn;
        repmat(i, numberOfRightPoints, 1)
    ];

    edgeSideColumn = [
        edgeSideColumn;
        repmat("Right", numberOfRightPoints, 1)
    ];

    yPositionColumnNM = [
        yPositionColumnNM;
        edgeYPositions{i}
    ];

    rawEdgePositionColumnNM = [
        rawEdgePositionColumnNM;
        edgePositionsRight{i}
    ];

    fittedEdgePositionColumnNM = [
        fittedEdgePositionColumnNM;
        edgeFitsRight{i}
    ];

    edgeDeviationColumnNM = [
        edgeDeviationColumnNM;
        edgeDeviationsRight{i}
    ];
end

edgeDeviationTable = table( ...
    regionColumn, ...
    edgeSideColumn, ...
    yPositionColumnNM, ...
    rawEdgePositionColumnNM, ...
    fittedEdgePositionColumnNM, ...
    edgeDeviationColumnNM, ...
    'VariableNames', { ...
        'Region', ...
        'EdgeSide', ...
        'YPositionNM', ...
        'RawEdgePositionNM', ...
        'FittedEdgePositionNM', ...
        'EdgeDeviationNM' ...
    });


%% Create simple one-column tables for copying into plotting software

allEdgeDeviationsTable = table( ...
    allEdgeDeviationsNM, ...
    'VariableNames', {'EdgeDeviationNM'});

allLeftEdgeDeviationsTable = table( ...
    allLeftEdgeDeviationsNM, ...
    'VariableNames', {'LeftEdgeDeviationNM'});

allRightEdgeDeviationsTable = table( ...
    allRightEdgeDeviationsNM, ...
    'VariableNames', {'RightEdgeDeviationNM'});


%% Export CSV files

if exportEdgeDeviationCSV

    [~, imageBaseName, ~] = fileparts(imageFile);

    edgeDeviationFile = ...
        [imageBaseName, '_local_edge_deviations.csv'];

    pooledDeviationFile = ...
        [imageBaseName, '_pooled_edge_deviations.csv'];

    lineSummaryFile = ...
        [imageBaseName, '_line_summary.csv'];

    writetable( ...
        edgeDeviationTable, ...
        edgeDeviationFile);

    writetable( ...
        allEdgeDeviationsTable, ...
        pooledDeviationFile);

    writetable( ...
        lineSummaryTable, ...
        lineSummaryFile);

    fprintf('\nExported files:\n');
    fprintf('%s\n', edgeDeviationFile);
    fprintf('%s\n', pooledDeviationFile);
    fprintf('%s\n', lineSummaryFile);
end


%% Finish figure

hold off;
