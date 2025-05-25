
function [ct, cst] = HI_substructure(ct, cst, resultGUI, doseThreshold, structureColor, structureColor2)   %, HU_highlight)
    % Inputs:
    %   ct: MatRad CT structure
    %   cst: MatRad CST structure
    %   resultGUI: MatRad resultGUI (contains dose cube)
    %   doseThreshold: Dose threshold (Gy) for hotspots (default: 5)
    %   structureColor: RGB color for hotspot ROI (default: [1 0 0] = red)
    %   HU_highlight: HU value to assign hotspots in CT (default: 1000)

    %% Set defaults
    if ~exist('doseThreshold','var') || isempty(doseThreshold)
        doseThreshold = 5;
    end
    if ~exist('structureColor','var') || isempty(structureColor)
        structureColor = [1 0 0];
    end
    if ~exist('structureColor2','var') || isempty(structureColor2)
        structureColor2 = [0 1 0];
    end
    % if ~exist('HU_highlight','var') || isempty(HU_highlight)
    %     HU_highlight = 1000;
    % end

    %% Find PTV
    ixPTV = find(strcmp(cst(:,2), 'PTV'));
    if isempty(ixPTV)
        error('PTV not found in CST');
    end
    ptvVoxels = cst{ixPTV,4}{1};

    %% Create initial hotspot mask (threshold only)
    doseCube = resultGUI.physicalDose;
    thresholdMask = false(ct.cubeDim);
    thresholdMask(ptvVoxels) = doseCube(ptvVoxels) > doseThreshold;

    %% Create filled mask (no holes)
    filledMask = imfill(thresholdMask, 'holes'); % Fill all holes
    
    % Restrict to PTV volume
    ptvMask = false(ct.cubeDim);
    ptvMask(ptvVoxels) = true;
    filledMask = filledMask & ptvMask;

    %% Create substructures
    % Original threshold-based structure
    ixThreshold = size(cst, 1) + 1;
    cst{ixThreshold,1} = ixThreshold - 1;
    cst{ixThreshold,2} = sprintf('PTV_>%dGy', doseThreshold);
    cst{ixThreshold,3} = 'EXTERNAL';
    cst{ixThreshold,4} = {find(thresholdMask)};
    cst{ixThreshold,5} = struct(...
        'TissueClass', 1, ...
        'alphaX', 0.1, ...
        'betaX', 0.05, ...
        'Priority', 1, ...
        'Visible', 1, ...
        'visibleColor', structureColor, ...
        'sfudOptimization', 0);

    % Filled structure (no holes)
    ixFilled = size(cst, 1) + 1;
    cst{ixFilled,1} = ixFilled - 1;
    cst{ixFilled,2} = sprintf('PTV_>%dGy_filled', doseThreshold);
    cst{ixFilled,3} = 'EXTERNAL';
    cst{ixFilled,4} = {find(filledMask)};
    cst{ixFilled,5} = struct(...
        'TissueClass', 1, ...
        'alphaX', 0.1, ...
        'betaX', 0.05, ...
        'Priority', 2, ...
        'Visible', 1, ...
        'visibleColor', structureColor2, ... % Darker shade
        'sfudOptimization', 0);

    % %% Update CT (optional)
    % if exist('HU_highlight','var')
    %     ct.cubeHU{1}(filledMask) = HU_highlight;
    % end
end


















































% function cst = HI_substructure(ct, cst, resultGUI, doseThreshold, structureColor, HU_highlight)
%     % Inputs:
%     %   ct: MatRad CT structure
%     %   cst: MatRad CST structure
%     %   resultGUI: MatRad resultGUI (contains dose cube)
%     %   doseThreshold: Dose threshold (Gy) for hotspots (default: 5)
%     %   structureColor: RGB color for hotspot ROI (default: [1 0 0] = red)
%     %   HU_highlight: HU value to assign hotspots in CT (default: 1000)
% 
%     %% Set defaults
%     if ~exist('doseThreshold','var') || isempty(doseThreshold)
%         doseThreshold = 5; % 5 Gy threshold
%     end
%     if ~exist('structureColor','var') || isempty(structureColor)
%         structureColor = [1 0 0]; % Red
%     end
%     if ~exist('HU_highlight','var') || isempty(HU_highlight)
%         HU_highlight = 1000; % Bright white in CT
%     end
% 
%     %% Find PTV
%     ixPTV = find(strcmp(cst(:,2), 'PTV'));
%     if isempty(ixPTV)
%         error('PTV not found in CST');
%     end
%     ptvVoxels = cst{ixPTV,4}{1}; % Linear indices
% 
%     %% Extract hotspot voxels (dose > threshold)
%     dosePTV = resultGUI.physicalDose(ptvVoxels);
%     hotspotVoxels = ptvVoxels(dosePTV > doseThreshold);
% 
%     if isempty(hotspotVoxels)
%         warning('No hotspots > %.1f Gy found in PTV.', doseThreshold);
%         return;
%     end
% 
%     %% Update CT to highlight hotspots
%     % ct.cubeHU{1}(hotspotVoxels) = HU_highlight; % Assign new HU value
% 
%     %% Add hotspot ROI to CST
%     ixHotspot = size(cst, 1) + 1;
%     cst{ixHotspot,1} = ixHotspot - 1; % ROI number
%     cst{ixHotspot,2} = sprintf('PTV_>%dGy', doseThreshold); % Name
%     cst{ixHotspot,3} = 'TARGET'; % ROI type
%     cst{ixHotspot,4} = {hotspotVoxels}; % Voxel indices
% 
%     % Set visualization properties
%     cst{ixHotspot,5} = struct(...
%         'TissueClass', 1, ...
%         'alphaX', 0.1, ...
%         'betaX', 0.05, ...
%         'Priority', 1, ... % High priority for visibility
%         'Visible', 1, ...
%         'visibleColor', structureColor, ...
%         'sfudOptimization', 0);
% 
%     %% (Optional) Add a shell around hotspots
%     % hotspotMask = false(ct.cubeDim);
%     % hotspotMask(hotspotVoxels) = true;
%     % shellMask = imdilate(hotspotMask, strel('sphere',3)) & ~hotspotMask;
%     % cst{end+1,4} = {find(shellMask)};
% end
% 
% 
% 
