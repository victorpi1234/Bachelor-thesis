function [ct, cst] = HI_parallelepiped_substructure(ct, cst, resultGUI, doseThreshold, structureColor, includeOriginalSubstruct, manualReferenceDose)
    % Inputs:
    %   ct: MatRad CT structure
    %   cst: MatRad CST structure
    %   resultGUI: MatRad resultGUI (contains dose cube)
    %   doseThreshold: Dose threshold (Gy) for hotspots
    %   structureColor: RGB color for hotspot ROI
    %   includeOriginalSubstruct: true/false to include original dose-based substructure
    %   manualReferenceDose: optional manual override for reference dose (Gy)

    %% Set defaults
    if ~exist('doseThreshold','var') || isempty(doseThreshold)
        doseThreshold = 5;
    end
    if ~exist('structureColor','var') || isempty(structureColor)
        structureColor = [1 0 0];
    end
    if ~exist('includeOriginalSubstruct','var') || isempty(includeOriginalSubstruct)
        includeOriginalSubstruct = true;
    end
    if ~exist('manualReferenceDose','var') || isempty(manualReferenceDose)
        manualReferenceDose = doseThreshold; % Default to doseThreshold if not specified
    end

    %% Find PTV
    ixPTV = find(strcmp(cst(:,2), 'PTV'));
    if isempty(ixPTV)
        error('PTV not found in CST');
    end
    ptvVoxels = cst{ixPTV,4}{1};
    
    %% Create dose-based mask
    doseMask = false(ct.cubeDim);
    doseMask(ptvVoxels) = resultGUI.physicalDose(ptvVoxels) > doseThreshold;
    
    if ~any(doseMask(:))
        warning('No voxels > %.1f Gy found in PTV', doseThreshold);
        return;
    end

    %% Create original dose-based substructure with objectives
    if includeOriginalSubstruct
        ixDose = size(cst, 1) + 1;
        cst{ixDose,1} = ixDose - 1;
        cst{ixDose,2} = sprintf('PTV_>%dGy', doseThreshold);
        cst{ixDose,3} = 'TARGET';
        cst{ixDose,4} = {find(doseMask)};
        cst{ixDose,5} = struct(...
            'TissueClass', 1, ...
            'alphaX', 0.1, ...
            'betaX', 0.05, ...
            'Priority', 1, ...
            'Visible', 1, ...
            'visibleColor', structureColor*0.7, ...
            'sfudOptimization', 0);
        
        % Add dose objective using either manualReferenceDose or doseThreshold
        cst{ixDose,6}{1} = struct(...
            'type', 'square deviation', ...
            'penalty', 100, ...
            'dose', manualReferenceDose, ...
            'EUD', manualReferenceDose, ...
            'parameters', {num2cell(manualReferenceDose)}); 
    end

    %% Create 3D parallelepiped with objectives
    [y,x,z] = ind2sub(ct.cubeDim, find(doseMask));
    
    % Get extremes with margin
    margin = 1; % voxels
    x_min = max(1, min(x)-margin);
    x_max = min(ct.cubeDim(2), max(x)+margin);
    y_min = max(1, min(y)-margin);
    y_max = min(ct.cubeDim(1), max(y)+margin);
    z_min = max(1, min(z)-margin);
    z_max = min(ct.cubeDim(3), max(z)+margin);

    paraMask = false(ct.cubeDim);
    paraMask(y_min:y_max, x_min:x_max, z_min:z_max) = true;

    ixPara = size(cst, 1) + 1;
    cst{ixPara,1} = ixPara - 1;
    cst{ixPara,2} = sprintf('HotspotBox_>%dGy', doseThreshold);
    cst{ixPara,3} = 'TARGET';
    cst{ixPara,4} = {find(paraMask)};
    cst{ixPara,5} = struct(...
        'TissueClass', 1, ...
        'alphaX', 0.1, ...
        'betaX', 0.05, ...
        'Priority', 2, ...
        'Visible', 1, ...
        'visibleColor', structureColor, ...
        'sfudOptimization', 0);
    
    % Add dose objective to parallelepiped
    cst{ixPara,6}{1} = struct(...
        'type', 'square deviation', ...
        'penalty', 100, ...
        'dose', manualReferenceDose, ...
        'EUD', manualReferenceDose, ...
        'parameters', {num2cell(manualReferenceDose)}); 
end