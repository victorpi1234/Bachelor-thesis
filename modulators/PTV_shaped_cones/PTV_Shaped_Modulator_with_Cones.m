function [ct, cst] = PTV_Shaped_Modulator_with_Cones(ct, cst, modulatorDepth, yLocation, coneBaseRadius, coneHeight, coneSpacingX, coneSpacingZ, HU, cutout_start, width_margingsYXZ, structureColor, varargin)
    % Parameters with defaults
    if ~exist('modulatorDepth','var') || isempty(modulatorDepth), modulatorDepth = 20; end
    if ~exist('HU','var') || isempty(HU), HU = 132; end
    if ~exist('coneBaseRadius','var') || isempty(coneBaseRadius), coneBaseRadius = 3; end
    if ~exist('coneHeight','var') || isempty(coneHeight), coneHeight = 10; end
    if ~exist('coneSpacingX','var') || isempty(coneSpacingX), coneSpacingX = 2*coneBaseRadius; end
    if ~exist('coneSpacingZ','var') || isempty(coneSpacingZ), coneSpacingZ = 2*coneBaseRadius; end
    if ~exist('cutout_start','var') || isempty(cutout_start), cutout_start = 0; end
    if ~exist('width_margingsYXZ','var') || isempty(width_margingsYXZ), width_margingsYXZ = [5 5 5]; end
    if ~exist('structureColor','var') || isempty(structureColor), structureColor = [0 1 1]; end

    %% Find PTV to determine width
    ixPTV = find(strcmp(cst(:,2), 'PTV'));
    if isempty(ixPTV), error('PTV structure not found in CST'); end
    
    PTV_voxels = cst{ixPTV,4}{1};
    [ptv_y, ptv_x, ptv_z] = ind2sub(ct.cubeDim, PTV_voxels);
    
    %% Calculate dimensions with margin
    x_margin = round(width_margingsYXZ(2)/ct.resolution.x);
    z_margin = round(width_margingsYXZ(3)/ct.resolution.z); 
    x_size = (max(ptv_x) - min(ptv_x)) + 2*x_margin;
    z_size = (max(ptv_z) - min(ptv_z)) + 2*z_margin;
    y_size = round(modulatorDepth/ct.resolution.y);
    
    %% Create modulator structure
    ixModulator = size(cst, 1) + 1; 
    cst{ixModulator,1} = ixModulator-1;
    cst{ixModulator,2} = 'modulator';
    cst{ixModulator,3} = 'EXTERNAL';
    cst{ixModulator,5} = struct('TissueClass',1, 'alphaX',0.1, 'betaX',0.05,...
                               'Priority',2, 'Visible',1, 'visibleColor',structureColor,...
                               'sfudOptimization',0);
    ct.structureColor(ixModulator, :) = structureColor;
    
    %% Create modulator volume
    center = round([ct.cubeDim./2]);
    Modulator_center = [center(1)-yLocation, center(2), center(3)];
    halfSize = [y_size, x_size, z_size]/2;
    
    % Create base mask
    [Y,X,Z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));
    modulatorMask = (abs(X - Modulator_center(2)) <= halfSize(2)) & ...
                   (abs(Y - Modulator_center(1)) <= halfSize(1)) & ...
                   (abs(Z - Modulator_center(3)) <= halfSize(3));
    
    %% Create PTV cutout FIRST
    cutout_start_vox = round(cutout_start/ct.resolution.y);
    top_of_modulator = Modulator_center(1) - halfSize(1);
    new_y = top_of_modulator + cutout_start_vox + (ptv_y - min(ptv_y));
    
    valid = new_y >= 1 & new_y <= ct.cubeDim(1) & ...
            new_y <= (Modulator_center(1) + halfSize(1));
    moved_PTV = sub2ind(ct.cubeDim, new_y(valid), ptv_x(valid), ptv_z(valid));
    
    modulatorMask(moved_PTV) = 0;  % Remove PTV from base FIRST
    
    %% Find actual bottom surface AFTER PTV removal
    % Get the lowest remaining y-value for each (x,z) position in modulator
    bottom_surface = zeros(ct.cubeDim(2), ct.cubeDim(3));
    for x = 1:ct.cubeDim(2)
        for z = 1:ct.cubeDim(3)
            y_values = find(modulatorMask(:,x,z));
            if ~isempty(y_values)
                bottom_surface(x,z) = max(y_values); % Highest y = bottom in MATLAB coordinates
            end
        end
    end
    
    %% Create downward-pointing cones starting from actual bottom surface
    x_positions = Modulator_center(2)-halfSize(2):coneSpacingX:Modulator_center(2)+halfSize(2);
    z_positions = Modulator_center(3)-halfSize(3):coneSpacingZ:Modulator_center(3)+halfSize(3);
    modulatorCones = zeros(ct.cubeDim);
    
    for x_idx = 1:length(x_positions)
        for z_idx = 1:length(z_positions)
            x = x_positions(x_idx);
            z = z_positions(z_idx);
            
            % Get actual bottom y-value at this (x,z) position
            y_base = bottom_surface(x,z);
            if y_base == 0, continue; end % Skip if no base at this position
            
            % Grow downward (increasing y) from the actual bottom
            for y_offset = 0:coneHeight
                y = y_base + y_offset;
                currentRadius = coneBaseRadius * (1 - y_offset/coneHeight);
                
                for dx = -ceil(currentRadius):ceil(currentRadius)
                    for dz = -ceil(currentRadius):ceil(currentRadius)
                        if dx^2 + dz^2 <= currentRadius^2
                            voxelX = x + dx;
                            voxelY = y;
                            voxelZ = z + dz;
                            
                            if voxelX >= 1 && voxelX <= ct.cubeDim(2) && ...
                               voxelY >= 1 && voxelY <= ct.cubeDim(1) && ...
                               voxelZ >= 1 && voxelZ <= ct.cubeDim(3)
                                modulatorCones(voxelY, voxelX, voxelZ) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Combine base and cones
    finalMask = (modulatorMask | modulatorCones);
    
    %% Finalize modulator
    cst{ixModulator,4}{1} = find(finalMask);
    ct.cubeHU{1}(cst{ixModulator,4}{1}) = HU;
end