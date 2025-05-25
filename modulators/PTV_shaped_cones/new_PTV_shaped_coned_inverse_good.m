function [ct, cst] = new_PTV_shaped_coned_inverse_good(ct, cst, modulatorDepth, yLocation, coneBaseRadius, maxConeHeight, coneSpacingX, coneSpacingZ, HU, cutout_start, width_margingsYXZ, structureColor, varargin)
    % Parameters with defaults
    if ~exist('modulatorDepth','var') || isempty(modulatorDepth), modulatorDepth = 20; end
    if ~exist('HU','var') || isempty(HU), HU = 132; end
    if ~exist('coneBaseRadius','var') || isempty(coneBaseRadius), coneBaseRadius = 3; end
    if ~exist('maxConeHeight','var') || isempty(maxConeHeight), maxConeHeight = 10; end
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
    
    %% Calculate PTV height profile (for cone length adjustment)
    % Create a map of PTV depth at each (x,z) position
    PTV_depth_map = zeros(ct.cubeDim(2), ct.cubeDim(3));
    for i = 1:length(ptv_x)
        x = ptv_x(i);
        z = ptv_z(i);
        PTV_depth_map(x,z) = max(PTV_depth_map(x,z), ptv_y(i));
    end
    
    % Normalize depth map (0 = no PTV, 1 = deepest PTV)
    max_depth = max(PTV_depth_map(:));
    min_depth = min(PTV_depth_map(PTV_depth_map>0));
    PTV_depth_map_norm = (PTV_depth_map - min_depth) / (max_depth - min_depth);
    PTV_depth_map_norm(isnan(PTV_depth_map_norm)) = 0;
    
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
    
    %% Create PTV cutout from BOTTOM
    cutout_start_vox = round(cutout_start/ct.resolution.y);
    bottom_of_modulator = Modulator_center(1) + halfSize(1);
    new_y = bottom_of_modulator - cutout_start_vox - (max(ptv_y) - ptv_y) - 1; % -1, in order to not cutout the lowest pixel
    
    valid = new_y >= 1 & new_y <= ct.cubeDim(1) & ...
            new_y >= (Modulator_center(1) - halfSize(1));
    moved_PTV = sub2ind(ct.cubeDim, new_y(valid), ptv_x(valid), ptv_z(valid));
    
    modulatorMask(moved_PTV) = 0;  % Remove PTV from base
    
    %% Create PTV projection map for cone height adjustment
    % Project PTV shape onto XZ plane to determine cone lengths
    PTV_projection = zeros(ct.cubeDim(2), ct.cubeDim(3));
    for i = 1:length(ptv_x)
        PTV_projection(ptv_x(i), ptv_z(i)) = max(PTV_projection(ptv_x(i), ptv_z(i)), ptv_y(i));
    end
    
    %% Find actual TOP surface after PTV removal
    top_surface = zeros(ct.cubeDim(2), ct.cubeDim(3));
    for x = 1:ct.cubeDim(2)
        for z = 1:ct.cubeDim(3)
            y_values = find(modulatorMask(:,x,z));
            if ~isempty(y_values)
                top_surface(x,z) = min(y_values); % Lowest y = top
            end
        end
    end
    
        %% Calculate base thickness map
    % First create a map of how much base material exists at each (x,z) position
    base_thickness = zeros(ct.cubeDim(2), ct.cubeDim(3));
    for x = 1:ct.cubeDim(2)
        for z = 1:ct.cubeDim(3)
            y_values = find(modulatorMask(:,x,z));
            if ~isempty(y_values)
                base_thickness(x,z) = length(y_values); % Thickness in voxels
            end
        end
    end
    
    % Normalize thickness (0 = no base, 1 = full thickness)
    min_thickness = min(base_thickness(base_thickness>0));
    max_thickness = max(base_thickness(:));
    thickness_factor = 1 - (base_thickness - min_thickness)/(max_thickness - min_thickness);
    thickness_factor(isnan(thickness_factor)) = 0;

%% Create cones with inverse length relationship to base thickness
% x_positions = linspace(Modulator_center(2) - halfSize(2), Modulator_center(2) + halfSize(2), round(2*halfSize(2)/coneSpacingX)+1);
% z_positions = linspace(Modulator_center(3) - halfSize(3), Modulator_center(3) + halfSize(3), round(2*halfSize(3)/coneSpacingZ)+1);

x_start = Modulator_center(2) - halfSize(2);
x_end = Modulator_center(2) + halfSize(2);
z_start = Modulator_center(3) - halfSize(3);
z_end = Modulator_center(3) + halfSize(3);

x_positions = x_start : coneSpacingX : x_end;
z_positions = z_start : coneSpacingZ : z_end;

    modulatorCones = zeros(ct.cubeDim);
    
    for x_idx = 1:length(x_positions)
        for z_idx = 1:length(z_positions)
            x = round(x_positions(x_idx));
            z = round(z_positions(z_idx));

            % Clamp to valid range
            x = max(1, min(x, size(top_surface, 1)));
            z = max(1, min(z, size(top_surface, 2)));

            y_base = top_surface(x, z);
            if y_base <= 0
                continue;
            end

            % Calculate cone length - longer where base is thinner
            if x >= 1 && x <= size(thickness_factor,1) && z >= 1 && z <= size(thickness_factor,2)
                cone_length = maxConeHeight * thickness_factor(x,z);
            else
                cone_length = maxConeHeight * 0.5;
            end
            
            % Ensure minimum length
            min_length = maxConeHeight * 0.3;
            cone_length = max(cone_length, min_length);
            
            % Create the cone
            for y_offset = 0:cone_length
                y = y_base - y_offset;
                if y < 1, break; end
                
                currentRadius = coneBaseRadius * (1 - y_offset/cone_length);
                
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