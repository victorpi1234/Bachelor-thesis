function [ct, cst] = working_prototype(ct, cst, modulatorDepth, yLocation, coneBaseRadius, coneLengthFactor, coneSpacingX, coneSpacingZ, HU, cutout_start, width_margingsYXZ, structureColor, varargin)
    % Parameters with defaults
    if ~exist('modulatorDepth','var') || isempty(modulatorDepth), modulatorDepth = 20; end
    if ~exist('HU','var') || isempty(HU), HU = 132; end
    if ~exist('coneBaseRadius','var') || isempty(coneBaseRadius), coneBaseRadius = 3; end
    if ~exist('coneLengthFactor','var') || isempty(coneLengthFactor), coneLengthFactor = 0.8; end % 80% of PTV depth
    if ~exist('coneSpacingX','var') || isempty(coneSpacingX), coneSpacingX = 2*coneBaseRadius; end
    if ~exist('coneSpacingZ','var') || isempty(coneSpacingZ), coneSpacingZ = 2*coneBaseRadius; end
    if ~exist('cutout_start','var') || isempty(cutout_start), cutout_start = 0; end
    if ~exist('width_margingsYXZ','var') || isempty(width_margingsYXZ), width_margingsYXZ = [5 5 5]; end
    if ~exist('structureColor','var') || isempty(structureColor), structureColor = [0 1 1]; end

  %% Find PTV and get voxels
    ixPTV = find(strcmp(cst(:,2), 'PTV'));
    if isempty(ixPTV), error('PTV structure not found in CST'); end
    PTV_voxels = cst{ixPTV,4}{1};
    [ptv_y, ptv_x, ptv_z] = ind2sub(ct.cubeDim, PTV_voxels);
    
    %% Calculate PTV center (median point)
    PTV_center = [median(ptv_y), median(ptv_x), median(ptv_z)];
    
    %% Create PTV height map for the "front" half (Y > center)
    front_half = ptv_y > PTV_center(1);
    front_y = ptv_y(front_half);
    front_x = ptv_x(front_half);
    front_z = ptv_z(front_half);
    
    % Create projection map of max Y values (height) for front half
    front_height_map = zeros(ct.cubeDim(2), ct.cubeDim(3));
    for i = 1:length(front_x)
        x = front_x(i);
        z = front_z(i);
        front_height_map(x,z) = max(front_height_map(x,z), front_y(i));
    end
    
    % Normalize height map (0-1 range)
    min_height = min(front_height_map(front_height_map>0));
    max_height = max(front_height_map(:));
    if max_height == min_height
        normalized_height = ones(size(front_height_map));
    else
        normalized_height = (front_height_map - min_height)/(max_height - min_height);
    end
    normalized_height(isnan(normalized_height)) = 0;

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
    new_y = bottom_of_modulator - cutout_start_vox - (max(ptv_y) - ptv_y);
    
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

    %% Create cones with lengths mirroring the front half PTV shape
    x_positions = Modulator_center(2)-halfSize(2):coneSpacingX:Modulator_center(2)+halfSize(2);
    z_positions = Modulator_center(3)-halfSize(3):coneSpacingZ:Modulator_center(3)+halfSize(3);
    modulatorCones = zeros(ct.cubeDim);
    
    for x_idx = 1:length(x_positions)
        for z_idx = 1:length(z_positions)
            x = round(x_positions(x_idx));
            z = round(z_positions(z_idx));
            
            % Ensure indices are within bounds
            if x < 1 || x > size(top_surface,1) || z < 1 || z > size(top_surface,2)
                continue;
            end

            y_base = top_surface(x,z);
            if y_base == 0, continue; end
            
            % Get normalized height from front half PTV
            if x >= 1 && x <= size(normalized_height,1) && z >= 1 && z <= size(normalized_height,2)
                height_factor = normalized_height(x,z);
            else
                height_factor = 0;
            end
            
            % Calculate cone length (mirroring front shape)
            max_cone_length = modulatorDepth * coneLengthFactor;
            cone_length = height_factor * max_cone_length;
            
            % Minimum length to ensure some modulation
            cone_length = max(cone_length, coneBaseRadius*1.5);
            
            % Create the cone
            for y_offset = 0:cone_length
                y = y_base - y_offset;
                if y < 1 || y > ct.cubeDim(1), break; end
                
                currentRadius = coneBaseRadius * (1 - y_offset/cone_length);
                radius_ceil = ceil(currentRadius);
                
                for dx = -radius_ceil:radius_ceil
                    for dz = -radius_ceil:radius_ceil
                        if dx^2 + dz^2 <= currentRadius^2
                            voxelX = x + dx;
                            voxelY = y;
                            voxelZ = z + dz;
                            
                            % Final boundary check
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