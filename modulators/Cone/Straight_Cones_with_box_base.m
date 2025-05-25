function [ct, cst] = Straight_Cones_with_box_base(ct, cst, modulatorDepth_mm, yLocation_px, coneBaseRadius_px, cone_length_px, coneSpacingX_px, coneSpacingZ_px, HU, cutout_start, width_margingsYXZ, structureColor, varargin)
    % Parameters with defaults
    if ~exist('modulatorDepth_mm','var') || isempty(modulatorDepth_mm), modulatorDepth_mm = 20; end
    if ~exist('yLocation_px','var') || isempty(yLocation_px), yLocation_px = 0; end
    if ~exist('HU','var') || isempty(HU), HU = 132; end
    if ~exist('coneBaseRadius_px','var') || isempty(coneBaseRadius_px), coneBaseRadius_px = 3; end
    if ~exist('cone_length_px','var') || isempty(cone_length_px), cone_length_px = 10; end
    if ~exist('coneSpacingX_px','var') || isempty(coneSpacingX_px), coneSpacingX_px = 2*coneBaseRadius_px + 1; end
    if ~exist('coneSpacingZ_px','var') || isempty(coneSpacingZ_px), coneSpacingZ_px = 2*coneBaseRadius_px + 1; end
    if ~exist('cutout_start','var') || isempty(cutout_start), cutout_start = 0; end
    if ~exist('width_margingsYXZ','var') || isempty(width_margingsYXZ), width_margingsYXZ = [5 5 5]; end
    if ~exist('structureColor','var') || isempty(structureColor), structureColor = [0 1 1]; end

    %% Find PTV to determine width
    ixPTV = find(strcmp(cst(:,2), 'PTV'));
    if isempty(ixPTV), error('PTV structure not found in CST'); end
    
    PTV_voxels = cst{ixPTV,4}{1};
    [ptv_y, ptv_x, ptv_z] = ind2sub(ct.cubeDim, PTV_voxels);
    
    %% Calculate PTV height profile (for cone length adjustment)
    % % Create a map of PTV depth at each (x,z) position
    % PTV_depth_map = zeros(ct.cubeDim(2), ct.cubeDim(3));
    % for i = 1:length(ptv_x)
    %     x = ptv_x(i);
    %     z = ptv_z(i);
    %     PTV_depth_map(x,z) = max(PTV_depth_map(x,z), ptv_y(i));
    % end
    % 
    % % Normalize depth map (0 = no PTV, 1 = deepest PTV)
    % max_depth = max(PTV_depth_map(:));
    % min_depth = min(PTV_depth_map(PTV_depth_map>0));
    % PTV_depth_map_norm = (PTV_depth_map - min_depth) / (max_depth - min_depth);
    % PTV_depth_map_norm(isnan(PTV_depth_map_norm)) = 0;
    % 
    %% Calculate dimensions with margin
    x_margin = round(width_margingsYXZ(2)/ct.resolution.x);
    z_margin = round(width_margingsYXZ(3)/ct.resolution.z); 
    x_size = (max(ptv_x) - min(ptv_x)) + 2*x_margin;
    z_size = (max(ptv_z) - min(ptv_z)) + 2*z_margin;
    y_size = round(modulatorDepth_mm/ct.resolution.y);
    
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
    Modulator_center = [center(1)-yLocation_px, center(2), center(3)];
    halfSize = [y_size, x_size, z_size]/2;
    
    % Create base mask
    [Y,X,Z] = ndgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));
    modulatorMask = (abs(X - Modulator_center(2)) <= halfSize(2)) & ...
                   (abs(Y - Modulator_center(1)) <= halfSize(1)) & ...
                   (abs(Z - Modulator_center(3)) <= halfSize(3));
    
    %% Create PTV cutout from BOTTOM
    % cutout_start_vox = round(cutout_start/ct.resolution.y);
    % bottom_of_modulator = Modulator_center(1) + halfSize(1);
    % new_y = bottom_of_modulator - cutout_start_vox - (max(ptv_y) - ptv_y) - 1; % -1, in order to not cutout the lowest pixel
    % 
    % valid = new_y >= 1 & new_y <= ct.cubeDim(1) & ...
    %         new_y >= (Modulator_center(1) - halfSize(1));
    % moved_PTV = sub2ind(ct.cubeDim, new_y(valid), ptv_x(valid), ptv_z(valid));
    % 
    % modulatorMask(moved_PTV) = 0;  % Remove PTV from base
    
    %% Create PTV projection map for cone height adjustment
    % Project PTV shape onto XZ plane to determine cone lengths
    % PTV_projection = zeros(ct.cubeDim(2), ct.cubeDim(3));
    % for i = 1:length(ptv_x)
    %     PTV_projection(ptv_x(i), ptv_z(i)) = max(PTV_projection(ptv_x(i), ptv_z(i)), ptv_y(i));
    % end
    
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
%% Create cones with inverse length relationship to base thickness
% x_positions = linspace(Modulator_center(2) - halfSize(2), Modulator_center(2) + halfSize(2), round(2*halfSize(2)/coneSpacingX)+1);
% z_positions = linspace(Modulator_center(3) - halfSize(3), Modulator_center(3) + halfSize(3), round(2*halfSize(3)/coneSpacingZ)+1);

x_start = Modulator_center(2) - halfSize(2);
x_end = Modulator_center(2) + halfSize(2);
z_start = Modulator_center(3) - halfSize(3);
z_end = Modulator_center(3) + halfSize(3);

x_positions = x_start : coneSpacingX_px : x_end;
z_positions = z_start : coneSpacingZ_px : z_end;

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

            % Set cone length to maxConeHeight_px for all cones
            cone_length = cone_length_px;
            
            % Create the cone
            for y_offset = 0:cone_length
                y = y_base - y_offset;
                if y < 1, break; end
                
                currentRadius = coneBaseRadius_px * (1 - y_offset/(cone_length));

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

    %% Check modulator lenghts
    % Initialize array to store distances
    cone_heights = zeros(size(top_surface));

    % Find maximum height of each cone
    for x = 1:size(modulatorCones,2)
        for z = 1:size(modulatorCones,3)
            if top_surface(x,z) > 0  % Only where cones exist
                % Find all y-positions with cone material
                y_positions = find(modulatorCones(:,x,z));

                if ~isempty(y_positions)
                    % Calculate distance (top_surface - highest cone point)
                    cone_top = min(y_positions);  % Highest point (smallest y)
                    cone_heights(x,z) = top_surface(x,z) - cone_top;
                end
            end
        end
    end
    if max(cone_heights(:)) ~= cone_length_px
        error('Cones failed to reach full length! Maximum height achieved: %d pixels (expected: %d)', ...
            max(cone_heights(:)), cone_length_px);
    end
    % % Display results
    % fprintf('Cone height analysis:\n');
    % fprintf('Maximum designed cone length: %d pixels\n', cone_length_px);
    % fprintf('Actual cone lengths:\n');
    % fprintf('  Min: %d pixels\n', min(cone_heights(cone_heights>0)));
    % fprintf('  Max: %d pixels\n', max(cone_heights(:)));
    % fprintf('  Mean: %.1f pixels\n', mean(cone_heights(cone_heights>0)));

    % % Visualize the height map
    % figure;
    % imagesc(cone_heights);
    % colorbar;
    % title('Distance from Top Surface to Cone Tips (pixels)');
    % xlabel('X Position');
    % ylabel('Z Position');
    % axis equal tight;
end
