function [ct, cst] = delete_Cone_mod_Copy(ct, cst, coneBaseRadius, coneHeight, yLocation, baseThickness, boxSize, coneSpacingX_Z, HU, varargin)
    % Parameter handling
    if ~exist('coneSpacingX_Z','var') || isempty(coneSpacingX_Z)
        coneSpacingX_Z = [2*coneBaseRadius 2*coneBaseRadius];
    end
    if ~exist('boxSize','var') || isempty(boxSize)
        boxSize = [20, 130, 30];
    end
    if ~exist('HU','var') || isempty(HU)
        HU = -1000; % Default to air
    end

    coneSpacingX = coneSpacingX_Z(1); 
    coneSpacingZ = coneSpacingX_Z(2);

    %% Find modulator structure
    ixModulator = find(strcmp(cst(:,2), 'modulator');
    if isempty(ixModulator)
        error('Modulator structure not found in CST');
    end
    
    %% Create deletion mask
    center = round([ct.cubeDim./2]);
    Modulator_center = [(center(1) - yLocation), center(2), center(3)];
    halfSize = boxSize / 2;
    deleteMask = false(ct.cubeDim);
    
    %% Mark base voxels for deletion (middle portion only)
    baseTopY = Modulator_center(1) + baseThickness - 1;
    if baseTopY > ct.cubeDim(1)
        baseTopY = ct.cubeDim(1);
    end
    
    % Only delete the middle portion (x and z dimensions)
    xRange = max(1,Modulator_center(2)-halfSize(2)):min(ct.cubeDim(2),Modulator_center(2)+halfSize(2));
    zRange = max(1,Modulator_center(3)-halfSize(3)):min(ct.cubeDim(3),Modulator_center(3)+halfSize(3));
    
    for x = xRange
        for y = Modulator_center(1):baseTopY
            for z = zRange
                deleteMask(y,x,z) = true;
            end
        end
    end
    
    %% Mark cone voxels for deletion (middle portion only)
    xConeRange = max(1,Modulator_center(2)-halfSize(2)):coneSpacingX:min(ct.cubeDim(2),Modulator_center(2)+halfSize(2));
    zConeRange = max(1,Modulator_center(3)-halfSize(3)):coneSpacingZ:min(ct.cubeDim(3),Modulator_center(3)+halfSize(3));
    
    for x = xConeRange
        for z = zConeRange
            for y = Modulator_center(1)-1:-1:Modulator_center(1)-coneHeight
                currentRadius = coneBaseRadius * (1 - (Modulator_center(1) - y)/coneHeight);
                for dx = -ceil(currentRadius):ceil(currentRadius)
                    for dz = -ceil(currentRadius):ceil(currentRadius)
                        if dx^2 + dz^2 <= currentRadius^2
                            voxelX = x + dx;
                            voxelY = y;
                            voxelZ = z + dz;
                            if voxelX >= 1 && voxelX <= ct.cubeDim(2) && ...
                               voxelY >= 1 && voxelY <= ct.cubeDim(1) && ...
                               voxelZ >= 1 && voxelZ <= ct.cubeDim(3)
                                deleteMask(voxelY,voxelX,voxelZ) = true;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Apply deletions
    % Get current modulator voxels
    modulatorVoxels = cst{ixModulator,4}{1};
    
    % Convert logical mask to linear indices
    deleteIndices = find(deleteMask);
    
    % Find intersection between modulator voxels and deletion mask
    voxelsToRemove = intersect(modulatorVoxels, deleteIndices);
    
    % Remove these voxels from the structure
    voxelsToKeep = setdiff(modulatorVoxels, voxelsToRemove);
    cst{ixModulator,4}{1} = voxelsToKeep;
    
    % Set HU values for deleted voxels
    ct.cubeHU{1}(voxelsToRemove) = HU;
    
    % Remove structure if empty
    if isempty(voxelsToKeep)
        cst(ixModulator,:) = [];
        ct.structureColor(ixModulator,:) = [];
    end
end