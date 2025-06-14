function [ct, cst] = makeCST_water(ct, RadiusPTV, RadiusOAR, PTV_center, HU_tissue, HU_PTV, HU_OAR, y_offset_water,PTV_Dose)
disp("makeCST function called!")
% PTV_center is optional [y x z] position - if not provided, uses center of CT

if ~exist('HU_tissue','var') || isempty(HU_tissue)
    HU_tissue = -1000;
end
if ~exist('HU_PTV','var') || isempty(HU_PTV)
    HU_PTV = 0;
end
if ~exist('HU_OAR','var') || isempty(HU_OAR)
    HU_OAR = -1000;
end
if ~exist('PTV_center','var') || isempty(PTV_center)
    PTV_center = round([ct.cubeDim./2]); % Default to center if not provided
end
% RadiusPTV = 40;
% RadiusOAR = 15;
if ~exist('y_offset_water','var') || isempty(y_offset_water)
    y_offset_water = 0;
end
if ~exist('PTV_Dose','var') || isempty(PTV_Dose)
    PTV_Dose = 11.6;
end

% Ensure structures fit in the CT volume
if 2*RadiusPTV > min(ct.cubeDim)
    error('PTV is too large for the CT volume.');
end
if 2*RadiusOAR > min(ct.cubeDim)
    error('OAR is too large for the CT volume.');
end

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
%In my case I only have left and right
ixTissue = 1;
ixPTV = 2;
ixOAR=3;
ixWater = 4;

cst{ixTissue,1} = 0;
cst{ixTissue,2} = 'tissue';
cst{ixTissue,3} = 'IGNORED';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'PTV';
cst{ixPTV,3} = 'TARGET';

cst{ixOAR,1} = 2;
cst{ixOAR,2} = 'OAR';
cst{ixOAR,3} = 'OAR';

cst{ixWater,1} = 3;
cst{ixWater,2} = 'water';
cst{ixWater,3} = 'EXTERNAL';

cst{ixTissue,5}.TissueClass  = 1;
cst{ixTissue,5}.alphaX       = 0.1000;
cst{ixTissue,5}.betaX        = 0.0500;
cst{ixTissue,5}.Priority     = 2;
cst{ixTissue,5}.Visible      = 1;
cst{ixTissue,5}.visibleColor = [0 0 0];
cst{ixTissue,5}.sfudOptimization=0;

cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,5}.visibleColor = [1 0 1];
cst{ixPTV,5}.sfudOptimization=0;
cst{ixPTV, 6}{1, 1}.parameters{1, 1}  = PTV_Dose;

cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 1;
cst{ixOAR,5}.Visible     = 1;
cst{ixOAR,5}.visibleColor = [1 1 0];
cst{ixOAR,5}.sfudOptimization=0;

cst{ixWater,5}.TissueClass = 1;
cst{ixWater,5}.alphaX = 0.1000;
cst{ixWater,5}.betaX  = 0.0500;
cst{ixWater,5}.Priority = 2;
cst{ixWater,5}.Visible = 1;
cst{ixWater,5}.visibleColor = [0 0 1];
cst{ixWater,5}.sfudOptimization=0;

% define objective as struct for compatibility with GNU Octave I/O
cst{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(5*10^3,11.6));

% cst{ixOAR,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(5*10^5,8.05));


%% For definition of either a cubic or spheric phantom
cubeHelper = zeros(ct.cubeDim);

% radiusTissue = xDim;
radiusTissue = ct.cubeDim(1,2);

%% Creation of tissue 
for x = 1:ct.cubeDim(1,2)
    for y = 1:ct.cubeDim(1,1)
        for z = 1:ct.cubeDim(1,3)
            currPost = [y x z] - round([ct.cubeDim./2]);
            if  sqrt(sum(currPost.^2)) <= radiusTissue
                cubeHelper(y,x,z) = 1;
            end
        end
    end
end

% extract the voxel indices and save it in the cst
cst{ixTissue,4}{1} = find(cubeHelper);

%% Creation of PTV at specified center
center = round([ct.cubeDim./2]);
cubeHelper1 = zeros(ct.cubeDim);
for x = 1:ct.cubeDim(1,2)
    for y = 1:ct.cubeDim(1,1)
        for z = 1:ct.cubeDim(1,3)
            currPost = [y x z] - center + PTV_center;
            if sqrt(sum(currPost.^2)) <= RadiusPTV
                cubeHelper1(y,x,z) = 1;
            end
        end
    end
end

%% Create OAR within bounds
min_y_pos = RadiusOAR + 5;  % Keep at least 5mm from edge
OAR_y_pos = center(1) - RadiusPTV - RadiusOAR - 20;
OAR_y_pos = max(OAR_y_pos, min_y_pos);  % Don't go below minimum

OAR_center = [OAR_y_pos, center(2), center(3)];

cubeHelper2 = zeros(ct.cubeDim);
for x = 1:ct.cubeDim(2)
    for y = 1:ct.cubeDim(1)
        for z = 1:ct.cubeDim(3)
            currPost = [y x z] - OAR_center;
            dist = sqrt(sum(currPost.^2));
            if dist <= RadiusOAR && all(currPost >= 0)  % Additional check
                % cubeHelper2(y,x,z) = 1;
            end
        end
    end
end
%% Put phantom in water pool
waterHelper = zeros(ct.cubeDim);
start_water = round(ct.cubeDim(1)/2) + y_offset_water;
% Fill ALL voxels below water level (not just y=water_level)
for y = 2:ct.cubeDim(1) - 1
    for x = 2:ct.cubeDim(2) - 1
        for z = 2:ct.cubeDim(3) -1
            if y >= start_water  % All voxels BELOW this y-level
                waterHelper(y,x,z) = 1;
            end
        end
    end
end

%% extract the voxel indices and save it in the cst
cst{ixPTV,4}{1} = find(cubeHelper1);
cst{ixOAR, 4}{1,1}=find(cubeHelper2);
cst{ixWater,4}{1} = find(waterHelper);

Ind_Total_PTV = cst{ixPTV,4}{1};
Ind_OAR = cst{ixOAR,4}{1};

ct.cubeHU{1}(cst{ixTissue,4}{1}) = HU_tissue;
ct.cubeHU{1}(Ind_Total_PTV) = HU_PTV;
ct.cubeHU{1}(Ind_OAR) = HU_OAR;
ct.cubeHU{1}(cst{ixWater,4}{1}) = 0; 

end
