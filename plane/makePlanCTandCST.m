function [ct, cst] = makePlanCTandCST(xDim, yDim, zDim, res)
disp('makePlanCTandCST function called!');
ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate

ct.resolution.x = res;
ct.resolution.y = res;
ct.resolution.z = res;
ct.numOfCtScen  = res;

% create an ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000; % assign HU of Air


%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target
%In my case I only have left and right
ixTissue = 1;
ixPTV = 2;
ixOAR=3;

cst{ixTissue,1} = 0;
cst{ixTissue,2} = 'tissue';
cst{ixTissue,3} = 'IGNORED';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'PTV';
cst{ixPTV,3} = 'TARGET';

cst{ixOAR,1} = 2;
cst{ixOAR,2} = 'OAR';
cst{ixOAR,3} = 'OAR';

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

cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 1;
cst{ixOAR,5}.Visible     = 1;
cst{ixOAR,5}.visibleColor = [1 1 0];
cst{ixOAR,5}.sfudOptimization=0;

% define objective as struct for compatibility with GNU Octave I/O
cst{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(5*10^3,11.6));

% cst{ixOAR,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(5*10^5,8.05));

%% For definition of either a cubic or spheric phantom
cubeHelper = zeros(ct.cubeDim);

radiusTissue = xDim;

for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPost = [y x z] - round([ct.cubeDim./2]);
            if  sqrt(sum(currPost.^2)) <= radiusTissue
                cubeHelper(y,x,z) = 1;
            end
        end
    end
end

% extract the voxel indices and save it in the cst
cst{ixTissue,4}{1} = find(cubeHelper);

cubeHelper1 = zeros(ct.cubeDim);
cubeHelper2 = zeros(ct.cubeDim);

RadiusPTV = 40;
RadiusOAR = 15;

center = round([ct.cubeDim./2]);
% Creation of the central target
for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPost = [y x z] - center;
            if  sqrt(sum(currPost.^2)) <= RadiusPTV
                cubeHelper1(y,x,z) = 1;
            end
        end
    end
end

OAR_center =[(center(1)-RadiusPTV-RadiusOAR-20.0) center(2) center(3)];
% Creation of the OAR
for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPost = [y x z] - OAR_center;
            if  sqrt(sum(currPost.^2)) <= RadiusOAR
                cubeHelper2(y,x,z) = 1;
            end
        end
    end
end

% extract the voxel indices and save it in the cst
cst{ixPTV,4}{1} = find(cubeHelper1);
cst{ixOAR, 4}{1,1}=find(cubeHelper2);

Ind_Total_PTV = cst{ixPTV,4}{1};
Ind_OAR = cst{ixOAR,4}{1};

ct.cubeHU{1}(cst{ixTissue,4}{1}) = 0;
ct.cubeHU{1}(Ind_Total_PTV) = 0;
ct.cubeHU{1}(Ind_OAR) = 0;

%%
% ct_no_mod = ct; % Create a copy of ct
% cst_no_mod = cst; % Create a copy of cst
% 
% ixModulator = 4;
% 
% % define general VOI properties
% cst{ixModulator,1} = 1;
% cst{ixModulator,2} = 'modulator';
% cst{ixModulator,3} = 'EXTERNAL';
% 
% % define optimization parameter for both VOIs
% cst{ixModulator,5}.TissueClass  = 1;
% cst{ixModulator,5}.alphaX       = 0.1000;
% cst{ixModulator,5}.betaX        = 0.0500;
% cst{ixModulator,5}.Priority     = 2;
% cst{ixModulator,5}.Visible      = 1;
% cst{ixModulator,5}.visibleColor = [0 1 1];
% cst{ixModulator,5}.sfudOptimization=0;
% ct.structureColor(ixModulator, :) = [0 1 1];
end