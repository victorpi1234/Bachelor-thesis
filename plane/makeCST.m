function [ct, cst] = makeCST(ct, RadiusPTV, RadiusOAR, HU_tissue, HU_PTV, HU_OAR) % Warning! Use a square CT!
disp("makeCST function called!")

if ~exist('HU_tissue','var') || isempty(HU_tissue)
    HU_tissue = 0;
end
if ~exist('HU_PTV','var') || isempty(HU_PTV)
    HU_PTV = 0;
end
if ~exist('HU_OAR','var') || isempty(HU_OAR)
    HU_OAR = 0;
end
% RadiusPTV = 40;
% RadiusOAR = 15;

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

% radiusTissue = xDim;
radiusTissue = ct.cubeDim(1,2);

% for x = 1:xDim
%     for y = 1:yDim
%         for z = 1:zDim
%             currPost = [y x z] - round([ct.cubeDim./2]);
%             if  sqrt(sum(currPost.^2)) <= radiusTissue
%                 cubeHelper(y,x,z) = 1;
%             end
%         end
%     end
% end

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

cubeHelper1 = zeros(ct.cubeDim);
cubeHelper2 = zeros(ct.cubeDim);


center = round([ct.cubeDim./2]);
% Creation of the central target
% for x = 1:xDim
%     for y = 1:yDim
%         for z = 1:zDim
%             currPost = [y x z] - center;
%             if  sqrt(sum(currPost.^2)) <= RadiusPTV
%                 cubeHelper1(y,x,z) = 1;
%             end
%         end
%     end
% end

for x = 1:ct.cubeDim(1,2)
    for y = 1:ct.cubeDim(1,1)
        for z = 1:ct.cubeDim(1,3)
            currPost = [y x z] - center;
            if  sqrt(sum(currPost.^2)) <= RadiusPTV
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



% OAR_center =[(center(1)-RadiusPTV-RadiusOAR-20.0) center(2) center(3)];
% Creation of the OAR
% for x = 1:xDim
%     for y = 1:yDim
%         for z = 1:zDim
%             currPost = [y x z] - OAR_center;
%             if  sqrt(sum(currPost.^2)) <= RadiusOAR
%                 cubeHelper2(y,x,z) = 1;
%             end
%         end
%     end
% end

% Creation of the OAR
% for x = 1:ct.cubeDim(1,2)
%     for y = 1:ct.cubeDim(1,1)
%         for z = 1:ct.cubeDim(1,3)
%             currPost = [y x z] - OAR_center;
%             if  sqrt(sum(currPost.^2)) <= RadiusOAR
%                 cubeHelper2(y,x,z) = 1;
%             end
%         end
%     end
% end

% extract the voxel indices and save it in the cst
cst{ixPTV,4}{1} = find(cubeHelper1);
cst{ixOAR, 4}{1,1}=find(cubeHelper2);

Ind_Total_PTV = cst{ixPTV,4}{1};
Ind_OAR = cst{ixOAR,4}{1};

ct.cubeHU{1}(cst{ixTissue,4}{1}) = HU_tissue;
ct.cubeHU{1}(Ind_Total_PTV) = HU_PTV;
ct.cubeHU{1}(Ind_OAR) = HU_OAR;
end