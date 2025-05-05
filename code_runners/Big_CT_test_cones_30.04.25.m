%%
clear all
close all
clc
%% set matRad configuration
addpath(genpath('C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\matRad-dev'));
matRad_rc;
%%
addpath(genpath('C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\MymatRad\Project'));
%%
ct_no_mod = makeCT(200,200,100, 1);

% Check for resume state
% resume = false;
% if isfile('resume_state.mat')
%     load('resume_state.mat');  % Loads `currentState`
%     resume = true;
% end

%% stf


clearvars ct cst
[ct,cst] = new_PTV_shaped_coned_inverse_good(ct_no_mod,cst_no_mod,40, 60,3,6,8,6,132);
pln = specific_pln_maker(ct, cst, 1);

stf = calc_STF_Single_Energy_Layer_Beam(ct, cst, pln);


% Calculation without optimization
w = ones(stf.totalNumOfBixels,1);
resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);
resultGUI_no_mod = matRad_calcDoseForward(ct_no_mod,cst_no_mod, stf, pln, w);



ixProfileX = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
ixProfileX = ixProfileX(1)-8:ixProfileX(1)+8;
ixProfileY = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
for ixProfileY = ixProfileY(2):ixProfileY(2)+8
    slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
    slice = slice(3);
    profile = resultGUI.physicalDose(:,ixProfileY,slice);
    % Store both profiles for comparison
    
    figure,plot(profile,'LineWidth',2),grid on,hold on,
    plot(profile,'LineWidth',2),legend({'original profile'}),
    xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
end



%%
%matRadGUI

% To fix: boxSize_y doesn't work with odd numbers
% Make also, that the names of every file, meaning  title_text = sprintf('Cone-width=3cm Cone-spacing=5cm \n s-water-to-phantom = 35cm HU = 132');