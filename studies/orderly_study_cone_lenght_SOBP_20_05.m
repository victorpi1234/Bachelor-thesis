%%
clear all
close all
clc 
%% set matRad configuration
addpath(genpath('C:\Users\victo\bachelor_thesis\dkfz\matRad-dev'));
matRad_rc;  
%%
addpath(genpath('C:\Users\victo\bachelor_thesis\dkfz\MymatRad\Project'));

%% Make CT, CST
ct_no_mod = makeCT(800,300,300, 0.25);
[ct_no_mod, cst_no_mod] = makeCST_water(ct_no_mod, 120, 20, [-250 0 0], [],[],[],-10, 40);
%% create pln
pln = specific_pln_maker(ct, cst, 1); %CT and CST only needed for Isocenter
pln.propDoseCalc.engine =    'MCsquare'; %HongPB

%% Cone Modulator
% for ConeHeight = 0: 
clearvars ct cst % ct_no_mod cst_no_mod
[ct, cst] = Straight_Cones_with_box_base(ct_no_mod,cst_no_mod,0.25 , 200,5,40,12, 12,600,[],[1 1 1]);

%% stf
stf = calc_STF_Single_Energy_Layer_Beam(ct, cst, pln);
% stf = STF_Specific_Single_Energy_Layer(ct, cst, pln);
%% Calculation without optimization
w = 10^4*ones(stf.totalNumOfBixels,1);
resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);
% resultGUI_no_mod = matRad_calcDoseForward(ct_no_mod,cst_no_mod, stf, pln, w);

%% Data analysis 
title_text = sprintf(["Cone Modulator 1px = 0.25mm, tooth lenght 40px \n tooth width 5px; 200px from Isocenter"]);
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
data_analysis_one_Dose(ct, cst, resultGUI, pln, title_text, slice(3)-3, "on", [0 50],15);
%% Homogeniety index
[ct,cst] = HI_parallelepiped_substructure(ct,cst, resultGUI, 40);
%%
cst{2, 6}{1, 1}.parameters{1, 1}  = 35;
qi = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose);
%% Bragg peak N_of_pics -number of pictures pictures
ixProfileX = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
ixProfileX = ixProfileX(1)-8:ixProfileX(1)+1;
ixProfileY = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
N_of_pics = 0;
for ixProfileY = ixProfileY(2):ixProfileY(2) + N_of_pics
    slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
    slice = slice(3);
    %ixProfileY = 140;
    profile = resultGUI.physicalDose(624,:,slice);
    % Store both profiles for comparison
    
    figure,plot(profile,'LineWidth',2),grid on,hold on,
    plot(profile,'LineWidth',2),legend({'original profile'}),
    xlabel('mm'),ylabel('Gy(RBE)'),title('profile plot')
end

%% Data analysis dose comparer  ixProfileY
title_text1 = 'Just range compensator';
title_text2 = sprintf(['Cone-width=1.25mm Cone-spacing=2.75mm \n Max-cone-height=95mm HU = 600']);
% title_text = sprintf('Comb Modulator, HU=1527(Aluminum), Tooth size = 4mm \n PTV-radius=40 mm');
data_analysis_Dose_comparer(ct_comb, cst_comb,ct, cst, resultGUI_comb, resultGUI, pln, title_text1,title_text2, [], slice(3)-2, [],[0 60],[-50 30]);
%% Save png
folderpng = 'C:\Users\victo\bachelor_thesis\dkfz\DKFZ_BachelorThesis\Thesis\Pictures\05_22\';
% saveas(gcf, [folderpng,'Cone_modulator_Difference_Cones_Lenght_80_width_3_HU_600_Z_axis_plus_3.png'])
saveas(gcf, [folderpng,'Dose_Distr_x_Comb_modulator_lenght_40px.png'])%Dose_Distr_y_

%%
dvh = matRad_calcDVH(cst, resultGUI.physicalDose);
% figure;
% matRad_showDVH(dvh, cst, pln);
%%
show_DVH_custom(dvh, cst, pln,'MaxDose',25);
title('Dose-Volume Histogram (DVH)');
%%
matRadGUI
