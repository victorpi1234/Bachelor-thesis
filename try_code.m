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
% [ct_no_mod, cst_no_mod] = makePlanCTandCST(300,300,100,1);

% ct_no_mod = makeCT(650,650,350,10);
ct_no_mod = makeCT(300,300,100,10);
[ct_no_mod, cst_no_mod] = makeCST_water(ct_no_mod, 40, 20, [-100 0 0]);
% [ct_no_mod, cst_no_mod] = makeCST_water(ct_no_mod, 40, 20, [-100 0 0]);
% [ct_no_mod, cst_no_mod] = bu(ct_no_mod, 40, 20, [-10 0 0]);
% ct = makeCT(900,900,100,1);
% [ct_no_mod, cst_no_mod] = makeCST(ct_no_mod, 40, 20);

% ct = makeCT(650,650,250,1);
[ct_no_mod, cst_no_mod] = makeCST(ct_no_mod, 40, 20);




ct = makeCT(300,300,100,10);
[ct, cst] = makeCST_water(ct, 40, 20, [-100 0 0]);




% [ct, cst] = makeCST_water(ct, 40, 20, [-10 0 0]);

% % [ct_no_mod, cst_no_mod] = makeCST(ct_no_mod, 40, 20);
% [ct, cst] = makeCST_modified(ct, 40, 20, [-100 0 0]);
%%
clearvars ct cst
% [ct,cst] = new_PTV_shaped_coned_inverse_good(ct_no_mod,cst_no_mod,40, 200,4,80,6,6,132);
% [ct,cst] = working_PTV_shaped_cone_mod(ct_no_mod,cst_no_mod,40,200,2,4,[],[],1527, 0);

% [ct, cst] = PTV_shaped_coned_inverse(ct_no_mod,cst_no_mod,40, 100,4,30,6,6,132 );
% [ct, cst] = PTV_Shaped_Modulator_with_Cones(ct_no_mod,cst_no_mod,40, 100,4,30,6,6,132 );
% [ct, cst] = PTV_Shaped_Modulator(ct_no_mod,cst_no_mod,60, 100);

[ct,cst] = Box_modulator(ct_no_mod,cst_no_mod, [40 80 50], 100, [4 4],[],[],132);

% %%
% clearvars ct cst
% x=40;
% [ct, cst] = Cone_modulator_upside_down(ct_no_mod,cst_no_mod, 3,30, 110, x ,20 ,[20, x, 30]);
% [ct, cst] = Cone_modulator_upside_down_Copy(ct,cst, 3,30, 300, -x ,20 ,[20, x, 30]);
% [ct, cst] = Cone_modulator_upside_down_Copy(ct,cst, 3,30, 280, 0 ,20 ,[20, x, 30]);

% [ct, cst] = delete_Cone_mod_Copy(ct,cst, 3,30, 300, 20 , [20, 150, 30],[], -1000);
% [ct, cst] = Cone_modulator_upside_down(ct,cst, 3,30, 280, 20 , [20, 150, 30]);
%%
% [ct,cst] = Box_modulator(ct_no_mod,cst_no_mod,[20, 100, 30],100, [5 5], [100 100], [], 300 ,[0 0] );
% [ct, cst] = Cone_modulator(ct,cst, 3,30, 300, [20, 300, 30]);
%% create pln
pln = specific_pln_maker(ct, cst, 10);

%% stf
stf = calc_STF_Single_Energy_Layer_Beam(ct, cst, pln);

%% Calculation without optimization
w = ones(stf.totalNumOfBixels,1);
resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,w);
resultGUI_no_mod = matRad_calcDoseForward(ct_no_mod,cst_no_mod, stf, pln, w);
%%
% resultGUI = plan(ct, cst, pln);
% [resultGUI_no_modulator, resultGUI, pln] = plan(ct_no_mod, cst_no_mod, ct, cst, 1);
%%
title_text = sprintf('PTV-shaped-mod-with-fancy-cones \n PTV-Radius=100; HU = 1527; Cone Radius= 2');
title_text = sprintf('Comb Modulator, HU=1527(Aluminum), Tooth size = 4mm \n PTV-radius=40 mm');
data_analysis_Dose_comparer(ct, cst,ct_no_mod, cst_no_mod, resultGUI, resultGUI_no_mod, pln, title_text, [], 50);
%% Save png
folderpng = 'C:\Users\victo\OneDrive\Работен плот\dskalo1\Heidelberg Physik\bachelor_thesis\dkfz\MymatRad\Project\new_code\pics\pics_presentation\';
saveas(gcf, [folderpng,'Cone_modulator_Alum.png'])
%%
matRadGUI



