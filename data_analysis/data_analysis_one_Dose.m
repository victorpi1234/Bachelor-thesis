function [] = data_analysis_one_Dose(ct_mod,cst_mod, resultGUI, pln, modulator_type, slice, pics_on_or_off, Doserange, Fontsize)
% if nargin < 9 || isempty(rows)
%     rows = 1;
% end

if ~exist('slice','var') || isempty(slice)
    slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct_mod);
    slice = slice(3);
end
if ~exist('pics_on_or_off','var') || isempty(pics_on_or_off)
    pics_on_or_off = "off";
end
if ~exist('Doserange','var') || isempty(Doserange)
    Doserange = []; % Default: automatic color limits
end
if ~exist('Fontsize','var') || isempty(Fontsize)
    Fontsize = 20; % Default: automatic color limits
end
disp('Results function called!');
% %figure, matRad_plotSliceWrapper(gca,ct_mod_comb,cst_mod_comb,1,resultGUI.physicalDose - resultGUI_comb.physicalDose ,3,slice);
% %figure, matRad_plotSliceWrapper(gca,ct_no_mod,cst_no_mod,1,resultGUI_no_mod.physicalDose,3,slice);
% figure, matRad_plotSliceWrapper(gca,ct_mod,cst_mod,1,resultGUI.physicalDose,3,slice);
% clim([-5e-3 5e-3]); % Dose axis
% ylabel(colorbar, 'Dose (Gray)'); % Add a label to the colorbar
% title('No modulator') %With no modulator res=1
% %font size change
% fontsize = 20;
% set(gca, 'FontSize', fontsize); % Set font size for axes text (ticks, labels, etc.)
% set(get(gca, 'Title'), 'FontSize', fontsize); % Set font size for the title
% set(get(gca, 'XLabel'), 'FontSize', fontsize); % Set font size for the x-axis label
% set(get(gca, 'YLabel'), 'FontSize', fontsize); % Set font size for the y-axis label
% set(cb, 'FontSize', fontsize); % Set font size for the colorbar text
% set(get(cb, 'Label'), 'FontSize', fontsize); % Set font size for the colorbar label

%% Plot
figure
set(gcf, 'Visible', pics_on_or_off); 
% subplot(rows,1,1);
matRad_plotSliceWrapper(gca, ...
    ct_mod, cst_mod, 1, resultGUI.physicalDose, ...
    3, slice, [], [], [], ...     % thresh, alpha, contourColorMap
    [], Doserange, [], [], 'Dose [Gy]', false); 
%clim([0 5e-3]);
ylabel(colorbar, 'Dose [Gy]');


%% Apply font settings AFTER plotting
set(gca, 'FontSize', Fontsize);
set(get(gca, 'Title'), 'FontSize', Fontsize);
set(get(gca, 'XLabel'), 'FontSize', Fontsize);
set(get(gca, 'YLabel'), 'FontSize', Fontsize);

% Handle colorbar separately (must exist first)
cb = colorbar;
set(cb, 'FontSize', Fontsize);
set(get(cb, 'Label'), 'FontSize', Fontsize);
ylabel(cb, 'Dose [Gy]');



title([modulator_type]);
end