function [] = plotModelIdentifiability(modSpace, saveDir)
% Plot model identifiability plot for both the LME winner method: 'model_ident_classification'
% and the BMC RFX model comparison
% Based on hessetal_spirl_analysis/job_runner_plot_results.m, lines
% 768-810
%
% INPUT
%   modSpace            struct              model space
%
%   saveDir             char array          base output directory
%-----------------------------------------------------------------------------
%
% Copyright (C) 2024 Anna Yurova, Irene Sophia Plank, LMU University Hospital;
%               based on the hessetal_spirl_analysis toolbox by Alex Hess (2024), TNU, ETHZ
%
%
% This file is part of the EMBA_HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. For further details, see the file LICENSE or <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------

rec = load(fullfile(saveDir, 'sim', 'model_identifiability'));

%% REC: plot model ident (Classification)
% title
t2 = sprintf('(chance threshold (0.90-CI) = %.2f)', rec.class.chancethr);
t3 = sprintf('Balanced Acc = %.2f', rec.class.balacc);
comp_tit = {'Classification'; t2; t3};
% axis labels
mod_names = {};
mod_nr = {};
for m = 1:size(modSpace,2)
    mod_names{m} = modSpace(m).name;
    mod_nr{m} = ['m', num2str(m)];
end

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
bwr = @(n)interp1([1 2], [0.9 0.9 0.9; 0.3 0.5 1], linspace(1, 2, n), 'linear');
imagesc(rec.class.percLMEwinner)
colormap(bwr(64));
colorbar;
set(gca, 'clim', [0 1])
set(gca, 'Fontsize', 14);
title(comp_tit)
xlabel('Recovered');
ylabel('Simulated');
ax = gca;
ax.XTick = [1:size(modSpace,2)];
ax.XTickLabel = mod_names;
ax.YTick = [1:size(modSpace,2)];
ax.YTickLabel = mod_names;

pos=get(gca,'position');
[rows,cols]=size(rec.class.percLMEwinner);
nums = flip(rec.class.percLMEwinner',2);
width=pos(3)/(cols);
height=pos(4)/(rows);
for i=1:cols
      for j=1:rows               
        annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(j-1),width,height], ...
        'string',num2str(nums(i,j), '%.2f'),'LineStyle','none','HorizontalAlignment','center',...
        'VerticalAlignment','middle');
      end
end

figdir = fullfile(saveDir, 'figures', 'model_identifiability');
if ~exist(figdir, 'dir')
   mkdir(figdir)
end
print([figdir filesep 'model_ident_lme'], '-dpng');
print([figdir filesep 'model_ident_lme'], '-dsvg');
close;

%% Plot BMC

% plot model comparison
% title
comp_tit = {'rfx bmc pxp'};

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
bwr = @(n)interp1([1 2], [0.9 0.9 0.9; 0.3 0.5 1], linspace(1, 2, n), 'linear');
imagesc(rec.bmc.rfx.pxp)
colormap(bwr(64));
colorbar;
set(gca, 'clim', [0 1])
set(gca, 'Fontsize', 14);
title(comp_tit)
xlabel('Recovered');
ylabel('Simulated');
ax = gca;
ax.XTick = [1:size(modSpace,2)];
ax.XTickLabel = mod_names;
ax.YTick = [1:size(modSpace,2)];
ax.YTickLabel = mod_names;

pos=get(gca,'position');
[rows,cols]=size(rec.bmc.rfx.pxp);
nums = flip(rec.bmc.rfx.pxp',2);
width=pos(3)/(cols);
height=pos(4)/(rows);
for i=1:cols
      for j=1:rows               
        annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(j-1),width,height], ...
        'string',num2str(nums(i,j), '%.2f'),'LineStyle','none','HorizontalAlignment','center',...
        'VerticalAlignment','middle');
      end
end

figdir = fullfile(saveDir, 'figures', 'model_identifiability',...
    ['model_ident_bmc']);
print(figdir, '-dpng');
print(figdir, '-dsvg');
close;

end
