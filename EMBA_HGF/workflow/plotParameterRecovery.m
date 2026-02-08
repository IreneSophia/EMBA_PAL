function [] = plotParameterRecovery(modSpace, saveDir)
% Plot parameter recovery plots: 'param_rec_m'
% Based on hessetal_spirl_analysis/job_runner_plot_results.m, lines
% 704-761
%
% !! CAUTION !!: This is HARDCODED for the models used in the EMBA project. 
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

res = struct();
res.rec = load(fullfile(saveDir, 'sim', 'recovery_analysis'));

for m = 1:size(modSpace,2)
    
    if contains(modSpace(m).name, "HGF")
    
        names_response   = {'$\beta_0$', '$\beta_1$', '$\beta_2$', ...
            '$\beta_3$', '$\beta_4$', '$\beta_5$'};
        names_perception = {'$\omega_2$', '$\omega_3$', '$\alpha$'};

    elseif contains(modSpace(m).name, "RW")

        names_response   = {'log$\zeta$'};
        names_perception = {'logit\alpha'};

    end

    npars = length(modSpace(m).prc_idx) + length(modSpace(m).obs_idx);
    nx = ceil(sqrt(npars));
    ny = round(sqrt(npars));
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    for j = 1:npars
        subplot(nx, ny, j)
        if j > length(modSpace(m).prc_idx) %% all obs param
            k = j-length(modSpace(m).prc_idx);
            scatter(res.rec.param.obs(m).sim(:,k), res.rec.param.obs(m).est(:,k), 15, 'k', 'filled');
            pcc = res.rec.param.obs(m).pcc(k);
            lsline;
            if k == size(modSpace(m).obs_idx, 2)
                title('response model: log($\Sigma$)', 'Interpreter','Latex')
            else
                title(['response model: ' names_response{k}], 'Interpreter','Latex')
            end
            sig = res.rec.param.obs(m).pval(k) < 0.05;
        else
            scatter(res.rec.param.prc(m).sim(:,j), res.rec.param.prc(m).est(:,j), 15, 'k', 'filled');
            pcc = res.rec.param.prc(m).pcc(j);
            lsline;
            title(['perception model: ' names_perception{j}], 'Interpreter','Latex')
            sig = res.rec.param.prc(m).pval(j) < 0.05;
        end        
        xlabel('Simulated values');
        ylabel('Recovered values');
        if sig 
            str = sprintf('r = %1.2f*', pcc);
        else
            str = sprintf('r = %1.2f', pcc);
        end
        textXpos = min(get(gca, 'xlim')) + (max(get(gca, 'xlim')) - min(get(gca, 'xlim')))*0.05;
        textYpos = max(get(gca, 'ylim')) - (max(get(gca, 'ylim')) - min(get(gca, 'ylim')))*0.05;
        T = text(textXpos, textYpos, str);
        set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    end
    figdir = fullfile(saveDir, 'figures', 'parameter_recovery');
    if ~exist(figdir, 'dir')
       mkdir(figdir)
    end
    print([figdir filesep 'param_rec_' convertStringsToChars(modSpace(m).name)], '-dpng');
    print([figdir filesep 'param_rec_' convertStringsToChars(modSpace(m).name)], '-dsvg');
    close;

end


end
