function emba_ehgf_binary_pu_tbt_plotTraj(r, plotsd)
%
% !!! CAUTION !!!: Adjusted from tapas_ehgf_binary_config.m for the EMBA 
% project by Anna Yurova, LMU Munich. 
%
% Here is the documentation from TAPAS:
%
% Plots the estimated or generated trajectories for the binary HGF perceptual model
% Usage example:  est = tapas_fitModel(responses, inputs); tapas_hgf_binary_plotTraj(est);
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2020 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Set up display
scrsz = get(0,'screenSize');
outerpos = [0.2*scrsz(3),0.2*scrsz(4),0.8*scrsz(3),0.8*scrsz(4)];
figure(...
    'OuterPosition', outerpos,...
    'Name', 'HGF trajectories');

% Time axis
t = ones(1,size(r.u,1));

ts = cumsum(t);
ts = [0, ts];

% Number of levels
lvl = r.c_prc.n_levels;

%% Upper levels
for j = 1:lvl-1

    % Subplots
    subplot(lvl,1,j);

    if plotsd == true
        % Prior +- sigma prior as the starting point
        upperprior = r.p_prc.mu_0(lvl-j+1) +sqrt(r.p_prc.sa_0(lvl-j+1));
        lowerprior = r.p_prc.mu_0(lvl-j+1) -sqrt(r.p_prc.sa_0(lvl-j+1));
        plot(0, upperprior, 'ob', 'LineWidth', 1);
        hold all;
        plot(0, lowerprior, 'ob', 'LineWidth', 1);
        
        % Then, move on to filling the trajectory of the mu +- sigma
        upper = [upperprior; r.traj.mu(:,lvl-j+1)+sqrt(r.traj.sa(:,lvl-j+1))];
        lower = [lowerprior; r.traj.mu(:,lvl-j+1)-sqrt(r.traj.sa(:,lvl-j+1))];
        fill([ts, fliplr(ts)], [(upper)', fliplr((lower)')], ...
             'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
    end

    % Plot the values of mu
    plot(ts, [r.p_prc.mu_0(lvl-j+1); r.traj.mu(:,lvl-j+1)], 'b', 'LineWidth', 2);
    hold all;
    % Plot the prior at 0
    plot(0, r.p_prc.mu_0(lvl-j+1), 'ob', 'LineWidth', 2);

    % Add some limits and title
    xlim([0 ts(end)]);
    title(['Posterior expectation of $x_' num2str(lvl-j+1) '$'], 'FontWeight', 'bold', 'interpreter','latex');
    ylabel(['$\mu_', num2str(lvl-j+1) '$'], 'interpreter','latex');

end

%% Input level
subplot(lvl,1,lvl);

% Plotting s(mu_2) instead of mu_1. 
plot(ts, [tapas_sgm(r.p_prc.mu_0(2), 1); tapas_sgm(r.traj.mu(:,2), 1)], 'r', 'LineWidth', 2);
hold all;

% Prior at 0
plot(0, tapas_sgm(r.p_prc.mu_0(2), 1), 'or', 'LineWidth', 2); 

% Add the actual input to the graph
plot(ts(2:end), r.u(:,1), '.', 'Color', [0 0.6 0]); % inputs

% Implied learning rate 
plot(ts(2:end), r.traj.wt(:,1), 'k') 

% add irregular responses
if ~isempty(find(strcmp(fieldnames(r),'irr'), 1))
    plot(ts(r.irr),  1.08.*ones([1 length(r.irr)]), 'x', 'Color', [1 0.7 0], 'Markersize', 4, 'LineWidth', 2); % irregular responses
    plot(ts(r.irr), -0.08.*ones([1 length(r.irr)]), 'x', 'Color', [1 0.7 0], 'Markersize', 4, 'LineWidth', 2); % irregular responses
end

% % Plot the normalised responses
% y = (r.y - min(r.y)) / ( max(r.y) - min(r.y) );
% plot(ts(2:end), y, 'Color', [0 0 0.7]);

% Add a title and label
title(['Input $u$ (green), PE weight (black) and posterior expectation of input $s(\mu_2)$ ', ...
       '(red) for $\omega_2$ = ', num2str(r.p_prc.om(2)), ' and $\omega_3$ = ' num2str(r.p_prc.om(3))], ...
  'FontWeight', 'bold', 'interpreter', 'latex');
ylabel('$u$, $s(\mu_2)$, learning rate', 'interpreter','latex');
axis([0 ts(end) -0.15 1.15]);

plot(ts(2:end), 0.5, 'k');
xlabel('Trial number', 'interpreter','latex');
hold off;

end