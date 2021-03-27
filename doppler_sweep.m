%% Doppler Sweep
% For the 3D MOT setup, measure the velocities and verify that the Doppler
% T limit is obeyed.
%%
% Perform the simulations

linewidth = 6; % MHz
detunings = -([0.2:0.05:0.45 0.5:0.1:2])*linewidth;
% detunings = -0.5*linewidth;
results = cell(1,length(detunings));

for i=1:length(detunings)
    detuning = detunings(i);
    
    % Write a configuration file for the simulation
    fid = fopen(fullfile('atomecs', 'doppler.json'), 'w');
    params = struct('detuning', detuning, 'number_of_steps', int32(3000));
    fprintf(fid, '%s', jsonencode(params));
    fclose(fid);
    
    wd = run_example('doppler_sweep');
    
    output_v = utils.read_output(fullfile(wd,'vel.txt'));
    velocities = {output_v(:).vec};
    vSq = cellfun(@(v) mean(sum(v.^2,2)), velocities);
    vx2 = cellfun(@(v) mean(v(:,1).^2), velocities);
    vy2 = cellfun(@(v) mean(v(:,2).^2), velocities);
    vz2 = cellfun(@(v) mean(v(:,3).^2), velocities);
    
    % convert to uK
    amu = 1.66e-27;
    kB = 1.38e-23;
    T = (amu * 87 * vSq / kB / 3);
    Tx = (amu * 87 * vx2 / kB) * 1e6;
    Ty = (amu * 87 * vy2 / kB) * 1e6;
    Tz = (amu * 87 * vz2 / kB) * 1e6;
    T = T * 1e6;
    results{i} = struct('T', T, 'Tx', Tx, 'Ty', Ty, 'Tz', Tz);
end

save('doppler_sweep.mat', 'results', 'detunings');

%% Plotting
% Plot the comparison graphs.
% 
% * The results vector holds all of the trajectories versus time.

clf;
cmap = hot(10);
getColor = @(delta) interp1(linspace(0,1,10), cmap, abs(delta)/max(abs(detunings)), 'linear');

% Doppler limit for Rb
gamma = 6e6 * 2 * pi;
hb = 6.63e-34 / (2*pi);
doppler = (hb * gamma / 2) / kB;

finalT = zeros(length(detunings),1);
for i=1:length(detunings)
    T = results{i}.T;
    finalT(i) = mean(T(floor(length(T)*[0.8 1.0])));
    plot(10*(1:length(T)), T, 'Color', getColor(detunings(i)));
    hold on;
%     plot(10*(1:length(T)), results{i}.Tx, '-', 'Color', getColor(detunings(i)));
%     plot(10*(1:length(T)), results{i}.Ty, ':', 'Color', getColor(detunings(i)));
%     plot(10*(1:length(T)), results{i}.Tz, '--', 'Color', getColor(detunings(i)));
end
plot(xlim, [1 1 ] * doppler * 1e6, '--k')
set(gcf, 'Color', 'w');
xlabel('$t$ ($\mu$s)', 'interpreter', 'latex');
ylabel('T ($\mu$K)', 'interpreter', 'latex');
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'latex');
grid on
set(gca, 'GridLineStyle', ':');
ylim([ 0 max(ylim) ]);

%%
% Plot a graph of Doppler T versus detuning
finalT = zeros(length(detunings),1);
for i=1:length(detunings)
    T = results{i}.T;
    finalT(i) = mean(T(floor(length(T)*[0.8 1.0])));
end
dopplerLim = @(delta) (1 / 4) * (1 + (2*delta/gamma).^2) / (2*abs(delta/gamma)) * (hb * gamma * 1e6 / kB);
plot(detunings/linewidth, finalT, '.k', 'MarkerSize', 10); hold on;
analytic = arrayfun(dopplerLim, 2*pi*detunings*1e6);
plot(detunings/linewidth, analytic, 'k-'); hold off;

xlabel('detuning $\delta/\Gamma$', 'interpreter', 'latex');
ylabel('T ($\mu$K)', 'interpreter', 'latex');
set(gcf, 'Color', 'w');
grid on;
set(gca, 'GridLineStyle', ':');
%%
% Produce a two panel plot with the data

clf;
set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 7.5 ]);

shown = [ 7 12 17 22 ];
cmap = parula(10);
getColor = @(delta) interp1(linspace(0,1,10), cmap, abs(delta)/16, 'linear');

% Doppler limit for Rb
gamma = 6e6 * 2 * pi;
hb = 6.63e-34 / (2*pi);
doppler = dopplerLim(0.5*gamma);

p_detunings = detunings(shown);
finalT = zeros(length(p_detunings),1);
for i=1:length(p_detunings)
    T = results{shown(i)}.T;
    hold on;
    plot(10*(1:length(T)), T, '-', 'Color', getColor(p_detunings(i)));
end
plot(xlim, [1 1 ] * doppler, '--k')
set(gcf, 'Color', 'w');
xlabel('$t$ ($\mu$s)', 'interpreter', 'latex');
ylabel('T ($\mu$K)', 'interpreter', 'latex');
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'latex');
grid on
set(gca, 'GridLineStyle', ':');
ylim([ 120 600 ]);
box on;

% Plot doppler limit panel
axes('Units', 'centimeters', 'Position', [ 5.02 4.3 3 2.5 ]);
box on;
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'latex');

for i=1:length(detunings)
    T = results{i}.T;
    finalT = mean(T(floor(length(T)*[0.8 1.0])));
    plot(detunings(i)/linewidth, finalT, '.', 'MarkerSize', 10, 'Color', getColor(detunings(i))); hold on;
end
dopplerLim = @(delta) (1 / 4) * (1 + (2*delta/gamma).^2) / (2*abs(delta/gamma)) * (hb * gamma * 1e6 / kB);
analytic = arrayfun(dopplerLim, 2*pi*detunings*1e6);
plot(detunings/linewidth, analytic, 'k-'); hold off;
xlabel('$\delta/\Gamma$', 'Interpreter', 'latex');
ylabel('$T_\textrm{D}$ ($\mu$K)', 'Interpreter', 'latex');
grid on
set(gca, 'GridLineStyle', ':');

set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
w = pos(3);
h = pos(4);
p = 0.01;
set(gcf,...
    'PaperUnits','centimeters',...
    'PaperPosition',[p*w p*h w h],...
    'PaperSize',[w*(1+2*p) h*(1+2*p)]);
set(gcf, 'Renderer', 'painters')
saveas(gcf, 'doppler_sweep.pdf')