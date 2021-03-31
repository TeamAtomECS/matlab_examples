%% Recoil Limit Sweep
%

%%
% Perform the simulations. 
detunings = -400:25:-50;
i_over_isat = 1.9;

results = struct('detuning', {}, 'IoverIsat', {}, 'T', {});

for detuning=detunings
    detuning
    simulate(detuning, i_over_isat);

    % Examine the atoms that are still in the simulation at the last frame -
    % these are the ones we will measure the temperature of.
    output_v = utils.read_output(fullfile('atomecs','vel.txt'));
    trapped = output_v(end).id;

    % Loop through and measure the V^2 of the atoms which remain trapped.
    trapped_velocities = arrayfun(@(frame) frame.vec(ismember(frame.id, trapped), :), output_v, 'UniformOutput', 0);
    vSq = cellfun(@(v) mean(sum(v.^2,2)), trapped_velocities);

    % Convert to uK temperature using equipartition theorem.
    amu = 1.66e-27;
    kB = 1.38e-23;
    T = 1e6 * (amu * 87 * vSq / kB / 3);

    dT = 2e-6;
    file_step = 50;
    clf; set(gcf, 'Color', 'w');
    plot(file_step*dT*(1:length(T)), T);
    xlabel('time (s)');
    ylabel('T (uK)');

    % Take a fraction of the trajectory to determine temperature.
    % * Don't take the end, because an atom may be falling from the trap.
    % * Don't take the start, there may be transients due to sag.
    interval_start = floor(0.5 * length(T));
    interval_end = floor(0.8 * length(T));
    meanT = mean(T(interval_start:interval_end));
    
    results(end+1) = struct('detuning', detuning, 'IoverIsat', i_over_isat, 'T', meanT);
end

%%
% Save results
save('recoil_sweep_S1.9.mat', 'results', 'detunings', 'i_over_isat');

%%
% Plot a graph of simulated recoil T versus detuning
% 
% source: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjvkouV5trvAhWg_7sIHTvvD4wQFjABegQIAxAD&url=https%3A%2F%2Fhal.archives-ouvertes.fr%2Fhal-01068704%2Fdocument&usg=AOvVaw1vUJ1oBvdf27RD8Ev37Yzo
kB = 1.38e-23;
hb = (6.63e-34 / (2*pi));
Gamma = 2 * pi * 7e3;
dopplerT = @(delta, IOverISat) (hb * Gamma) / 2 * (1+6*IOverISat+(2*delta./Gamma).^2)./(4*abs(delta)./Gamma) / kB * 1e6;
clf; set(gcf, 'Color', 'w');
plot([results.detuning], [results.T], '.', 'MarkerSize', 10); hold on;
delta = [results.detuning]*2*pi*1e3;
plot([results.detuning], dopplerT(delta, 60));
xlabel('detuning (kHz)', 'interpreter', 'latex');
ylabel('T ($\mu$K)', 'interpreter', 'latex');
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'latex');
legend('sim', 'Doppler limit');

%%
% Plot a graph of simulated recoil T versus detuning
% 
% source: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjvkouV5trvAhWg_7sIHTvvD4wQFjABegQIAxAD&url=https%3A%2F%2Fhal.archives-ouvertes.fr%2Fhal-01068704%2Fdocument&usg=AOvVaw1vUJ1oBvdf27RD8Ev37Yzo
kB = 1.38e-23;
hb = (6.63e-34 / (2*pi));
Gamma = 2 * pi * 7e3;
dopplerT = @(delta, IOverISat) (hb * Gamma) / 2 * (1+6*IOverISat+(2*delta./Gamma).^2)./(4*abs(delta)./Gamma) / kB * 1e6;
clf; set(gcf, 'Color', 'w');
load('recoil_sweep_S1.9.mat');
plot([results.detuning], [results.T], '-', 'MarkerSize', 10); hold on;
delta = [results.detuning]*2*pi*1e3;
%plot([results.detuning], dopplerT(delta, 60));
load('recoil_sweep_S60.mat');
plot([results.detuning], [results.T], '-', 'MarkerSize', 10); hold on;
xlabel('detuning (kHz)', 'interpreter', 'latex');
ylabel('T ($\mu$K)', 'interpreter', 'latex');
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'latex');
legend('S=1.9', 'S=60');

%%

output_v = utils.read_output(fullfile(wd,'vel.txt'));
velocities = {output_v(:).vec};
vSq = cellfun(@(v) mean(sum(v.^2,2)), velocities);

% convert to uK
amu = 1.66e-27;
kB = 1.38e-23;
T = (amu * 87 * vSq / kB / 3);
T = T * 1e6;
hb = (6.63e-34 / (2*pi));
Gamma = 2 * pi * 7e3;
% Gamma = 2 * pi * 6e6;
TD = hb * Gamma / (2 * kB) * 1e6;

clf;
plot(10*(1:length(T)), T)
fprintf('Mean T=%.2f uK\n', mean(T))
hold on
plot(xlim, [1 1 ] * TD, '--k')
set(gcf, 'Color', 'w');
xlabel('$t$ ($\mu$s)', 'interpreter', 'latex');
ylabel('T ($\mu$K)', 'interpreter', 'latex');
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'latex');
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'latex');
grid on
set(gca, 'GridLineStyle', ':');
ylim([ 0 max(ylim) ]);

% Render to file
set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 7.5 ]);

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
saveas(gcf, 'doppler.pdf')


%%
rate = utils.read_output(fullfile(wd, 'rate.txt'), 'Format', '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,');
rates = cat(1,rate(2:end).vec);
mean(rates(:,1))

%%
utils.animate('File', fullfile(wd, 'pos.txt'), 'AxisView', [ 0 0 ], 'SimulationRegion', [ -150 150; -150 150; -150 150 ]*1e-6);

%%
utils.animate('File', fullfile(wd, 'pos.txt'), 'SimulationRegion', 1e-3*[ -1.0 1.0; -1.0 1.0; -1.0 1.0 ]);


%%
function simulate(detuning, IOverIsat)

% write config file
s = struct('detuning', detuning, ...
    'number_of_steps', 50000, ...
    'i_over_isat', IOverIsat, ...
    'quad_grad', 4);
jsonstr = jsonencode(s);
fH = fopen(fullfile('atomecs', 'recoil.json'),'w');
fprintf(fH, '%s', jsonstr);
fclose(fH);

wd = run_example('recoil_limit');

end