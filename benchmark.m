%% Benchmark the simulation performance
% This can take a long time to run!

thread_numbers = 1:1:6;
repeats = 10;
thread_atom_numbers = 10.^[2:0.5:6];
steps = 5e2;

%%
% Run once to force compilation.
bench(10, 1, steps);

thread_results = {};
for atom_number=thread_atom_numbers
    atom_results = {};
    for thread_number=thread_numbers
        repeat_results = {};
        for repeat=1:repeats
            time = bench(thread_number, atom_number, steps);
            repeat_results{end+1} = struct('threads', thread_number, 'atoms', atom_number, 'time', time);
        end
        atom_results{end+1} = cat(3, repeat_results{:});
    end
    atom_results = cat(2,atom_results{:});
    thread_results{end+1} = atom_results;
end
thread_results = cat(1,thread_results{:});
save('bench.mat', 'thread_results', 'steps');

%%
% Perform averaging
averaged_results = struct('threads', {}, 'atoms', {}, 'time', {});
std_results = struct('threads', {}, 'atoms', {}, 'time', {});
for i=1:size(thread_results, 1)
    for j=1:size(thread_results, 2)
        sim = thread_results(i,j,1);
        sim.time = mean([thread_results(i,j,:).time]);
        averaged_results(i,j) = sim;
        sim.time = std([thread_results(i,j,:).time]) / sqrt(repeats);
        std_results(i,j) = sim;
    end
end

%%
% Plot a graph showing the results

set(gcf, 'Units', 'centimeters');
pos = get(gcf, 'Position');
set(gcf, 'Position', [ pos(1) pos(2) 9 11 ]);
clf;
axes('Units', 'centimeters', 'Position', [ 1.2 4.7 7.5 6 ]);

c1 = [ 0.2 0.85 1.0 ];
c0 = 0*[ 0.0118 0.0196 0.1176 ];
c2 = [0.0830    0.6576    0.8615];
get_color = @(n) interp1([1; 2; 8], [ c1; c2; c0 ], log10(n), 'linear', true);

% cm = parula(10);
% norm = @(x, range) (x - min(range)) ./ (max(range)- min(range));
% get_color = @(n) interp1((length(cm):-1:1)'/length(cm), cm, norm(log10(n), -0.5:6.8));

set(gcf, 'Color', 'w');
h = [];
for i=1:1:size(averaged_results, 1)
   h(i) = plot([averaged_results(i,:).threads], 1e6*[averaged_results(i,:).time]./([averaged_results(i,:).atoms].*steps), '.-', 'Color', get_color(averaged_results(i,1).atoms)); hold on;
end
xlabel('', 'interpreter', 'latex', 'FontSize', 11);
% tau = total wall time per atom, per thread
ylabel('$\tau$ ($\mu$s) ', 'Interpreter', 'latex', 'FontSize', 11);
grid on;
set(gca, 'GridLineStyle', ':');
xlim([min(thread_numbers) max(thread_numbers)]);
set(get(gca, 'XAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(gca, 'YAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
% set(gca, 'YScale', 'log');
set(gca, 'XTick', 1:6, 'XTickLabel', {});
% set(gca, 'YScale', 'log');
xlim([1 6]);
labels = arrayfun(@(x) [num2str(x) ' atoms'], thread_atom_numbers, 'UniformOutput', 0);
selected = [ 1 5 length(labels) ];
legend(h(selected),labels{selected}, 'Interpreter', 'Latex');

% fit Amdahl's law and show on the graph
ax2 = axes('Units', 'centimeters', 'Position', [ 1.2 1.2 7.5 3.3 ]);
f = @(A, p, x) A*(1-p+p./x);
ft = fittype(f);
fitResult = fit([averaged_results(end,:).threads]', [averaged_results(end,:).time]', ft);
plot([averaged_results(end,:).threads]', [averaged_results(end,:).time]'/steps, 'o', 'Color', get_color(1e6)); hold on;
%errorbar([averaged_results(end,:).threads]', [averaged_results(end,:).time]'/steps, [std_results(end,:).time]'/steps, 'o', 'Color', c1); hold on;
xs = [averaged_results(end,:).threads]';
xs = linspace(min(xs), max(xs), 1000);
plot(xs, fitResult(xs)/steps, 'k--'); hold off;
grid on;
set(ax2, 'GridLineStyle', ':');
set(gca, 'XTick', 1:6, 'YTick', [ 0.5, 1.0 ]);
set(get(ax2, 'XAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
set(get(ax2, 'YAxis'), 'TickLabelInterpreter', 'Latex', 'FontSize', 10);
xlabel('number of threads, $N_t$', 'interpreter', 'latex', 'FontSize', 11);
% tau = total wall time per atom, per thread
ylabel('$T_\textrm{sim} / N_a$ (s) ', 'Interpreter', 'latex', 'FontSize', 11);
ylim([ 0 1.5 ]);
annotation('textbox', 'String', '(a)', 'Interpreter', 'Latex', 'FontSize', 11, 'Units', 'centimeters', 'LineStyle', 'none', 'Position', [-0.2 9.0 2 2])
annotation('textbox', 'String', '(b)', 'Interpreter', 'Latex', 'FontSize', 11, 'Units', 'centimeters', 'LineStyle', 'none', 'Position', [-0.2 2.8 2 2])

% Render to file
pos = get(gcf, 'Position');
w = pos(3); 
h = pos(4);
p = 0.01;
set(gcf,...
  'PaperUnits','centimeters',...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);
set(gcf, 'Renderer', 'painters')
saveas(gcf, 'bench.pdf')

function loop_time = bench(thread_number, atom_numbers, steps)
   
config = struct('n_threads', int32(thread_number), 'n_steps',  int32(steps), 'n_atoms',  int32(atom_numbers));
fH = fopen('atomecs/benchmark.json', 'w');
oc = onCleanup(@() fclose(fH));
fprintf(fH, '%s', jsonencode(config));
clear oc;

wd = run_example('benchmark');

fH = fopen(fullfile(wd,'benchmark_result.txt'), 'r');
oc = onCleanup(@() fclose(fH));
simOutput = jsondecode(fgets(fH));
loop_time = simOutput.time;

end