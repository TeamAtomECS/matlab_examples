%% Cross section
% 

wd = run_example('cross_section');

output_v = utils.read_output(fullfile(wd,'vel.txt'));
vs = output_v(end).vec;
output_s = utils.read_output(fullfile(wd,'scattered.txt'), 'Format', '%f,');
scat = output_s(end).vec;

% Convert velocity to units of Gamma.
gamma = 2*pi*6e6; % MHz
lambda = 780e-9; %m, Rubidium D2
k = 2*pi/lambda;
doppler = vs(:,1) * k;


% timestep is 1us, so number of photons per timestep = rate in MHz.
plot(doppler/gamma, 1e6*scat(:,1)/gamma, 'k-');

xlabel('$\delta/\Gamma$', 'Interpreter', 'latex');
ylabel('scattering rate$/\Gamma$', 'Interpreter', 'latex');
box on;
set(gcf, 'Color', 'w');
grid on;
set(gca, 'GridLineStyle', ':');
ylim([0, 0.55]);

hold on;
% Calculate expression from Chris' book

IoverIsat = 1;
rate_th = 0.5 * IoverIsat ./ (1 + IoverIsat + 4*(doppler/gamma).^2);
plot(doppler/gamma, rate_th, '.r', 'MarkerSize', 7);
hold off;

%%
% 
hb = 6.63e-34/(2*pi);
c = 3e8;
omega = (2*pi*c/lambda);

Isat = hb * omega^3 * gamma / (8 * pi^2 * c^2)