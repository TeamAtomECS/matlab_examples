function wd = run_example(example)
%RUN_EXAMPLE Runs the given example from the atomecs repository
system(sprintf('cd atomecs & cargo run --example %s --release', example));
wd = fullfile(cd, 'atomecs');

end