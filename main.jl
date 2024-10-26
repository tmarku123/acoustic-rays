include("functions.jl");

## Define domain ##
x_max = 200;
z_max = 200;
no_layers = 7;
scaling_factor = 1;
domain = create_domain(x_max,z_max,no_layers,scaling_factor);

## Define Simulations parameters ## 
dt = 0.1;
tmax = 50;
time_vector = collect(0:dt:Int(round(tmax)));
parameters = Parameters(time_vector);
no_rays = 100;

## Set Initial Conditions ## 
position = [z_max/2,x_max/2];
amplitude = 25.0;

## Single Ray ## 
#angle = 10.0;
#initial_conditions = set_initial_conditions(position,angle,amplitude);
#rays = initialise_rays(no_rays,initial_conditions);
#simulation = initialise_simulation(parameters,domain,rays);
#simulation = compute_ray_simulation(simulation::Simulation);
#animation = animate_rays(simulation.data,domain,10)

## Multiple Rays (Form a wavefront) ## 
dθ = 0.25;
angles = collect(0:dθ:360)[1:end-1];
simulations = Vector{Simulation}(undef,size(angles,1));
for i = 1:size(angles)[1] 
    initial_conditions = set_initial_conditions(position,angles[i],amplitude);
    rays = initialise_rays(no_rays,initial_conditions);
    simulation = initialise_simulation(parameters,domain,rays);
    simulation = compute_ray_simulation(simulation::Simulation);
    simulations[i] = simulation
end

wavefronts = obtain_wavefronts(simulations)
animate_wavefronts(simulations,wavefronts,8)







