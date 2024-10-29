include("functions/cartesian_simulation.jl");
include("functions/visuals.jl")

## Define domain ##
x_max = 400;
z_max = 400;
no_layers = 3;
scaling_factor = 1;
domain = create_domain(x_max,z_max,no_layers,scaling_factor);

## Define Simulations parameters ## 
dt = 0.1;
tmax = 50.0;
time_vector = collect(0:dt:Int(round(tmax)));
parameters = Parameters(time_vector);
no_rays = 20;

## Set Initial Conditions ## 
pos = [z_max/2,x_max/2];
amp = 10^2; 


## Single Ray ## 
# Run this at 0 degrees to get a good no_rays to run for multiple case 
angle = 0;
initial_conditions = set_initial_conditions(pos,angle,amp);
rays = initialise_rays(no_rays,initial_conditions);
simulation = initialise_simulation(parameters,domain,rays);
simulation = compute_ray_simulation(simulation::Simulation);
#animation = animate_rays(simulation.data,domain,7,no_layers)
no_rays = get_number_of_rays_used(simulation)
println("No. rays used in Multiple Rays: ", no_rays)
## Multiple Rays (Form a wavefront) ## 
dθ = 0.5;
angles = collect(0:dθ:360)[1:end-1];
simulations = Vector{Simulation}(undef,size(angles,1));
for i = 1:size(angles)[1] 
    angle = Float64(angles[i])
    initial_conditions = set_initial_conditions(pos,angle,amp);
    rays = initialise_rays(no_rays,initial_conditions);
    simulation = initialise_simulation(parameters,domain,rays);
    simulation = compute_ray_simulation(simulation::Simulation);
    simulations[i] = simulation
end

wavefronts = obtain_wavefronts(simulations);
animation_time = 10;
darkmode_animate_wavefronts(simulations,wavefronts,animation_time,no_layers)
