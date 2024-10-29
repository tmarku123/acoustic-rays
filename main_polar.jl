include("functions/polar_simulation.jl");

r_max = 1.0;
theta_max = 360;
no_layers = 3;
min_amp = 0;

velocity = collect(range(1,2,no_layers))
density = collect(range(3,1,no_layers))
properties = Vector{Property}(undef,no_layers)

for i = 1:no_layers
    properties[i] = Property(velocity[i],density[i])
end

polar_domain = PolarDomain(r_max,theta_max,no_layers,properties,min_amp)

dt = 0.001;
tmax = 1;
time_vector = collect(0:dt:Int(round(tmax)));
parameters = Parameters(time_vector);
no_rays = 50;

## Set Initial Conditions ## 
pos = [0, r_max-0.1];
theta = -90;
amp = 1000; 

## Singular Ray 
initial_conditions = set_initial_conditions(pos,theta,amp);
rays = initialise_rays(no_rays,initial_conditions);
simulation = initialise_simulation(parameters,polar_domain,rays);
simulation = compute_polar_rays(simulation);
#animate_polar_rays(simulation.data,polar_domain,5,no_layers)
no_rays = get_number_of_rays_used(simulation)
println("No. rays used in Multiple Rays: ", no_rays)

## Multiple Rays (Form a wavefront) ## 
dθ = 1
angles = collect(0:dθ:theta_max)[1:end-1];
simulations = Vector{Simulation}(undef,size(angles,1));
for i = 1:size(angles)[1] 
    angle = Float64(angles[i])
    initial_conditions = set_initial_conditions(pos,angle,amp);
    rays = initialise_rays(no_rays,initial_conditions);
    simulation = initialise_simulation(parameters,polar_domain,rays);
    simulation = compute_polar_rays(simulation::Simulation);
    simulations[i] = simulation
end


scaling = 125
wavefronts = obtain_polar_wavefronts(simulations,scaling,blur=true);

include("functions/visuals.jl")
animation_time = 7;
darkmode_animate_polar_wavefronts(simulations,wavefronts,animation_time,no_layers,scaling)


