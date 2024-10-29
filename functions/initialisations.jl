include("utilities.jl")

function create_domain(x_max::Int64,z_max::Int64,no_layers::Int64,scaling_factor::Int64)
    nx = x_max*scaling_factor;
    nz = z_max*scaling_factor;
    changes = collect(range(0,z_max,no_layers+1))
    current_domain = Matrix{Float64}(undef,nz,nx);

    for i = 1:no_layers
        a = Int(round(clamp(changes[i]*scaling_factor,1,z_max*scaling_factor)))
        b = Int(round(clamp(changes[i+1]*scaling_factor,1,z_max*scaling_factor)))
        current_domain[a:b,:] .= i
    end

    zones = Vector{Zone}(undef,no_layers)
    Δv = 1;
    Δρ = 1;

    for i = 1:no_layers
        properties = MaterialProperties(Δv*i,Δρ*i)
        zone = Zone(i,properties)
        zones[i] = zone
    end

    boundary_indcs = compute_boundaries(current_domain)
    normals = Matrix{Vector{Float64}}(undef,nz,nx)

    for i = 1:nz, j = 1:nx
        if boundary_indcs[i,j] == 1
            normals[i,j] = [1, 0]
        else
            normals[i,j] = [0, 0]
        end
    end

    boundaries = Boundaries(boundary_indcs,normals)
    domain = Domain(zones,boundaries,current_domain,scaling_factor)
    return domain
end;

function set_initial_conditions(position::Vector{},angle,amplitude)
    initial_conditions = InitialConditions(position,angle,amplitude)
    return initial_conditions
end

function initialise_rays(no_rays::Int64,
    initial_conditions::InitialConditions,)

    pos = initial_conditions.position
    angle = initial_conditions.angle
    dir = [cosd(angle), sind(angle)]
    amp = initial_conditions.amplitude
    initial_ray = Ray(pos,dir,amp,0,initial_conditions,"active")
    rays = Vector{Ray}(undef, no_rays)
    rays[1] = initial_ray
    dummy_IC = InitialConditions([NaN,NaN],NaN,NaN)
    dummy_ray = Ray([NaN,NaN],[NaN,NaN],NaN,NaN,dummy_IC,"inactive")

    for i = 2:no_rays
        rays[i] = dummy_ray
    end

    return rays
end

function initialise_simulation(parameters::Parameters,domain::Union{Domain,PolarDomain},rays::Vector{Ray})
    nt = size(parameters.time,1)
    no_rays = size(rays,1)
    data = Array{Float64}(undef,nt,3,no_rays).*NaN
    simulation = Simulation(parameters,domain,rays,data)
    return simulation
end
