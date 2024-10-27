include("structs.jl")
using LinearAlgebra
using Plots 
using Images

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

function initialise_simulation(parameters::Parameters,domain::Domain,rays::Vector{Ray})
    nt = size(parameters.time,1)
    no_rays = size(rays,1)
    data = Array{Float64}(undef,nt,3,no_rays).*NaN
    simulation = Simulation(parameters,domain,rays,data)
    return simulation
end

function get_zone_id(ray::Ray,domain::Domain)
    (z,x) = ray.position
    (nz,nx) = size(domain.domain)
    scaling_factor = domain.scaling_factor
    x = clamp(Int(round(x*scaling_factor)),1,nx);
    z = clamp(Int(round(z*scaling_factor)),1,nz);
    id = domain.domain[z,x]
    return Int(id)
end;

function get_velocity(ray,domain::Domain)
    zone_id = get_zone_id(ray,domain)
    velocity = domain.zones[zone_id].properties.wave_velocity
    return velocity
end;

function check_boundary(ray::Ray,domain::Domain)
    updated_ray = update_ray(ray,domain)
    vel1 = get_velocity(ray,domain)
    vel2 = get_velocity(updated_ray,domain)
    if vel1 != vel2
        return true
    end
end;

function compute_boundaries(domain::Matrix{Float64})
    nz,nx = size(domain)
    boundaries = zeros(nz,nx)
    for i = 2:nz-1, j = 2:nx-1
        submatrix = domain[(i-1:i+1),(j-1:j+1)]
        if sum(submatrix) != domain[i,j]*9
            boundaries[i,j] = 1
        end
    end
    return boundaries
end;

function is_in_domain(ray::Ray,domain::Domain)
    position = ray.position
    scaling_factor = domain.scaling_factor
    (nz,nx) = size(domain.domain)./scaling_factor
    (z,x) = position 
    if (0 <= z < nz) && (0 <= x < nx)
        return true 
    else
        return false
    end
end;

function get_active_rays(rays::Vector{Ray})
    status = [ray.status for ray in rays] .== "active"
    return status
end;

function get_empty_ray(rays::Vector{Ray})
    status = [ray.status for ray in rays] .== "inactive"
    indc = findfirst(status .== 1)
    return indc 
end;

function complete_ray(ray::Ray)
    pos = ray.position
    dir = ray.direction
    amp = ray.amplitude
    length = ray.length
    IC = ray.initial_conditions
    completed_ray = Ray(pos,dir,amp,length,IC,"completed")
    return completed_ray
end;

function update_ray(ray::Ray,domain::Domain)
    pos = ray.position
    dir = ray.direction
    vel = get_velocity(ray,domain)
    ds = dir.*vel*dt 
    new_pos = pos .+ ds 
    length = ray.length + norm(ds)
    IC = ray.initial_conditions
    a0 = IC.amplitude
    new_amp = a0*sqrt(norm(ds)/length)

    return Ray(new_pos,dir,new_amp,length,IC,"active")
end;

function reflect_ray(ray::Ray,domain::Domain)
    pos = ray.position 
    dir = ray.direction
    amp = ray.amplitude
    n = [1, 0]
    outgoing = dir - 2*dot(dir,n)*n 
    angle = acosd(dir[1])
    IC = InitialConditions(pos,angle,amp)
    length = 0
    damping = 0.8
    new_amp = amp*damping
    reflected_ray = Ray(pos,outgoing,new_amp,length,IC,"active")
    return reflected_ray
end;

function save_ray(ray::Ray,data::Array{},ray_index::Int64,time_index::Int64)
    (z,x) = ray.position 
    amp = ray.amplitude
    data[time_index,:,ray_index] .= [z,x,amp]
    return data
end;

function snells_law(ray::Ray,domain::Domain)
    pos = ray.position
    dir = ray.direction
    amp = ray.amplitude
    vel1 = get_velocity(ray,domain)
    updated_ray = update_ray(ray,domain)
    vel2 = get_velocity(updated_ray,domain)
    (dz,dx) = dir
    hyp = sqrt(dx^2 + dz^2)
    sinθ1 = dx/hyp
    sinθ2 = (vel2/vel1)*sinθ1

    if -1 <= sinθ2 <= 1
        θ2 = asind(sinθ2)
        new_direction = [cosd(θ2)*sign(dir[1]), sind(θ2)]
        length = updated_ray.length
        a0 = ray.initial_conditions.amplitude
        new_pos = updated_ray.position
        new_amp = updated_ray.amplitude
        IC = InitialConditions(new_pos,θ2,a0)
        snells_ray = Ray(new_pos,new_direction,new_amp,length,IC,"active")

        return snells_ray
    end
    return complete_ray(ray)
end;

function compute_ray_simulation(simulation::Simulation)
    parameters = simulation.parameters
    time = parameters.time
    dt = time[2]-time[1]
    domain = simulation.domain
    rays = simulation.rays
    data = simulation.data

    first_ray = rays[1]
    data = save_ray(first_ray,data,1,1)

    for i = 1:size(time)[1]-1
        status = get_active_rays(rays)
        for j = 1:size(status)[1]
            if status[j] == 0
                continue
            end
            ray = rays[j]
            if is_in_domain(ray,domain)
                if check_boundary(ray,domain) == true 
                    empty_indc = get_empty_ray(rays)
                    if !isnothing(empty_indc)
                        reflected_ray = reflect_ray(ray,domain)
                        data = save_ray(reflected_ray,data,empty_indc,i)
                        reflected_ray = update_ray(reflected_ray,domain)
                        data = save_ray(reflected_ray,data,empty_indc,i+1)
                        rays[empty_indc] = reflected_ray
                    end
                    snells_ray = snells_law(ray,domain)
                    data = save_ray(snells_ray,data,j,i+1)
                    rays[j] = snells_ray 
                else
                    updated_ray = update_ray(ray,domain)
                    data = save_ray(updated_ray,data,j,i+1)
                    rays[j] = updated_ray
                end
            else
                completed_ray = complete_ray(ray)
                rays[j] = completed_ray
            end
        end
    end
    return Simulation(parameters,domain,rays,data)
end;

function get_number_of_rays_used(simulation::Simulation)
    data = simulation.data 
    for i = 1:size(data,3)
        data_slice = data[:,:,i]
        if maximum(.!isnan.(data_slice)) == false
            return i 
        end
    end
end

function obtain_wavefronts(simulations::Vector{Simulation})
    time = simulations[1].parameters.time
    nt = size(time,1)
    scaling = 1
    (nz,nx) = Int.(size(simulations[1].domain.domain).*scaling)
    wavefronts = Vector{Wavefront}(undef,nt)    
    for i = 1:nt
        t = time[i]
        snapshot = zeros(nz,nx)
        for j = 1:size(simulations)[1]
            data = simulations[j].data[i,:,:]
            indc = .!isnan.(data)
            data = data[indc]
            if size(data,1) > 0 
                cols = Int(size(data,1)/3)
                data = reshape(data,(3,cols))
                z = Int.(round.(data[1,:].*scaling))
                x = Int.(round.(data[2,:].*scaling))
                amp = data[3,:]
                x = clamp.(x,1,nx)
                z = clamp.(z,1,nz)
                for k = 1:cols
                    if isnan(snapshot[z[k],x[k]])  
                        snapshot[z[k],x[k]] = 0
                    else
                        snapshot[z[k],x[k]] += amp[k]
                    end
                end
            end
        end
        wavefront = Wavefront(snapshot,t)
        wavefronts[i] = wavefront
    end

    return wavefronts
end

function animate_rays(data,domain::Domain,animation_time::Int64,no_layers::Int64)
    frame_rate = 30
    total_frames = animation_time*frame_rate
    max_steps = size(data,1)
    (z_max,x_max) = size(domain.domain)
    changes = collect(range(0,z_max,no_layers+1))
    changes[1] = NaN; changes[end] = NaN
    x_line = ones(x_max,size(changes,1)).*collect(1:x_max)
    y_line = ones(x_max,size(changes,1)).*changes'

    x = data[:,2,:]
    z = data[:,1,:]

    (nz,nx) = size(domain.domain)
    x_line = ones(x_max,size(changes,1)).*collect(1:x_max)
    y_line = ones(x_max,size(changes,1)).*changes'
    anim = @animate for i = 1:total_frames
        indc = Int(round((i/total_frames)*max_steps))
        plot(x_line,y_line,labels=false,color=:black,yflip=true,xlims=(0,x_max),ylims=(0,z_max))
        plot!(
            x[1:indc,:],
            z[1:indc,:],
            yflip = true,
            labels=false
            )

    end
    display(anim)
    gif(anim, "ray_animation.gif", fps = frame_rate); 
end;

function animate_wavefronts(simulations::Vector{Simulation},wavefronts::Vector{Wavefront},animation_time::Int64,no_layers::Int64)
    domain_map = simulations[1].domain.domain
    (z_max,x_max) = size(domain_map)
    changes = collect(range(0,z_max,no_layers+1))
    changes[1] = NaN; changes[end] = NaN
    x_line = ones(x_max,size(changes,1)).*collect(1:x_max)
    y_line = ones(x_max,size(changes,1)).*changes'
    nt = size(wavefronts,1)
    desired_framerate = 30
    frames_req = animation_time*desired_framerate
    colors = cgrad([:transparent, RGB(0.0, 0.0, 0.0)], [0, 1])  # Change RGB values for different dark colors

    anim = @animate for i = 1:frames_req
        i = clamp(Int(round((i/frames_req)*nt)), 1, nt)
        #heatmap(domain_map,alpa=0.1,yflip=true)
        plot(x_line,y_line,labels=false,color=:black,yflip=true)
        heatmap!(wavefronts[i].snapshot,yflip=true,colorbar=false,clims=(0,5)
        ,color=colors,title = "Time: $(wavefronts[i].time)")
    end

    display(anim)
    gif(anim, "wavefront_animation.gif", fps=desired_framerate)
end


function darkmode_animate_wavefronts(simulations::Vector{Simulation}, wavefronts::Vector{Wavefront}, animation_time::Int64, no_layers::Int64)
    domain_map = simulations[1].domain.domain
    (z_max, x_max) = size(domain_map)
    changes = collect(range(0, z_max, no_layers + 1))
    changes[1] = NaN
    changes[end] = NaN
    x_line = ones(x_max, size(changes, 1)) .* collect(1:x_max)
    y_line = ones(x_max, size(changes, 1)) .* changes'
    nt = size(wavefronts, 1)
    desired_framerate = 30
    frames_req = animation_time * desired_framerate
    colors = cgrad([RGBA(0.15, 0.15, 0.15, 0.0), RGBA(1.0, 1.0, 1.0)], [0, 1])  # Heatmap with transparent dark to bright white

    anim = @animate for i = 1:frames_req
        i = clamp(Int(round((i / frames_req) * nt)), 1, nt)

        plot(
            x_line, y_line,
            labels = false,
            color = :white,
            yflip = true,
            background_color = RGB(0.15, 0.15, 0.15),  
            xlims = (0, x_max),
            ylims = (0, z_max)
        )

        heatmap!(
            wavefronts[i].snapshot,
            yflip = true,
            colorbar = false,
            clims = (0, 5),
            color = colors,
            alpha = 0.8, 
            title = "Time: $(wavefronts[i].time)"
        )
    end

    display(anim)
    gif(anim, "wavefront_animation.gif", fps = desired_framerate)
end

