using Plots 
using Images
include("utilities.jl")


## Cartesian animations
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
    gif(anim, "polar_wavefront_animation.gif", fps = desired_framerate)
end


## Polar animations

function draw_polar_domain(domain::PolarDomain,scaling::Int64)
    R = domain.r_max
    theta = domain.theta_max
    no_layers = domain.no_layers
    r = collect(range(0,R,no_layers+1))[2:end]
    r[end] = R;
    no_divisions = 100
    theta = collect(range(0,theta,no_divisions))
    x = zeros(no_divisions,no_layers)
    z = zeros(no_divisions,no_layers)
    if scaling == 1
        for j = 1:no_layers
            rad = r[j]
            for i = 1:no_divisions
                x[i,j] = rad*cosd(theta[i]) + R
                z[i,j] = rad*sind(theta[i]) + R
            end
        end
    else
        nx = Int(round(r_max*2*scaling))
        for j = 1:no_layers
            rad = r[j]
            for i = 1:no_divisions
                x1 = rad*cosd(theta[i])
                x1 = get_scaled_points(x1,R,nx)
                x[i,j] = x1
                z1 = rad*sind(theta[i])
                z1 = get_scaled_points(z1,R,nx)
                z[i,j] = z1 
            end
        end
    end
    fig = plot(x,z,label=false,color=:black,yflip=false)
    return fig 
end

function animate_polar_rays(data,domain::PolarDomain,
    animation_time::Int64,no_layers::Int64)
    frame_rate = 30
    total_frames = animation_time*frame_rate
    max_steps = size(data,1)

    x = data[:,2,:]
    z = data[:,1,:]

    fig = draw_polar_domain(domain,1)

    anim = @animate for i = 1:total_frames
        indc = Int(round((i/total_frames)*max_steps))
        plot(fig,
            x[1:indc,:],
            z[1:indc,:],
            yflip = true,
            labels=false
            )

    end
    display(anim)
    gif(anim, "ray_animation.gif", fps = frame_rate); 
end;

function animate_polar_wavefronts(simulations::Vector{Simulation},
    wavefronts::Vector{Wavefront},animation_time::Int64,no_layers::Int64,
    scaling::Int64)
    domain = simulations[1].domain
    nt = size(wavefronts,1)
    desired_framerate = 30
    frames_req = animation_time*desired_framerate
    my_cmap = cgrad([RGBA(1, 1, 1,0), RGBA(0.15,0.15,0.15,0.9)], scale=:linear)
    anim = @animate for i = 1:frames_req
        i = clamp(Int(round((i/frames_req)*nt)), 1, nt)
        draw_polar_domain(domain,scaling)
        heatmap!(wavefronts[i].snapshot,yflip=false,colorbar=false,clims=(0,5),color=my_cmap)
    end

    display(anim)
    gif(anim, "polar_wavefront_animation.gif", fps=desired_framerate)
end

function darkmode_animate_polar_wavefronts(simulations::Vector{Simulation},
    wavefronts::Vector{Wavefront},animation_time::Int64,no_layers::Int64,
    scaling::Int64)
    domain = simulations[1].domain
    r_max = domain.r_max
    min_amp = domain.minimum_amplitude
    nx = size(wavefronts[1].snapshot,1)
    nt = size(wavefronts,1)
    desired_framerate = 30
    frames_req = animation_time*desired_framerate
    my_cmap = cgrad([RGBA(0.3,0.3,0.3,0.7), RGBA(1,1,1,0.8)], scale=:linear)
    tickss = collect(range(-r_max,r_max,5))
    anim = @animate for i = 1:frames_req
        i = clamp(Int(round((i/frames_req)*nt)), 1, nt)
        fig = draw_polar_domain(domain, scaling)

        # Add the heatmap for the current frame
        heatmap!(fig, wavefronts[i].snapshot, yflip=false, colorbar=false,
        nan_color=:white,
        clims=(min_amp, 100), color=my_cmap,
        xticks=(range(1, nx, length=5),tickss),
        yticks=(range(1, nx, length=5),tickss),
        framestyle=:box,               # Sets a box frame to enclose the plot
         tickfontsize=10,               # Increases font size for readability
         margin=5Plots.mm   
)
    end
    display(anim)
    gif(anim, "polar_wavefront_animation.gif", fps=desired_framerate)
end
