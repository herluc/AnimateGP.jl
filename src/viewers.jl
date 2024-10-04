using Plots

function loop_frames_2d(frames,x,y,delay=0.01)
    f = size(frames)[2]
    for i in collect(range(1,f))
        c = frames[:,i]
        grid_data = reshape(c, (length(unique(x)), length(unique(y))))
        plt = heatmap(x,y,grid_data,clim=(-3,3),color=cgrad(:default))

        display(plt)
        sleep(delay)
    
    end
end
