module AnimateGP

using LinearAlgebra
using KernelFunctions

export compute_frames

# Write your package code here.
function sample_standard_normal_frames(d,n)
    x = randn((1,d)) #starting sample
    r = sqrt(sum(x.^2))
    x = x./r #project onto sphere
    t = randn((1,d)) #sample tangent direction
    t = t - ((t'*x) * x')' #orthogonalize by Gram-Schmidt
    t = t ./ sqrt(sum(t.^2)) # standardise
    s = collect(range(0,2*pi,n+1))
    s = s[1:end-1] # space to span
    t = s * t # span linspace in direction of t
    X = r .* exp_map(x,t) # project onto sphere, re-scale

    return X
end #function


function exp_map(mu,E)
    D = size(E',1)
    theta = sqrt.(sum((E'.^2),dims=1))
    M = mu' * cos.(theta) .+ (E' .* repeat(sin.(theta)./theta, D, 1))
    if any(x->abs.(x)<=(1e-7), theta)
        a_vec = findall(x->abs.(x)<=(1e-7), theta)
        for a in a_vec
            M[:, a] = mu
        end # for
    end # if
    # M (:,abs (theta) <= 1e-7) = mu;
    return M
end #function


function compute_frame(L, mean_vec, standards, i)
    c = mean_vec + L'*standards[:,i]

    return c
end #function


function compute_frames(cov_mat, mean_vec, n_frames)
    d = size(mean_vec)[1]
    standards = sample_standard_normal_frames(d,n_frames)
    frames = deepcopy(standards)
    L = cholesky(cov_mat + (I*1.0e-8)).U;
    for i in collect(range(1,n_frames))
        frames[:,i] .= compute_frame(L,mean_vec,standards,i)
    end
    return frames
end #function


end #module
