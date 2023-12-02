using LinearAlgebra
using Plots

function hamiltonian(kx, t, delta, tri)
    return 2*t*cos(kx/2)*[0 1; 1 0] - 2*delta*sin(kx/2)*[0 -im; im 0] + tri*[1 0; 0 -1]
end

function berry_connection(kx, t, delta, tri)
    h = hamiltonian(kx, t, delta, tri)
    eigs = eigvecs(h)
    deigs = [gradient(eigs[:,i]) for i in 1:size(eigs)[2]]
    return eigs'*deigs
end

function berry_curvature(kx_range, t, delta, tri)
    A = [berry_connection(kx, t, delta(l), tri(l)) for (kx,l) in zip(kx_range,linspace(0.,2*pi,length=length(kx_range)))]
    B = [gradient(A[i], kx_range)[3] for i in 1:length(A)]
    return B
end

function zak_phase(lambdas::Array{Float64}, t::Float64)
    kx_range = range(0., stop=2*pi,length=1000)
    delta(x) = -0.4*abs(t)*sin(x)
    tri(x) = -0.4*abs(t)*cos(x)
    B = [sum(berry_curvature(kx_range=kx_range,t=t,delta=delta(l),tri=tri(l)))/(2*pi) for l in lambdas]
    return B
end

t = -1
lambdas = range(0., stop=2*pi,length=1000)
zak_phases = zak_phase(lambdas,t)

plot(lambdas,zak_phases,xlabel="lambda",ylabel="Zak phase",title="Zak phase vs lambda")
