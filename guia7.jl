#=using LinearAlgebra
using Plots


sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
iden =[1. 0.; 0. 1.]

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
    A = [berry_connection(kx, t, delta, tri) for (kx,l) in zip(kx_range,linspace(0.,2*pi,length=length(kx_range)))]
    B = [gradient(A[i], kx_range)[3] for i in 1:length(A)]
    return B
end

function Zak(vec,t)
	pasos=100
	sol=[]
	kx_range=collect(0:2*pi/pasos:2*pi)
	delta(x) = -0.4*abs(t)*sin(x)
    tri(x) = -0.4*abs(t)*cos(x)
	for x in vec
		suma_berry=sum(berry_curvature(kx_range=kx_range,t=t,delta=delta(x),tri=tri(x)))
		sum=sum/(2*pi)
		append!(sol,(x,suma_berry))
	end
	return sol
end



H_max=2*pi
I=100
t=-1
lambdas=collect(0:H_max/I:H_max)
Zak_phases=Zak(lambdas,t)
=#
using LinearAlgebra
using SparseArrays


Nx = 10 # Cantidad de celdas en la dirección x
a = 1 # parámetro de red
masa = 1 # Según el enunciado debo fijar m=1
discretizacion = 100 # es una discretización de ky para graficar la energía en función de ky

# Hamiltoniano
auxiliar(ky, m) = (-2 + m + cos(ky))*[0 0;0 1] + sin(ky)*[0 -im;im 0]
listaPerturbacionC = [0.2*[1 0;0 0], 0.2*[0 1;0 0], zeros(2,2)]

function primer(ky, m, perturbacionC)
	zer=zeros(Nx,Nx)
	for i in 1:Nx
		zer[i,i]=1
	end
    return kron(zer, kron([1 0;0 1], auxiliar(ky, m)) + kron([0 1;0 0], perturbacionC))
end

function segundo(ky, m)
    return 0.5*kron(sparse(Band([1,2]), [1 for i in range(1,Nx-1)]), kron([1 0;0 -1], [0 1;1 0]))
end

function tercer(ky, m)
    return 0.5im*kron(sparse(Band([1,2]), [im for i in range(1,Nx-1)]), kron([1 0;0 -1], [0 -im;im 0]))
end

function segundoConjugado(ky, m)
    return segundo(ky,m)'
end

function tercerConjugado(ky, m)
    return tercer(ky,m)'
end

function H(ky, m, perturbacionC)
    return primer(ky,m,perturbacionC) + segundo(ky,m) + tercer(ky,m) + segundoConjugado(ky,m) + tercerConjugado(ky,m)
end

eigs = eigvals(H(pi, masa, listaPerturbacionC[1]))
plot(eigs)