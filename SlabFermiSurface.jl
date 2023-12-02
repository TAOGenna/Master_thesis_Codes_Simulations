using LinearAlgebra
using Plots
using LaTeXStrings

sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
iden =[1. 0.; 0. 1.]
# Define Hamiltonian
function H(gam, k0, kx, ky, kz, m, tx, t, L)
    H=zeros(2*L,2*L)
    mat=zeros(L,L)
    for i in range(1,L)
        mat[i,i]=1
        H=H+kron(mat,gam*(cos(kx)-cos(k0))*iden-(m*(2-cos(kz))+2*tx*(cos(kx)-cos(k0)))*sig_x-2*t*sin(kz)*sig_z)
        mat[i,i]=0
    end
    for i in range(1,L-1)
        mat[i,i+1]=1
        H=H+kron(mat,(m/2)*sig_x+im*t*sig_z)
        mat[i,i+1]=0
        mat[i+1,i]=1
        H=H+kron(mat,(m/2)*sig_x-im*t*sig_z)
        mat[i+1,i]=0
    end
    return eigvals(H)
end


H_max=1
plot(xlims=(-H_max,H_max),ylims=(-H_max,H_max),framestyle=:origin)
title!(L"$FermiSurface\quad E=0.2*t\quad \gamma=0$")
xlabel!(L"$k_x$")
ylabel!(L"$k_z$")
I=70
vals_kx = collect(-H_max:H_max/I:H_max)
vals_kz = collect(-H_max:H_max/I:H_max)
z = zeros(length(vals_kx), length(vals_kz), 2)
#parametros para el calculo del Hamiltoniano
gamma=0
k0=pi/2
ky=0
t=1
tx=1 #tx=t
m=2 #m=2*t
Num_Sitios=20
eps=10^(-9)
Energia=0
for i in 1:length(vals_kx)
    for j in 1:length(vals_kz)
        espec=H(gamma, k0, vals_kx[i], ky, vals_kz[j], m, tx, t, Num_Sitios)
        #display(espec)
        for xx in 1:length(espec)
            if espec[xx]<Energia+eps && espec[xx]>Energia-eps
                display(scatter!([vals_kx[i]], [vals_kz[j]], markercolor = "blue", markerstrokecolor = "blue", markersize = 0.5,legend=false))
            end
        end
    end
end
display(scatter!([k0], [0.], markercolor = RGB(0.5, 1, 0), markerstrokecolor = RGB(0.5, 1, 0), markersize = 4,legend=false))
display(scatter!([-k0], [0.],markercolor = RGB(0.5, 1, 0), markerstrokecolor = RGB(0.5, 1, 0), markersize = 4,legend=false))