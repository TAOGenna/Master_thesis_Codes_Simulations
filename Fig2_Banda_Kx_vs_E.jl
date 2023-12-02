using LinearAlgebra
using Plots
using ColorSchemes


sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
iden =[1. 0.; 0. 1.]
# Define Hamiltonian
function H(gam, k0, kx, ky, kz, m, tx, t, gamx)
    return eigvals(gam*(cos(kx)-cos(k0))*iden-(m*(2-cos(ky)-cos(kz))
        +2*tx*(cos(kx)-cos(k0)))*sig_x-2*t*sin(ky)*sig_y-2*t*sin(kz)*sig_z
        -gamx*(cos(3*kx)-cos(3*k0))*sig_x )
end


H_max=pi
plot(xlims=(-H_max,H_max),framestyle=:origin)
#title!(L"$\gamma=2t$")
xlabel!(L"$k_x$")
ylabel!(L"$E$")
I=300
#z = zeros(length(x), length(y), 2)
vals_kx = collect(-H_max:H_max/I:H_max)
vals_kz = [0.]
#parametros para el calculo del Hamiltoniano
k0=pi/2
ky=0
t=1
gamma=t*1.5 # 0,t,1.5t
gammax=t/2.
tx=1 #tx=t
m=2 #m=2*t
#-----------------------------------------------------
for i in 1:length(vals_kx)
    for j in 1:length(vals_kz)
        espec=H(gamma, k0, vals_kx[i], ky, vals_kz[j], m, tx, t, gammax)
        display(scatter!([vals_kx[i]], [espec[1]], markercolor = "blue", markerstrokecolor = "blue", markersize = 1,legend=false))
        display(scatter!([vals_kx[i]], [espec[2]], markercolor = "red", markerstrokecolor = "red", markersize = 1,legend=false))
    end
end
#Plot the Weyl nodes at -k and k
display(scatter!([k0], [0.], markercolor = RGB(0.5, 1, 0), markerstrokecolor = RGB(0.5, 1, 0), markersize = 4,legend=false))
display(scatter!([-k0], [0.],markercolor = RGB(0.5, 1, 0), markerstrokecolor = RGB(0.5, 1, 0), markersize = 4,legend=false))