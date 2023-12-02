using LinearAlgebra
using Plots
using ColorSchemes


sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
iden =[1. 0.; 0. 1.]
# Define Hamiltonian
function H(gam, k0, kx, ky, kz, m, tx, t)
    return eigvals(gam*(cos(kx)-cos(k0))*iden-(m*(2-cos(ky)-cos(kz))+2*tx*(cos(kx)-cos(k0)))*sig_x-2*t*sin(ky)*sig_y-2*t*sin(kz)*sig_z)
end


H_max=pi
plot(xlims=(-H_max,H_max))
title!(L"$apilando\quad bandas,\quad\gamma=0,\quad k_y=0,\quad k_x=\{0,2.5,\pi\}$")
xlabel!(L"$k_z$")
ylabel!(L"$Energia$")
I=20
vals_kx = [0,2.5,pi]
vals_kz = collect(-H_max:H_max/I:H_max)
#parametros para el calculo del Hamiltoniano
gamma=0
k0=pi/2
ky=0
t=1
tx=1 #tx=t
m=2 #m=2*t
Num_Sitios=20
eps=10^(-1)
colores=["blue","red","green"]
curvas=[L"$k_x=0$",L"$k_x=2.5$",L"$k_x=\pi$"]
#-----------------------------------------------------
eps=10^(-6)
for indx in 1:length(vals_kx)
    legend_added=false
    for valz in vals_kz
       espec=H(gamma, k0, vals_kx[indx], ky, valz, m, tx, t, Num_Sitios)
        for energy in espec
            if !legend_added
                display(scatter!([valz], [energy], markercolor = colores[indx], markerstrokecolor = colores[indx], markersize = 1,label=curvas[indx]))
                #display(scatter!([mu], [suma], markercolor = esc, markerstrokecolor = esc, markersize = tama[letsgo],label=curvas[letsgo]))
                legend_added=true
            else
                display(scatter!([valz], [energy], markercolor = colores[indx], markerstrokecolor = colores[indx], markersize = 1,label=""))
                #display(scatter!([mu], [suma], markercolor = esc, markerstrokecolor = esc, markersize = tama[letsgo],label=""))
            end
            #display(scatter!([valz], [energy], markercolor = colores[indx], markerstrokecolor = colores[indx], markersize = 0.5,legend=false))
        end
    end
end

#Plot the Weyl nodes at -k and k
display(scatter!([k0], [0.], markercolor = RGB(0.5, 1, 0), markerstrokecolor = RGB(0.5, 1, 0), markersize = 4,legend=false))
display(scatter!([-k0], [0.],markercolor = RGB(0.5, 1, 0), markerstrokecolor = RGB(0.5, 1, 0), markersize = 4,legend=false))