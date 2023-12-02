using LinearAlgebra
using Plots
using LaTeXStrings

sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
iden =[1. 0.; 0. 1.]

# Define Hamiltonian
function H(ky,m,L)
    mat=zeros(L,L)
    for i in range(1,L)
        mat[i,i]=1
    end
    primero=kron(mat,kron(iden,(-2+m+cos(ky))*sig_z+sin(ky)*sig_y))
    mat=zeros(L,L)
    for i in range(1,L-1)
    	mat[i,i+1]=1
    end
    segundo=kron(mat, (im/2)*kron(sig_z,sig_x)+(1/2)*kron(iden,sig_z))
    mat=zeros(L,L)
    for i in range(1,L-1)
    	mat[i+1,i]=1
    end
    tercero=kron(mat, (1/2)*kron(iden,sig_z)-(im/2)*kron(sig_z,sig_x))
    asd=primero+segundo+tercero
    return asd
end


H_max=pi
plot(xlims=(-H_max,H_max),ylims=(-H_max,H_max),framestyle=:origin)
#title!(L"$FermiSurface\quad E=0.2*t\quad \gamma=0$")
xlabel!(L"$k_y$")
ylabel!(L"$E$")
I=40
L=2
vals_ky = collect(-H_max:H_max/I:H_max)
vals_kz = collect(-H_max:H_max/I:H_max)
for i in 1:length(vals_ky)
    hammy=H(vals_ky[i],1.,L)
    espec=eigvals(hammy)
    espec_vec=eigvecs(hammy)
    ind_neg=1
    ind_pos=1
    inf_neg=-1000000000
    inf_pos=100000000
    for pos in 1:length(espec)
        if espec[pos]>infty_neg && espec[pos]<0
            infty_neg=espec[pos]
            ind_neg=pos
        end
        if espec[pos]<infty_pos && espec[pos]>0
            infty_pos=espec[pos]
            ind_pos=pos
        end
    end
    display(espec)
    display(espec_vec)
    #el único eigenfunction que nos importa es el de energías más bajas 
    #for x in espec
    #    display(x)
        #display(scatter!([vals_ky[i]], [x], markercolor = "blue", markerstrokecolor = "blue", markersize = 1,legend=false))
    #end
end