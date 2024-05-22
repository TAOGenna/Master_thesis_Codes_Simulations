using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings
using DataFrames
using IterTools
using Printf
using PrettyTables

global N=40+2 # always an even number +2 (+2 because of the quantum dot)
#-------------Values for the Kitaev Chain-----------
delta=0.2
t=1
mu=0
#-------------Values for the quantum dot------
t_prime=0.2
e_d=0
H_max=2.5
#----------------------------------------------
iterations=100
dt(x,y)= ==(x,y) # like a delta dirac
#---------------Plot Settings-----------------
myplot=plot(xlims=(-H_max,H_max))
title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
xlabel!(L"$e_d$")
ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------

H=zeros(N,N)
legend_added = false
for i in 3:N#     j'
   for j in 3:N# j
        i1=ceil(i/2)
        j1=ceil(j/2)
        if !iseven(i) && !iseven(j)#     <i',c|H|j,c>
             H[i,j]=dt(j1-1,i1)*(-t)+dt(j1+1,i1)*(-t)
        elseif iseven(i) && !iseven(j)#  <i',a|H|j,c>
             H[i,j]=dt(i1,j1-1)*(delta)-dt(i1,j1+1)*(delta)
        elseif !iseven(i) && iseven(j)#  <i',c|H|j,a>
             H[i,j]=dt(i1,j1+1)*(delta)-dt(i1,j1-1)*(delta)
        else#                            <i',a|H|j,a>
             H[i,j]=dt(i1,j1+1)*(t)+dt(i1,j1-1)*(t)
        end
    end
end

H[1,1]+=e_d                  # <0,c|H|0,c>
H[2,2]+=-e_d                 # <0,d|H|0,d>
H[3,1]+=-t_prime             # <1,c|H|0,c>=t_prime
H[4,2]+=t_prime              # <1,d|H|0,d>=-t_prime
H[1,3]+=-t_prime
H[2,4]+=t_prime
H[4,1]=H[1,4]=-delta	     # <1,d|H|0,c>=-delta
H[2,3]=H[3,2]=delta			# <1,c|H|0,d>=delta

# +=conj(eigenstates[1,i])*eigenstates[4,i]# val4 = < c+d c+1 >

#num_elements=200
#vec_temp = [i/(num_elements - 1) for i in range(num_elements)]

#for x in 1:iterations
espec=eigvals(H)
espec_vec=eigvecs(H)
ev=0
for i in 1:N
    if espec[i]<0.
        if abs(espec[i])<10^(-8) #we're using element-wise multiplication and addition
            global ev += conj(espec_vec[1,i]) * espec_vec[4,i]*(1/2)
        else 
            global ev += conj(espec_vec[1,i]) * espec_vec[4,i]
        end
    end
end #at the end of this loop we get the new expectation values, they are supposed to match the initial ones
#end
#result = abs.(ev_to_compare .- ev) .< epsilon_expec
display(ev)# 0.21543740472517486

# V*<d+ c+_1> = delta
# V = 0.2*0.21543740472517486
# V = 0.043087480945034972
