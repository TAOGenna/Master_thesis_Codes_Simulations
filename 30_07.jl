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
dt(x,y)= ==(x,y) # like a delta dirac
#---------------Plot Settings-----------------
#myplot=plot(xlims=(-H_max,H_max))
#title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
#xlabel!(L"$e_d$")
#ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------

function ini!(H::Matrix{Float64},val,vv)
    H[4,1]=H[1,4]=-vv*val         # <1,d|H|0,c>=-delta
    H[2,3]=H[3,2]=vv*val           # <1,c|H|0,d>=delta 
end
function create_array(number::Float64, length::Int, epsilon::Float64)
    # Determine the step size for the given range and length
    step = 2.0 * epsilon / (length - 1)
    
    # Generate the array
    array = [number - epsilon + i*step for i in 0:(length - 1)]
    
    return array
end

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
#H[4,1]=H[1,4]=-delta	     # <1,d|H|0,c>=-delta
#H[2,3]=H[3,2]=delta			# <1,c|H|0,d>=delta

# +=conj(eigenstates[1,i])*eigenstates[4,i]# val4 = < c+d c+1 >
epsilon_expec=1e-4
num_elements = 2000
iterations=200
#vec_temp = [i/(num_elements - 1) for i in 0:(num_elements - 1)]
#display(vec_temp)
vec_temp = create_array(0.21543740472517486,2000,1e-3)
for values in vec_temp
    values_temp=values
    for x in 1:iterations
        V=0.92834389763993045514561204877889
        ini!(H,values_temp,V)
        espec=eigvals(H)
        espec_vec=eigvecs(H)
        ev=0
        for i in 1:N
            if espec[i]<0.
                if abs(espec[i])<10^(-8) #we're using element-wise multiplication and addition
                    ev += conj(espec_vec[1,i]) * espec_vec[4,i]*(1/2)
                else 
                    ev += conj(espec_vec[1,i]) * espec_vec[4,i]
                end
            end
        end #at the end of this loop we get the new expectation values, they are supposed to match the initial ones
        values_temp=ev
    end
    if abs(values-values_temp)<epsilon_expec
        display("Found values that coincide with our initial guess for ev")
        display("Found values that coincide with our initial guess for ev")
        display(values)
        display(values_temp)
        display(abs(values-values_temp))
        break
    end
end

#result = abs.(ev_to_compare .- ev) .< epsilon_expec
#display(ev)# 0.21543740472517486
#V=0.92834389763993045514561204877889


"Found values that coincide with our initial guess for ev"
#guessed value = 0.21533785495028743
#self consitence value = 0.21543740472517156
#absolute difference = 9.954977488413341e-5









#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#"Found values that coincide with our initial guess for ev"
#0.15577889447236182

#"Found values that coincide with our initial guess for ev" parameters -> num_elements=700, iterations=150, epsilon_expec=1e-3
#0.15515515515515516
#0.15554429665146574
#0.0003891414963105855

"Found values that coincide with our initial guess for ev"
#epsilon_expec=1e-4
#num_elements = 2000
#iterations=200
#0.15557778889444723
#0.15554429665146574
#3.3492242981492115e-5

#"Found values that coincide with our initial guess for ev" -> using vec_temp = create_array(0.15557778889444723,2000,1e-4)
#0.15553431715857932
#0.15554429665146574
#9.979492886419417e-6

#"Found values that coincide with our initial guess for ev" -> vec_temp = create_array(0.15553431715857932,2000,1e-4)
#epsilon_expec=1e-6
#0.15554337168584295
#0.15554429665146574
#9.24965622789431e-7

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------