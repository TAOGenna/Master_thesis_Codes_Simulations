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
e_d=-0.01
H_max=2.5
#----------------------------------------------
dt(x,y)= ==(x,y) # like a delta dirac
string_array2 = [" | d , c >", " | d , a >"," | 1 , c >", " | 1 , a >", " | 2 , c >", " | 2 , a >", " | 3 , c >", " | 3 , a >", " | 4 , c >", " | 4 , a >", " | 5 , c >", " | 5 , a >", " | 6 , c >", " | 6 , a >", " | 7 , c >", " | 7 , a >", " | 8 , c >", " | 8 , a >", " | 9 , c >", " | 9 , a >", " | 10 , c >", " | 10 , a >", " | 11 , c >", " | 11 , a >", " | 12 , c >", " | 12 , a >", " | 13 , c >", " | 13 , a >", " | 14 , c >", " | 14 , a >", " | 15 , c >", " | 15 , a >", " | 16 , c >", " | 16 , a >", " | 17 , c >", " | 17 , a >", " | 18 , c >", " | 18 , a >", " | 19 , c >", " | 19 , a >", " | 20 , c >", " | 20 , a >"]
#---------------Plot Settings-----------------
#myplot=plot(xlims=(-H_max,H_max))
#title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
#xlabel!(L"$e_d$")
#ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------

function ini!(H::Matrix{Float64},val,vv,hopping_shift)
    H[4,1]=H[1,4]=-vv*val         # <1,d|H|0,c>=-delta
    H[2,3]=H[3,2]=vv*val           # <1,c|H|0,d>=delta
    #----------------------------------------------------
    H[1,3] += -hopping_shift*vv
    H[3,1] += -hopping_shift*vv
    H[2,4] +=  hopping_shift*vv
    H[4,2] +=  hopping_shift*vv
end
function create_array(number::Float64, length::Int, epsilon::Float64)
    # Determine the step size for the given range and length
    step = 2.0 * epsilon / (length - 1)
    
    # Generate the array
    array = [number - epsilon + i*step for i in 0:(length - 1)]
    
    return array
end
function separator()
	display("-------------------------------------------------------------------------")
	display("-------------------------------------------------------------------------")
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
iterations=40
#display(eigvecs(H))
#tmr=eigvecs(H)
#display(tmr[:,1])
#V=0.92834389763993045514561204877889
V=1.0
initial_ev = [0.21533785495028743,0.0,0.0,0.0,0.0]
initial_temp=initial_ev
minimun=-100000000.0
majorana_eigenstate=[]
for x in 1:iterations
   ini!(H,initial_ev[1],V,initial_ev[2])
   espec=eigvals(H)
   espec_vec=eigvecs(H)
   if x==iterations
   		display("Eigenenergies:")
   		display(espec)
   		separator()
   		#second_excitation=espec[3] # the second excitation it is not in [3] because it is not ordered, you can see
   		#display("energy of the second excited state = $second_excitation")# that by displaying the espectrum
   		#display(espec_vec[:,3])	
   		#separator()
	   	display("Now we print states whose site d has at least 0.2 as its eigenenergy:")
	   	for mk in 1:N
	   		energia=espec[mk]
	   		if abs(espec_vec[1,mk])>=0.2
	   			display("Eigenenergy = $energia")
	   			display("Eigenstate: ")
	   			df2 = DataFrame(strings = string_array2, states = espec_vec[:,mk])
	   			pretty_table(df2, header = ["eigenvector", "value"], tf = tf_markdown)
	   			#display(espec_vec[:,mk])
	   			separator()
	   		end
	   	end
   end
   ini!(H,initial_ev[1],V,-initial_ev[2])
   ev=[0.0,0.0,0.0,0.0,0.0]
   for i in 1:N
       if espec[i]<0.
       	if minimun<espec[i]
       		global minimun=max(minimun,espec[i])
       		global majorana_eigenstate=espec_vec[:,i]
       	end
          if abs(espec[i])<10^(-8) #we're using element-wise multiplication and addition
              ev[1] += conj(espec_vec[1,i]) * espec_vec[4,i]*(1/2) # < c+d c+1 >
              ev[2] += conj(espec_vec[1,i]) * espec_vec[3,i]*(1/2) # < c+d c_1 >
              ev[3] += conj(espec_vec[1,i]) * espec_vec[1,i]*(1/2) # < c+d c_d >
              ev[4] += conj(espec_vec[3,i]) * espec_vec[3,i]*(1/2) # < c+1 c_1 >
              ev[5] += conj(espec_vec[3,i]) * espec_vec[5,i]*(1/2) # < c+1 c_2 >
          else 
              ev[1] += conj(espec_vec[1,i]) * espec_vec[4,i]
              ev[2] += conj(espec_vec[1,i]) * espec_vec[3,i]
              ev[3] += conj(espec_vec[1,i]) * espec_vec[1,i] # = < c+d c_d >
              ev[4] += conj(espec_vec[3,i]) * espec_vec[3,i]
              ev[5] += conj(espec_vec[3,i]) * espec_vec[5,i]
          end
       end
   end #at the end of this loop we get the new expectation values, they are supposed to match the initial ones
   global initial_ev = ev
end
display("when V = $V")
minimun=abs(minimun)
display("Majorana fermion's energy = $minimun")
#df2 = DataFrame(strings = string_array2, states = majorana_eigenstate)	
#pretty_table(df2, header = ["eigenvector", "value"], tf = tf_markdown)
string_array = ["< c+_d c+_1 >", "< c+_d c_1 >", "< c+_d c_d >", "< c+_1 c_1 >", "< c+_1 c_2 >"]
df = DataFrame(strings = string_array, expectation_values = initial_temp, similar_values = initial_ev)
# Display the DataFrame
display("the bigger values are at the end, this is the last site")
pretty_table(df, header = ["equation", "initial expec val", "self-consisten value"], tf = tf_markdown)
	
