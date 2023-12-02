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
myplot=plot(xlims=(-H_max,H_max))
title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
xlabel!(L"$e_d$")
ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------

function ini!(H::Matrix{Float64},val1,val2,val3,val4,V)
    H[1,1]+=val1*V         # val1 = < c+1 c_1 > 
    H[2,2]+=-val1*V        
    H[3,3]+=val2*V         # val2 = < c+d c_d >
    H[4,4]+=-val2*V
    #-----------------------------
    H[1,3]+=-val3*V
    H[3,1]+=-val3*V # val3 = < c+d c_1 >
    H[4,1]+=val4*V
    H[1,4]+=val4*V  # val4 = < c+d c+1 >
    H[2,3]+=-val4*V
    H[3,2]+=-val4*V
    H[2,4]+=val3*V
    H[4,2]+=val3*V 
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

H[1,1]=e_d                  # <0,c|H|0,c>
H[2,2]=-e_d                 # <0,d|H|0,d>
H[3,1]=-t_prime             # <1,c|H|0,c>=t_prime
H[4,2]=t_prime              # <1,d|H|0,d>=-t_prime
H[1,3]=-t_prime
H[2,4]=t_prime
 
#--------set up for the search of the four expectarion values-------------
initial_expectation_values=[0.5000000000000027,0.4999999999999968,0.25415779315629816,0.15343176248119264]
V_ = 0.0000001
#while V_ < 2.1
vec=[]
for x in initial_expectation_values
	eps=0.002
	steps=10
	vec_temp = [(x - eps + i * 2 * eps / steps) for i in 0:steps]
	push!(vec,vec_temp)
end
#------------------------------------------------------------
# Generate all combinations of one element from each vector
combinations = Iterators.product(vec...)

# Print each combination
epsilon_expec=1e-3
V_temp=V_
id1=[3,1,1,1]
id2=[3,1,3,4]
ssize = prod(length.(vec))
for c in combinations
	ev=[c[1],c[2],c[3],c[4]] #initial guess 
	ev_to_compare=ev
	iterations=100
	for x in 1:iterations
	   	ini!(H,ev[1],ev[2],ev[3],ev[4],V_temp)
	   	espec=eigvals(H)
		espec_vec=eigvecs(H)
		ini!(H,-ev[1],-ev[2],-ev[3],-ev[4],V_temp) #this is done to clean the matrix
		ev=[0.0,0.0,0.0,0.0]
	    for i in 1:N
		    if espec[i]<0.
		        if abs(espec[i])<10^(-8) #we're using element-wise multiplication and addition
		            ev .+= conj(espec_vec[id1,i]) .* espec_vec[id2,i]*(1/2) #bloque 1
		        else 
		            ev .+= conj(espec_vec[id1,i]) .* espec_vec[id2,i]
		        end
		    end
		end #at the end of this loop we get the new expectation values, they are supposed to match the initial ones
	end # abajo de esto Bloque 2
	result = abs.(ev_to_compare .- ev) .< epsilon_expec
	if all(result)
		display("Found values that coincide with our initial guess for ev")
		display(ev)
		global initial_expectation_values=ev_to_compare
		break
	end
	dife=abs(ev_to_compare[1]-ev[1])
	#----------------------------PRINTABLE---------------------------------------------------------------
	println("left to do: $ssize")
	# Define the string array
	string_array = ["< c+_1 c_1 >", "< c+_d c_d >", "< c+_d c_1 >", "< c+_d c+_1 >"]
	#string_array2 = [" | d , c >", " | d , a >"," | 1 , c >", " | 1 , a >", " | 2 , c >", " | 2 , a >", " | 3 , c >", " | 3 , a >", " | 4 , c >", " | 4 , a >", " | 5 , c >", " | 5 , a >", " | 6 , c >", " | 6 , a >", " | 7 , c >", " | 7 , a >", " | 8 , c >", " | 8 , a >", " | 9 , c >", " | 9 , a >", " | 10 , c >", " | 10 , a >", " | 11 , c >", " | 11 , a >", " | 12 , c >", " | 12 , a >", " | 13 , c >", " | 13 , a >", " | 14 , c >", " | 14 , a >", " | 15 , c >", " | 15 , a >", " | 16 , c >", " | 16 , a >", " | 17 , c >", " | 17 , a >", " | 18 , c >", " | 18 , a >", " | 19 , c >", " | 19 , a >", " | 20 , c >", " | 20 , a >", " | 21 , c >", " | 21 , a >", " | 22 , c >", " | 22 , a >", " | 23 , c >", " | 23 , a >", " | 24 , c >", " | 24 , a >", " | 25 , c >", " | 25 , a >", " | 26 , c >", " | 26 , a >", " | 27 , c >", " | 27 , a >", " | 28 , c >", " | 28 , a >", " | 29 , c >", " | 29 , a >", " | 30 , c >", " | 30 , a >", " | 31 , c >", " | 31 , a >", " | 32 , c >", " | 32 , a >", " | 33 , c >", " | 33 , a >", " | 34 , c >", " | 34 , a >", " | 35 , c >", " | 35 , a >", " | 36 , c >", " | 36 , a >", " | 37 , c >", " | 37 , a >", " | 38 , c >", " | 38 , a >", " | 39 , c >", " | 39 , a >", " | 40 , c >", " | 40 , a >"]
	# Define the expectation value array
	expectation_value = ev_to_compare

	# Define another array with a similar structure to expectation_value
	similar_structure = ev

	#information about the difference between the inicial
	# and the final expectation value 
	_difference=abs.(ev_to_compare .- ev)

	# Create a DataFrame with the string array, the expectation value array, and the similar structure array
	df = DataFrame(strings = string_array, expectation_values = expectation_value, similar_values = similar_structure, other_values = _difference)
	#df2 = DataFrame(strings = string_array2, states = )
	# Display the DataFrame
	pretty_table(df, header = ["equation", "initial expec val", "self-consisten value", "difference"], tf = tf_markdown)
	#-------------------------------------------------------------------------------------------
	global ssize-=1
end
#	display(V_)
#    global V_ += 0.001
#end
#------------------------------------------------------
# bloque 1:
	#display(scatter!([e_d], [suma], markercolor = "magenta",label="", markerstrokecolor = "magenta", markersize = 1.0,markeralpha = 0.5))

	#ev[1]+=conj(espec[3,i])*espec[3,i] #= < c+1 c_1 > 0.5
	#ev[1]+=conj(espec[1,i])*espec[1,i] #= < c+d c_d > 0.5
	#ev[1]+=conj(espec[1,i])*espec[3,i] #= < c+d c_1 > 0.25+epsilon(it could be 0.05)
	#ev[1]+=conj(espec[1,i])*espec[4,i] #= < c+d c+1 > 0.15+epsilon(it could be 0.05)

	#ev[1]+=conj(espec_vec[3,i])*espec_vec[3,i] # = < c+1 c_1 > 0.5  --------------------------------------
	#ev[2]+=conj(espec_vec[1,i])*espec_vec[1,i] # = < c+d c_d > 0.5
	#ev[3]+=conj(espec_vec[1,i])*espec_vec[3,i] # = < c+d c_1 > 0.25+epsilon(it could be 0.05)
	#ev[4]+=conj(espec_vec[1,i])*espec_vec[4,i] # = < c+d c+1 > 0.15+epsilon(it could be 0.05)
# valores iniciales para un V=0
	# 0.5000000000000027
	# 0.4999999999999968
	# 0.25415779315629816
	# 0.15343176248119264
# Bloque 2: Cosas que podrian dar info sobre lo que sale de los expectation values en la marcha 
	#println("ev_modified -> ", join(["| $(@sprintf("%.7f", x)) " for x in ev], "|"))
	#println("ev_original -> ", join(["| $(@sprintf("%.7f", x)) " for x in ev_to_compare], "|"))
	#println("--------------------------left to do $ssize ------------------------------------")
	#global ssize-=1
	#display(ev)
	#display(ev_to_compare)
	#display("///////////////////////////////////////")
#----------------------------------------------------------
# computing the expectation values 
ev[1]+=conj(espec_vec[3,i])*espec_vec[3,i] # = < c+1 c_1 >
ev[2]+=conj(espec_vec[1,i])*espec_vec[1,i] # = < c+d c_d > 
ev[3]+=conj(espec_vec[1,i])*espec_vec[3,i] # = < c+d c_1 >
ev[4]+=conj(espec_vec[1,i])*espec_vec[4,i] # = < c+d c+1 >
