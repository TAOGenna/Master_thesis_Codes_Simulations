using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings

function f(x,k)
	eps2=0.004
	if x>=k && x<=k+eps2
		return true
	else
		return false
	end
end

sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
N=20+2 # SIEMPRE NUMERO PAR +2(por el quantum dot)
zer1=zeros(N,N)
zer2=zeros(N,N)
zer3=zeros(N,N)
#-------------Valores para la cadena-----------
delta=0.2
t=1
#-------------Valores para el quantum dot------
e_d=0
t_prime=0.5
#----------------------------------------------
H_max=2.5
I=200# previously "I" was my identity matrix, now I'm using it as a number
vals=collect(0:H_max/(I-1):H_max)#posibles valores que podría tomar \mu
dt(x,y)= ==(x,y)
#display(dt(3,2))
eps=10^(-1)
#myplot=plot(xlims=(0.,2.5), ylims=(-eps,eps)	) #creates a plot doing a close up at the given limits 
plot() # creates an empty plot
#myplot=plot(ylims=(0,1.2)	)
title!("Occupation number for the Quantum Dot, 10 sites.")
xlabel!(L"\mu")
ylabel!(L"$\langle c^\dagger_d c_1\rangle$")
for mu in vals
	H=zeros(N,N)
	for i in 3:N#     j'
		for j in 3:N# j
			i1=ceil(i/2)
			j1=ceil(j/2)
			if !iseven(i) && !iseven(j)#     <j',c|H|j,c>
				H[i,j]=dt(j1,i1)*(-mu)+dt(j1-1,i1)*(-t)+dt(j1+1,i1)*(-t)
			elseif iseven(i) && !iseven(j)#  <j',a|H|j,c>
				H[i,j]=dt(i1,j1-1)*(delta)+dt(i1,j1+1)*(-delta)
			elseif !iseven(i) && iseven(j)#  <j',c|H|j,a>
				H[i,j]=dt(i1,j1+1)*(delta)+dt(i1,j1-1)*(-delta)
			else#                            <j',a|H|j,a>
				H[i,j]=dt(i1,j1)*(mu)+dt(i1,j1+1)*t+dt(i1,j1-1)*t
			end
		end
	end
	#----completo datos para el quantum dot-------
	H[1,1]=e_d                  # <0,c|H|0,c>
	H[2,2]=-e_d                 # <0,d|H|0,d>
	H[3,1]=H[1,3]=t_prime       # <1,c|H|0,c>=t_prime
	H[4,2]=H[2,4]=-t_prime      # <1,d|H|0,d>=-t_prime
	#---------------------------------------------
	espec=eigvals(H)
	espec_vec=eigvecs(H)
	if mu==vals[1]
		display(espec_vec[1:N,2])
		display(espec_vec)
	end
	#------- now we only care for <c+_d c_d> so we only focus on terms B_v1, where v
	# stands for the corresponding eigenvalue -----
	suma=0
	for i in 1:N
		if mu==vals[1]
			display(espec_vec[1,i])
		end
		suma=espec_vec[1,i]*espec_vec[3,i]
	end 
	display(scatter!([mu], [suma], color = "blue", label = "", markersize = 1))


	#---------- cuenta de los zero modes cuando mu=0---------------------
	#=
	count_zero_modes=0
	if mu==vals[1]
		for asd in espec
			if asd<eps && asd>-eps
				count_zero_modes=count_zero_modes+1
			end
		end
		display("cuenta de los zero modes en mu=0:")
		display(count_zero_modes)
	end
	=#
	#----------------------------------------------------------------------
	
	#----------------------- bolds the lowest positive eigenvalue for a given \mu-----------------------
	#=
	lw_value=100
	lw_index=1
	index=1
	for asd in espec
		if asd>=0 && asd<lw_value# lowest positive eigenvalue for a given \mu
			lw_value=asd
			lw_index=index
		end
		display(scatter!([mu], [asd], color = "green", label = "", markersize = 0.3))
		index+=1
	end
	display(scatter!([mu], [lw_value], color = "blue", label = "", markersize = 1))
	=#
	#---------------------------------------------------------------------------------------------------	
end
#savefig(myplot,"file.png")
#savefig("com_delta23.png")