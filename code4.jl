using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings

global N=40+2 # SIEMPRE NUMERO PAR +2(por el quantum dot)
zer1=zeros(N,N)
zer2=zeros(N,N)
zer3=zeros(N,N)
#-------------Valores para la cadena-----------
delta=0.2
t=1
#-------------Valores para el quantum dot------
t_prime=0.2
#----------------------------------------------
H_max=2.5
cuut=0.6
I_menos=100
I=40# previously "I" was my identity matrix, now I'm using it as a number
vals=collect(-H_max:H_max/I:H_max)
sort!(vals)
dt(x,y)= ==(x,y)
eps=10^(-6)
#plot()
myplot=plot(xlims=(-H_max,H_max),ylims=(-0.1,0.1))
#title!(L"$e_d=0.0,\quad 0.5,\quad 1.0$")
title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0$")
xlabel!(L"$\epsilon_d$")
ylabel!(L"$Eigenvalues$")
#mu=0.01
#H=tmr(mu)
pruebas=[0.0] # diferentes valores para \mu
pruebas2=[1.5,2.0,2.5]
colores=["blue","magenta","orange","blue","red","black"]
#curvas=[L"$\mu=-0.5$",L"$\mu=0.0$",L"$\mu=0.5$",L"$\mu=1.0$",L"$\mu=2.0$",L"$\mu=2.5$"]
curvas=[L"$\mu=0.0$"]
curvas2=[L"$\mu=1.5$",L"$\mu=2.0$",L"$\mu=2.5$"]
tama=[1.2]

#for e_d in vals
#	esc=colores[rand(1:7)]
global letsgo=1
for mu in pruebas
	esc=colores[letsgo]
	H=zeros(N,N)
	for i in 3:N#     j'
		for j in 3:N# j
			i1=ceil(i/2)
			j1=ceil(j/2)
			if !iseven(i) && !iseven(j)#     <j',c|H|j,c>
				H[i,j]=dt(j1,i1)*(-mu)+dt(j1-1,i1)*(t)+dt(j1+1,i1)*(t)
			elseif iseven(i) && !iseven(j)#  <j',a|H|j,c>
				H[i,j]=dt(i1,j1-1)*(delta)+dt(i1,j1+1)*(-delta)
			elseif !iseven(i) && iseven(j)#  <j',c|H|j,a>
				H[i,j]=dt(i1,j1+1)*(delta)+dt(i1,j1-1)*(-delta)
			else#                            <j',a|H|j,a>
			H[i,j]=dt(i1,j1)*(mu)+dt(i1,j1+1)*(-t)+dt(i1,j1-1)*(-t)
			end
		end
	end
	legend_added = false
	for e_d in vals
		#----completo datos para el quantum dot-------
		H[1,1]=e_d                  # <0,c|H|0,c>
		H[2,2]=-e_d                 # <0,d|H|0,d>
		H[3,1]=H[1,3]=t_prime       # <1,c|H|0,c>=t_prime
		H[4,2]=H[2,4]=-t_prime      # <1,d|H|0,d>=-t_prime
		#---------------------------------------------
		espec=eigvals(H)
		eigenstates=eigvecs(H)
		for tmr in espec
			display(scatter!([e_d], [tmr], markercolor = esc, markerstrokecolor = esc, markersize = tama[letsgo],label="",markeralpha = 0.5))
		end
		#=
		suma=0
		for i in 1:N
			if espec[i]<0.
				if abs(espec[i])<10^(-8)
					suma=suma+conj(eigenstates[1,i])*eigenstates[1,i]*(1/2)
				else 
					suma=suma+conj(eigenstates[1,i])*eigenstates[1,i]# < c+d c_d >
				end
			end
		end
		
		if !legend_added
			display(scatter!([e_d], [suma], markercolor = esc, markerstrokecolor = esc, markersize = tama[letsgo],label=curvas[letsgo],markeralpha = 0.5))
			legend_added=true
		else 
			display(scatter!([e_d], [suma], markercolor = esc, markerstrokecolor = esc, markersize = tama[letsgo],label="",markeralpha = 0.5))
		end
		=#
	end
	global letsgo=letsgo+1 		
end
#end


#=

for ind in 1:size(pruebas)[1]
	local H=zeros(N,N)
	local mu=pruebas[ind]
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
	display("valor: ",mu)
	display(H)
	for e_d in vals
		#----completo datos para el quantum dot-------
		H[1,1]=e_d                  # <0,c|H|0,c>
		H[2,2]=-e_d                 # <0,d|H|0,d>
		H[3,1]=H[1,3]=t_prime       # <1,c|H|0,c>=t_prime
		H[4,2]=H[2,4]=-t_prime      # <1,d|H|0,d>=-t_prime
		#---------------------------------------------
		espec=eigvals(H)
		eigenstates=eigvecs(H)
		#------- now we only care for <c+_d c_d> so we only focus on terms B_v1, where v
		# stands for the corresponding eigenvalue -----
		suma=0
		for i in 1:N
			suma=eigenstates[1,i]^2
		end 
		display(scatter!([e_d], [suma], color = colores[ind], label = "", markersize = 1))
	end
end
x = range(0, 2.5, length=100)
y1 = sin.(x)
y2 = cos.(x)
display(scatter!(x, [y1 y2]))
=#