using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings

#global N=40+2 # SIEMPRE NUMERO PAR +2(por el quantum dot)
#-------------Valores para la cadena-----------
delta=0.2
t=1
#-------------Valores para el quantum dot------
t_prime=0.5
#----------------------------------------------
H_max=2.5
I=600# previously "I" was my identity matrix, now I'm using it as a number
vals=collect(0:H_max/I:H_max)
sort!(vals)
dt(x,y)= ==(x,y)
eps=10^(-6)
#plot()
myplot=plot(xlims=(0,H_max),ylims=(-0.1,0.1))

#title!(L"$e_d=0.0,\quad 0.5,\quad 1.0$")
title!(L"$\Delta=0.2,t=1,t'=0.5,\epsilon_d=0.5$")
xlabel!(L"$\mu$")
ylabel!(L"$Eigenvalues$")
#mu=0.01
#H=tmr(mu)
pruebas=[0.5,1.2] # diferentes valores para \e_d
pruebas2=[0.5,1.2]
colores=["blue","magenta","orange","blue","red","black"]
#curvas=[L"$\mu=-0.5$",L"$\mu=0.0$",L"$\mu=0.5$",L"$\mu=1.0$",L"$\mu=2.0$",L"$\mu=2.5$"]
curvas=[L"$\mu=0.0$"]
curvas2=[L"$\epsilon_d=0.5$",L"$\epsilon_d=1.2$",L"$\epsilon_d=1.0$",L"$\epsilon_d=1.5$"]
tama=[1.2,1.2,1.2,1.2]
matrix_size=[22,42,62]
curvas=[L"$L=10$",L"$L=20$",L"$L=30$",L"$L=60$"]

global letsgo=1
for sizes in matrix_size
	e_d=0.5
	N=sizes
	H=zeros(N,N)
	#----completo datos para el quantum dot-------
	H[1,1]=e_d                  # <0,c|H|0,c>
	H[2,2]=-e_d                 # <0,d|H|0,d>
	H[3,1]=H[1,3]=t_prime       # <1,c|H|0,c>=t_prime
	H[4,2]=H[2,4]=-t_prime      # <1,d|H|0,d>=-t_prime
	#---------------------------------------------
	legend_added = false
	for mu in vals
		esc=colores[letsgo]
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
		
		espec=eigvals(H)
		for tmr in espec
			if tmr>-0.1 && tmr <0.1
				if !legend_added
					display(scatter!([mu], [tmr], markercolor = esc, markerstrokecolor = esc, markersize = 0.8,label=curvas[letsgo],markeralpha = 0.5))
					legend_added=true
				else 
					display(scatter!([mu], [tmr], markercolor = esc, markerstrokecolor = esc, markersize = 0.8,label="",markeralpha = 0.5))
				end
			end
		end	
		
	end
	global letsgo=letsgo+1
end