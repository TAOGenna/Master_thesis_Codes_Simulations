using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings

#-----------------------Cosas Globales--------------------------------------------
#---------------------------------------------------------------------------------
global N=40+2 # SIEMPRE NUMERO PAR +2(por el quantum dot)
global iteraciones=100
#-------------Valores para la cadena-----------
delta=0.2
t=1
#-------------Valores para el quantum dot------
t_prime=-0.2
#----------------------------------------------
H_max=2.5
cuut=0.6
I_menos=100
V=2.8 #The goal is to reach aV=2.8
steps=100
devi=0.02
I=200# previously "I" was my identity matrix, now I'm using it as a number
#vals=collect(-H_max:H_max/I:H_max)
vals=[0.] # probamos un solo valor para el e_d para hacer el mean field
sort!(vals)
dt(x,y)= ==(x,y)
eps=10^(-6)
pruebas=[0.0] # diferentes valores para \mu
#-----------------------------------------------------------------------------------
function iniciales!(H::Matrix{Float64},val1,val2,val3,val4)
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

letsgo=1
for mu in pruebas #we are only intersted in \mu=0
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
                 H[i,j]=dt(i1,j1)*(mu)+dt(i1,j1+1)*(t)+dt(i1,j1-1)*(t)
            end
        end
    end
    #display(H)
    for e_d in vals # we are only intersted in e_d=0 because we observe a peak there
        #----completo datos para el quantum dot-------
        #Esto es incluyendo el termino de hopping dot-1er sitio tmb
        H[1,1]=e_d                  # <0,c|H|0,c>
        H[2,2]=-e_d                 # <0,d|H|0,d>
        H[3,1]=H[1,3]=-t_prime       # <1,c|H|0,c>=t_prime
        H[4,2]=H[2,4]=t_prime      # <1,d|H|0,d>=-t_prime
        #---------------------------------------------
        #espec=eigvals(H) # calculo de eigenenergies
        #display(H)
        #display(espec)
        #eigenstates=eigvecs(H) # calculo de eigenstates
        #meto valores_iniciales en H
        for x in max(0.0,(0.5-devi)):1/steps:min(0.5+devi,1.0)
            for y in max(0.0,(0.5-devi)):1/steps:min(0.5+devi,1.0)
                for z in max(0.0,(0.25-devi)):1/steps:min(0.25+devi,1.0)
                    for k in max(0.0,(0.14-devi)):1/steps:min(0.15+devi,1.0)
                    # do something with x, y, z, k
                        expeta = zeros(4)
                        diff = zeros(4) # aqui guardare la diferencia para ver la convergencia
                        expeta[1]=x
                        expeta[2]=y
                        expeta[3]=z
                        expeta[4]=k
                        temp = expeta
                        iniciales!(H,expeta[1],expeta[2],expeta[3],expeta[4])
                        #for( de 1 hasta 100 ){ iteramos las veces necesarias para que haya convergencia 
                        for ite in 1:iteraciones
                        #   ° calculo los expectation values->valen valores 
                            espec=eigvals(H) # calculo de eigenenergies
                            eigenstates=eigvecs(H) # calculo de los eigenstates para esos nuevos parametros del mean field
                        #   ° resto del hamiltoniano los valores que tenía inicialmente
                            iniciales!(H,-expeta[1],-expeta[2],-expeta[3],-expeta[4])
                            expeta = zeros(4)
                            for i in 1:N #corremos a traves de los eigenenergies cogiendo solos los <0
                                if espec[i]<0. 
                                    #calculo de los expectations values para el mean field - self consistency
                                    expeta[1]+=conj(eigenstates[3,i])*eigenstates[3,i]# val1 = < c+1 c_1 > 0.5
                                    expeta[2]+=conj(eigenstates[1,i])*eigenstates[1,i]# val2 = < c+d c_d > 0.5
                                    expeta[3]+=conj(eigenstates[1,i])*eigenstates[3,i]# val3 = < c+d c_1 > 0.25
                                    expeta[4]+=conj(eigenstates[1,i])*eigenstates[4,i]# val4 = < c+d c+1 > 0.14
                                end
                            end
                            converge=true
                            for tmr in 1:4
                                diff[tmr]=abs(expeta[tmr]-temp[tmr])
                                if diff[tmr]>0.05
                                    converge=false
                                end
                            end
                            if converge==true
                                display("encontroooo")
                                display(temp) #aquí estan los valores originales que probé 
                                display(expeta) # Valores luego de las iteraciones 
                                sleep(20)
                            end
                            if ite != iteraciones # mientras que no sea la ultima iteraciones meto nuevos valores, esto porque si meto la ultima iteracione entonces no hay function que lo borre después
                                #   ° meto los nuevos valores para los valores_iniciales
                                iniciales!(H,expeta[1],expeta[2],expeta[3],expeta[4])
                            end
                        end
                        #--------------- PRINT USEFUL INFORMATION-----------------------------------------------
                        display("-------------NUEVOS-------------------------------")
                        display(expeta)
                        #display("-------------Eigenvalues---Ori--------------------")
                        #display(temp)
                        #display("--------------------------------------------------")
                        display("--------------------------------------------------")
                        display("--------------------------------------------------")
                        chi=[x,y,z,k]
                        chi2=round.(chi,digits=7)
                        display(join(chi2, "|"))
                        #display(diff)
                        #----------------------------------------------------------------------------------------
                    end
                end
            end
        end
    end
    display("tmr")  
end