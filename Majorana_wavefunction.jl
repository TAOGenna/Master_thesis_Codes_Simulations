#=
This code consider the case of a Kitaev chain and a
quantum dot attached at one end. There is a repultion considered
between the quantum dot and the first site of the chain. 
This repultion is thought as a Hartree-Fock such that there 
will be shiftings in the occupation of site d and 1, as well
as shifting to the cooper pairs, etc. and all of this changes
are considered in the function ini!(...). The iteration
is because we are solving this problem with a self-consistent
approach. Finally, we plot the results to see how the majorana
fermions are affected. The plot involves seeing how the 
eigenenergies change with respect to the value of e_d, 
the occupation energy of the quantum dot. 
=#
#=
CODE:   What are we plotting?
        SITES VS COEFFICIENTS WEIGHT
=#
using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings
using DataFrames
using IterTools
using Printf
using PrettyTables
using Measures
gr()
global N=100+2 # always an even number +2 (+2 because of the quantum dot)
#-------------Values for the Kitaev Chain-----------
# delta=0.2 | t=1 | mu=0 
delta=0.2
t=1
mu=0
#-------------Values for the quantum dot------
t_prime=0.2
e_d=-1


H_max=2
#----------------------------------------------
dt(x,y)= ==(x,y) # like a delta dirac
string_array2 = [" | d , c >", " | d , a >"," | 1 , c >", " | 1 , a >", " | 2 , c >", " | 2 , a >", " | 3 , c >", " | 3 , a >", " | 4 , c >", " | 4 , a >", " | 5 , c >", " | 5 , a >", " | 6 , c >", " | 6 , a >", " | 7 , c >", " | 7 , a >", " | 8 , c >", " | 8 , a >", " | 9 , c >", " | 9 , a >", " | 10 , c >", " | 10 , a >", " | 11 , c >", " | 11 , a >", " | 12 , c >", " | 12 , a >", " | 13 , c >", " | 13 , a >", " | 14 , c >", " | 14 , a >", " | 15 , c >", " | 15 , a >", " | 16 , c >", " | 16 , a >", " | 17 , c >", " | 17 , a >", " | 18 , c >", " | 18 , a >", " | 19 , c >", " | 19 , a >", " | 20 , c >", " | 20 , a >"]
#---------------Plot Settings-----------------
#myplot=plot(xlims=(-H_max,H_max))
#title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
#xlabel!(L"$e_d$")
#ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------

function ini!(H::Matrix{Float64},val,hopping_shift,ocd,oc1,vv)
    H[4,1] += -vv*val
    H[1,4] += -vv*val
    H[2,3] += vv*val
    H[3,2] += vv*val
    # <1,d|H|0,c>=-delta
    # <1,c|H|0,d>=delta
    #----------------------------------------------------
    H[1,3] += -hopping_shift*vv
    H[3,1] += -hopping_shift*vv
    H[2,4] +=  hopping_shift*vv
    H[4,2] +=  hopping_shift*vv
    #----------------------------------------------------
    
    H[1,1] +=  oc1*vv
    H[2,2] += -oc1*vv
    H[3,3] +=  ocd*vv
    H[4,4] += -ocd*vv
    
end
function generate_array(length::Int,H_max)
    return LinRange(-H_max, H_max, length)
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
#through=generate_array(100,H_max)
through=[e_d]
#display(through)
H[3,1]+=-t_prime             # <1,c|H|0,c>=t_prime
H[4,2]+=t_prime              # <1,d|H|0,d>=-t_prime
H[1,3]+=-t_prime
H[2,4]+=t_prime
#H[4,1]=H[1,4]=-delta        # <1,d|H|0,c>=-delta
#H[2,3]=H[3,2]=delta            # <1,c|H|0,d>=delta
epsilon_expec=1e-4
num_elements = 2000
iterations=50
V=1.0
initial_ev = [0.133222,0.315101,0.531797,0.489767]
legend_added=false
#---------------PLOT SETTINGS---------------------------#
min_energy=0.001
#myplot=plot(xlims=(0,N+20))#,ylims=(-min_energy,min_energy))
#title!(L"$e_d=0.0,\quad 0.5,\quad 1.0$")
#title!(L"$N=40+2,\Delta=0.2,t=-1,t'=-0.2$")
#title!(L"$N=40+2$")
#xlabel!(L"$Site$")
#ylabel!(L"$coefficient$")
#xlabel!(L"$e_d$")
#ylabel!(L"$\langle c^\dagger_d c_d\rangle$")
#-------------------------------------------------------#
x_values = Int[]
y_values = Float64[] 
for values_ed in through # check what is "through" it might be just one number or an array
    H[1,1]=values_ed                  # <0,c|H|0,c>
    H[2,2]=-values_ed                 # <0,d|H|0,d>
    matrix_plot=[]
    for x in 1:iterations
       ini!(H,initial_ev[1],initial_ev[2],(initial_ev[3]-0.5),(initial_ev[4]-0.5),V)
       espec=eigvals(H)
       espec_vec=eigvecs(H)
       #display(H)
       if x==iterations
            matrix_plot=deepcopy(H)
       end
       ini!(H,-initial_ev[1],-initial_ev[2],-(initial_ev[3]-0.5),-(initial_ev[4]-0.5),V)
       ev=[0.0,0.0,0.0,0.0]
       for i in 1:N
           if espec[i]<0.
              if abs(espec[i])<10^(-8) #we're using element-wise multiplication and addition
                  ev[1] += conj(espec_vec[1,i]) * espec_vec[4,i]*(1/2) # < c+d c+1 >
                  ev[2] += conj(espec_vec[1,i]) * espec_vec[3,i]*(1/2) # < c+d c_1 >
                  ev[3] += conj(espec_vec[1,i]) * espec_vec[1,i]*(1/2) # < c+d c_d >
                  ev[4] += conj(espec_vec[3,i]) * espec_vec[3,i]*(1/2) # < c+1 c_1 >
              else 
                  ev[1] += conj(espec_vec[1,i]) * espec_vec[4,i]
                  ev[2] += conj(espec_vec[1,i]) * espec_vec[3,i]
                  ev[3] += conj(espec_vec[1,i]) * espec_vec[1,i] # = < c+d c_d >
                  ev[4] += conj(espec_vec[3,i]) * espec_vec[3,i]
              end
           end
       end #at the end of this loop we get the new expectation values, they are supposed to match the initial ones
       global initial_ev = ev
    end
    #-----------------plotting expectation---------------------------------------------
    #display(scatter!([values_ed], [initial_ev[3]], markercolor = "blue", markerstrokecolor = "blue", markersize = 1.5,label="",markeralpha = 0.5))
    #------------------------------------------------------------------------------------- 
    espec_vec_plot=eigvecs(matrix_plot)
    espec_plot=eigvals(matrix_plot)
    #display(espec_plot)
    #display(espec_plot[21])
    el_=espec_vec_plot[:,51]
    for num in 1:N
        push!(x_values, num)
        push!(y_values, el_[num])
    end
    #=for tmr in 1:N
        if espec_plot[tmr]>-min_energy && espec_plot[tmr]<min_energy && espec_plot[tmr]<0          
            el_elegido=espec_vec_plot[:,N//2]
            #display(el_elegido)
            for num in 1:N
                #push!(x_values, num)
                #push!(y_values, el_elegido[num])            
            end
            break
            
        end
    end=#  
end
#--------------------------------------------------------------------------
# What are we plotting?
# SITES VS COEFFICIENTS WEIGHT
#--------------------------------------------------------------------------
# Width in inches for a single-column format (3 3/8 inches)
width_inches = 3 + 3/8

# Desired DPI for print quality
dpi = 300

# Convert width to pixels
width_pixels = round(Int, width_inches * dpi)

# Determine the font size for the smallest capital letters and numerals (at least 2mm)
fontsize_mm = 1
fontsize_points = fontsize_mm * 2.83465 # Conversion from mm to points

# Determine the curve's linewidth (at least 0.18mm)
linewidth_mm = 0.18
linewidth_points = linewidth_mm * 2.83465 # Conversion from mm to points

# Specific value for the font size of labels and annotations inside the plot
label_fontsize=30
numbers_in_axis_fontsize=20
# Start with an empty scatter plot
scatter([], [],size = (900, 730),
        xtickfont=numbers_in_axis_fontsize, # Font size for x-axis ticks
        ytickfont=numbers_in_axis_fontsize, # Font size for y-axis ticks
        linewidth=linewidth_points,      # Linewidth for the curve
        xguidefontsize=label_fontsize,
        yguidefontsize=label_fontsize,
        legendfontsize=label_fontsize,
        framestyle=:box, # Enclose the plot in a box
        showaxis=:all,
        minorticks=5,
        tick_direction=:in,
        ticksize=5,
        grid=false,
        margin = 5mm,
        markercolor="blue",
        markerstrokecolor="blue",
        markersize=4, label="",
        markeralpha=0.5)

# Add the points
scatter!(x_values, y_values, markercolor="blue", markerstrokecolor="blue", markersize=6, label="", markeralpha=0.5)

# Connect the points with a red line
plot!(x_values, y_values, color="red",legend=false)
xticks!(0:10:100)
ylabel!("Coefficients")
xlabel!("Operator index")
# Annotate with model parameters
x_position=40
y_position=0.35
separa=0.050
annotate!(x_position, y_position, Plots.text(L"N=50", label_fontsize))
annotate!(x_position, y_position - separa, Plots.text(L"\Delta=0.2", label_fontsize))
annotate!(x_position, y_position - 2*separa, Plots.text(L"t'=0.2", label_fontsize))
annotate!(x_position, y_position - 3*separa, Plots.text(L"V=1", label_fontsize))
annotate!(x_position, y_position - 4*separa, Plots.text(L"\epsilon_d=-1", label_fontsize))
annotate!(x_position, y_position - 5*separa, Plots.text(L"\mu=0", label_fontsize))

# Display the plot
display(current())
savefig("wm1.pdf")
#------------------------------------------------------------------------------------------------
#---------------------------DISPLAYING BASIC INFORMATION-----------------------------------------
#=
display("Paramenters value -> delta=$delta | t=$t | V=$V | N=$N")
string_array = ["< c+_d c+_1 >", "< c+_d c_1 >", "< c+_d c_d >", "< c+_1 c_1 >"]
initial_ev=reverse(initial_ev)
string_array=reverse(string_array)
df = DataFrame(strings = string_array, similar_values = initial_ev)
# Display the DataFrame
#display("the bigger values are at the end, this is the last site")
pretty_table(df, header = ["equation", "self-consisten value"], tf = tf_markdown)
=#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
