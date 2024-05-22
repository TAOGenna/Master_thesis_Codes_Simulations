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

CODE: PLOTS ENERGY VS e_d | Change sites in variable N, range x in array "through", range y in variable "min_energy"

=#
using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings
using DataFrames
using IterTools
using Printf
using PrettyTables
using CSV
using Measures
gr()
global N=10+2 # always an even number +2 (+2 because of the quantum dot)
#-------------Values for the Kitaev Chain-----------
# delta=0.2 | t=1 | mu=0 
delta=0.2
t=1
mu=0
#-------------Values for the quantum dot------
t_prime=0.2
e_d=-0.01
H_max=1.5
V=0
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
    #display(vv)
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
through=generate_array(300,H_max)
#through=[0.0]
#display(through)
H[3,1]+=-t_prime             # <1,c|H|0,c>=t_prime
H[4,2]+=t_prime              # <1,d|H|0,d>=-t_prime
H[1,3]+=-t_prime
H[2,4]+=t_prime
#H[4,1]=H[1,4]=-delta        # <1,d|H|0,c>=-delta
#H[2,3]=H[3,2]=delta            # <1,c|H|0,d>=delta
epsilon_expec=1e-4
num_elements = 2000
iterations=40
initial_ev = [0.133222,0.315101,0.531797,0.489767]
legend_added=false
#---------------PLOT SETTINGS---------------------------#
#min_energy=1# aumentar esto para el caso de N=40
#myplot=plot(xlims=(-H_max,H_max),ylims=(-min_energy,min_energy))
#title!(L"$e_d=0.0,\quad 0.5,\quad 1.0$")
#title!(L"$N=40+2,\Delta=0.2,t=1,t'=0.2, V=0$")
#xlabel!(L"$e_d$")
#ylabel!(L"$Eigenenergies$")
#-------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
#---------------PLOT SETTINGS--------------------------------------------------------------------------------#
min_energy=0.6# aumentar esto para el caso de N=40
min_energy_=0.05
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
# Create a scatter plot with the specified size, font size for tick numbers, and linewidth
p=scatter(size = (900, 730),
        xlims=(-H_max,H_max),
        ylims=(-min_energy,min_energy),
        xticks=nothing,
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
        margin = 7mm
        #mirror=true
)

xlabel!("")
#xlabel!(L"\epsilon_d")
ylabel!("Energy")

# Annotate with model parameters
x_pos = -1
y_pos = 0.44
in_set_font = 40
separation = 0.076
#annotate!( x_pos, y_pos, Plots.text(L"N=20", in_set_font))
#annotate!( x_pos, y_pos - separation, Plots.text(L"\Delta=0.2", in_set_font))
#annotate!( x_pos, y_pos - separation*2, Plots.text(L"t=1", in_set_font))
#annotate!( x_pos, y_pos - separation*3, Plots.text(L"t'=0.2", in_set_font))
annotate!( x_pos, y_pos - separation, Plots.text(L"V=0", label_fontsize))
#annotate!( x_pos, y_pos - separation*5, Plots.text(L"\mu=0", in_set_font))
contour_thickness = 0.5
#---------------------------------------END OF PLOT SETTINGS-------------------------------------------------#
#------------------------------------------------------------------------------------------------------------#
information_to_csv=Tuple{Float64,Float64}[]
#-------------------------------------------------------------------------------------------------------------
for values_ed in through
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
            #display(matrix_plot[1:4,1:4])
        
            espec_plot=eigvals(matrix_plot)
            for tmr in espec_plot
                if tmr>-min_energy && tmr<min_energy
                    #push!(information_to_csv,(values_ed,tmr))
                    if !legend_added
                        display(scatter!([values_ed], [tmr], markercolor = "black", markerstrokecolor = "black", markersize = 1.9,label="",markeralpha = 1.9))
                        global legend_added=true
                    else 
                        display(scatter!([values_ed], [tmr], markercolor = "black", markerstrokecolor = "black", markersize = 1.9,label="",markeralpha = 1.9))
                    end
                end
            end  
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
    #-----------------plotting eigenenergies---------------------------------------------
    #display(matrix_plot[1:4,1:4])
    
    #-------------------------------------------------------------------------------------   
end
#df = DataFrame(col1 = first.(information_to_csv), col2 = last.(information_to_csv))
#filename_csv = "Energy vs e_d, caso V=$(V) y N=$(N).csv"
filename_pdf = "diaV0.pdf"

#CSV.write(filename_csv, df)
savefig(filename_pdf)
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
