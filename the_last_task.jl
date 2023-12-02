# PLOTTING KITAEV CHAIN WITHOUT THE QUANTUM DOT

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
global N=100 # always an even number +2 (+2 because of the quantum dot)
#-------------Values for the Kitaev Chain-----------
# delta=0.2 | t=1 | mu=0 
delta=0.2
t=1
#----------------------------------------------
dt(x,y)= ==(x,y) # like a delta dirac


function generate_array(length::Int,H_maxx)
    return LinRange(0, H_maxx, length)
end
H_max=2
through=generate_array(2,H_max)
#information_to_csv=Tuple{Float64,Float64}[]
x_values = Float64[]
y_values = Float64[] 
print("N=",N/2,"\n")
print("--------------------------------------------------------------------\n")
for vals_mu in [0.0,0.5,1.0,1.5,1.9,2.1,2.5]
    H=zeros(N,N)
    legend_added = false
    for i in 1:N#     j'
       for j in 1:N# j
            i1=ceil(i/2)
            j1=ceil(j/2)
            if !iseven(i) && !iseven(j)#     <i',c|H|j,c>
                 H[i,j]=dt(j1,i1)*(-vals_mu)+dt(j1-1,i1)*(-t)+dt(j1+1,i1)*(-t)
            elseif iseven(i) && !iseven(j)#  <i',a|H|j,c>
                 H[i,j]=dt(i1,j1-1)*(delta)-dt(i1,j1+1)*(delta)
            elseif !iseven(i) && iseven(j)#  <i',c|H|j,a>
                 H[i,j]=dt(i1,j1+1)*(delta)-dt(i1,j1-1)*(delta)
            else#                            <i',a|H|j,a>
                 H[i,j]=dt(i1,j1)*(vals_mu)+dt(i1,j1+1)*(t)+dt(i1,j1-1)*(t)
            end
        end
    end
    eigenvalues, eigenvectors = eigen(H)
    #display(eigenvalues)
    # Find the index of the smallest eigenvalue
    ground_state_index = argmin(abs.(eigenvalues))
    #display(eigenvectors[:,20])
    #display(eigenvectors[:,21])
    # Get the ground state eigenvector
    ground_state = eigenvectors[:, ground_state_index]
    el_=ground_state
    quan=0.0
    quan_single_term=0.0
    for num in 1:N
          if mod(num,2)==1
               quan += (el_[num]*conj(el_[num])-el_[num+1]*conj(el_[num+1]))
               quan_single_term += el_[num]*conj(el_[num])
          end
    end
    print("mu: ",vals_mu, " -> sum: ", quan, "\n")
    print("sum for a single term: ", quan_single_term, "\n")
    print("-------------------------------------------------------------------- \n")
    #push!(x_values, vals_mu)
    #push!(y_values, quan)
end


#=
####################### COSAS PLOT ####################333
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
label_fontsize=25
# Start with an empty scatter plot
scatter([], [],size = (width_pixels, width_pixels),
        xtickfont=label_fontsize-5, # Font size for x-axis ticks
        ytickfont=label_fontsize-5, # Font size for y-axis ticks
        linewidth=linewidth_points,      # Linewidth for the curve
        xguidefontsize=label_fontsize,
        yguidefontsize=label_fontsize,
        legendfontsize=20,
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
#plot!(x_values, y_values, color="red",legend=false)
#xticks!(0:10:100)
ylabel!(L"\sum_i (|D_{i\nu}|^2-|E_{i\nu}|^2)")
xlabel!(L"\mu")
=#