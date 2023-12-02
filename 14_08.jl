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
global N=40 # always an even number +2 (+2 because of the quantum dot)
#-------------Values for the Kitaev Chain-----------
# delta=0.2 | t=1 | mu=0 
delta=1
t=1
#----------------------------------------------
dt(x,y)= ==(x,y) # like a delta dirac
#---------------Plot Settings-----------------
#myplot=plot(xlims=(-H_max,H_max))
#title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
#xlabel!(L"$e_d$")
#ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------
H_max=2.5
upper_limit=3.0
ceily=3

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
axis_fontsize=20
# Create a scatter plot with the specified size, font size for tick numbers, and linewidth
scatter(size = (900, 600),
        xlims=(0,H_max),
        ylims=(-upper_limit,upper_limit),
        #xticks=nothing,
        xtickfont=axis_fontsize, # Font size for x-axis ticks
        ytickfont=axis_fontsize, # Font size for y-axis ticks
        #linewidth=linewidth_points,      # Linewidth for the curve
        xguidefontsize=label_fontsize,
        yguidefontsize=label_fontsize,
       # yformatter=:scientific,
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
xlabel!(L"$\mu$")
#xlabel!("")
ylabel!("Energy")
separa=0.33
annotation_fontsize = label_fontsize
#annotate!( 0.5, 0.000065, Plots.text(L"N=50", annotation_fontsize))
# annotate!( 0.25, 1.3-separa*2, Plots.text(L"\Delta=1", 27))
# annotate!( 0.25, 1.3-separa*3, Plots.text(L"t=1", 27))


function generate_array(length::Int,H_maxx)
    return LinRange(0, H_maxx, length)
end
through=generate_array(200,H_max)
#information_to_csv=Tuple{Float64,Float64}[]
for vals_mu in through
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

    espec_vec=eigvecs(H)
    espectrum=eigvals(H)
    for energies in espectrum
        if energies> - upper_limit && energies<upper_limit
            #push!(information_to_csv,(vals_mu,energies))
            display(scatter!([vals_mu], [energies], markercolor = "black", markerstrokecolor = "black", markersize = 1.3,label="",markeralpha = 1.3))
        end
    end
end
#df = DataFrame(col1 = first.(information_to_csv), col2 = last.(information_to_csv))
#CSV.write("Kita20", df)
savefig("Kita_classic.pdf")
