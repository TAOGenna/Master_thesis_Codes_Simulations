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
# Define the parameters
e_c = 0.0
t = 1.0
H_max = 5.0

# Define the range for e_d
e_d_values = range(-H_max, H_max, length=400)

# Define the functions
lambda_plus(e_d) = 0.5 * (e_d + e_c) + sqrt((e_d - e_c)^2 + 4*t^2)
lambda_minus(e_d) = 0.5 * (e_d + e_c) - sqrt((e_d - e_c)^2 + 4*t^2)

# Prepare data for plotting
lambda_plus_values = [lambda_plus(e_d) for e_d in e_d_values]
lambda_minus_values = [lambda_minus(e_d) for e_d in e_d_values]

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
# Create a scatter plot with the spe
# Plotting
plot(e_d_values, [lambda_plus_values, lambda_minus_values],
    size = (900, 700),
    xlims=(-H_max, H_max),
    ylims=(-8, 8),
    ylabel=L"E",
    xlabel=L"\epsilon_d",
    linecolor=:black,
    legend=false,
    xtickfont=numbers_in_axis_fontsize,
    ytickfont=numbers_in_axis_fontsize,
    linewidth=3,
    xguidefontsize=label_fontsize,
    yguidefontsize=label_fontsize,
    framestyle=:box,
    showaxis=:all,
    minorticks=5,
    tick_direction=:in,
    ticksize=5,
    grid=false,
    margin = 9.2mm
)

# Add a horizontal dashed line at y = 0
hline!([0], linestyle=:dash, color=:blue, linewidth=3)  # Adjust the linewidth as needed
#vline!([0], linestyle=:dash, color=:red, linewidth=3, yrange=(-2, 2))
savefig("molecula.pdf")
