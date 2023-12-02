using CSV
using DataFrames
using Plots
using Measures

H_max = 1.5
ceily = 0.7
width_inches = 3 + 3/8
dpi = 300
width_pixels = round(Int, width_inches * dpi)
fontsize_mm = 1
fontsize_points = fontsize_mm * 2.83465
linewidth_mm = 0.18
linewidth_points = linewidth_mm * 2.83465
label_fontsize = 27

# Read the CSV File for Datos1
#filename_csv1 = "Datos1.csv"
#df1 = DataFrame(CSV.File(filename_csv1))

# Read the CSV File for Datos2
filename_csv2 = "Energy vs e_d, caso V=1 y N=42.csv"
df_combined = DataFrame(CSV.File(filename_csv2))

# Merge the data frames for Datos1 and Datos2
#df_combined = vcat(df1, df2)
x_values = df_combined.col1
y_values = df_combined.col2

x_y_font = 20 
# Create a Scatter Plot for the Combined CSV Files with adjusted marker linewidth
scatter(x_values, y_values,
        size=(width_pixels, width_pixels*(2.0/3.0)),
        xlims=(-H_max, H_max),
        ylims=(-ceily, ceily),
        xtickfont=x_y_font,
        ytickfont=x_y_font,
        linewidth=linewidth_points,
        xguidefontsize=label_fontsize,
        yguidefontsize=label_fontsize,
       # yformatter=:scientific,
        legendfontsize=27,
        framestyle=:box,
        showaxis=:all,
        minorticks=5,
        tick_direction=:in,
        ticksize=5,
        grid=false,
        margin=5mm,
        mlinewidth=0.1,
        markersize=1,
        markercolor = "blue",
        markerstrokecolor = "blue",
        legend=:none)

# Labeling
xlabel!(L"$\mu$")
ylabel!("Energy")
x_pos = -1.2
y_pos = 0.41
in_set_font = 23
separation = 0.073
annotate!( x_pos, y_pos, Plots.text(L"N=20", in_set_font))
annotate!( x_pos, y_pos - separation, Plots.text(L"\Delta=0.2", in_set_font))
annotate!( x_pos, y_pos - separation*2, Plots.text(L"t=1", in_set_font))
annotate!( x_pos, y_pos - separation*3, Plots.text(L"t'=0.2", in_set_font))
annotate!( x_pos, y_pos - separation*4, Plots.text(L"V=1", in_set_font))
annotate!( x_pos, y_pos - separation*5, Plots.text(L"\mu=0", in_set_font))
# Optional: Save the Plot
savefig("v1n20.pdf")
