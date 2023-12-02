using CSV
using DataFrames
using Plots
using LaTeXStrings
using Measures

# Plotting------------------------------------
label_fontsize=30
axis_fontsize=20
upper_limit=7
H_max=3
scatter(size=(900, 600),
        xlims=(-H_max, H_max),
        ylims=(-upper_limit, upper_limit),
        xtickfont=axis_fontsize,
        ytickfont=axis_fontsize,
        xguidefontsize=label_fontsize,
        yguidefontsize=label_fontsize,
        #yformatter=:scientific,
        legendfontsize=label_fontsize,
        framestyle=:box,
        showaxis=:all,
        minorticks=5,
        tick_direction=:in,
        ticksize=5,
        grid=false,
        margin=7mm,
        markersize=1.6, # Make the dots smaller
        color="black" # Set the color of the dots
)
#---------------------------------------------------

#H_max=3.5
t=1
delta=1
mu=2.75*t
function generate_array(length::Int,H_maxx)
    return LinRange(-H_maxx, H_maxx, length)
end
through=generate_array(300,H_max)
for k in through
    foo=sqrt((2*t*cos(k)+mu)^2+4*(delta)^2*(sin(k))^2)
    display(scatter!([k], [ foo ], markercolor = "black", markerstrokecolor = "black", markersize = 2,label="",markeralpha = 2))
    display(scatter!([k], [-foo ], markercolor = "black", markerstrokecolor = "black", markersize = 2,label="",markeralpha = 2))
end

ylabel!(L"E")
#xlabel!("asd")
xlabel!(L"$k$")
#ot=5/100000
annotate!(0, 2.5, Plots.text(L"\mu=2.75t", label_fontsize))
#annotate!(0, 4  , Plots.text("Topológico", label_fontsize-10))
annotate!(0, 1.5  , Plots.text("Trivial", label_fontsize-10))
# Display the plot
filename_pdf = "kita_bound6.pdf"
savefig(filename_pdf)
#display(plot)