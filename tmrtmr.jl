using Plots

# Sample data for the plot
x = 1:0.1:10
y = sin.(x)

# Create the plot
p = plot(x, y, linewidth=2, framestyle=:box)

# Customize the tick marks on the x and y axes
xticks_pos = 1:2:9
yticks_pos = 0.5:0.5:1
plot!(xticks=xticks_pos, yticks=yticks_pos)

# Calculate the y-range and x-range
yrange = ylims(p)
xrange = xlims(p)

# Annotate tick marks on the top
for xtick in xticks_pos
    annotate!(p, xtick, yrange[2], Plots.text("", :center, :bottom))
    plot!(p, [xtick, xtick], [yrange[2], yrange[2] + 0.05*(yrange[2]-yrange[1])], line=:solid, color=:black)
end

# Annotate tick marks on the right
for ytick in yticks_pos
    annotate!(p, xrange[2], ytick, Plots.text("", :left, :center))
    plot!(p, [xrange[2], xrange[2] + 0.05*(xrange[2]-xrange[1])], [ytick, ytick], line=:solid, color=:black)
end

# Show the plot
display(p)
