using Plots

open("out.txt", "w") do file
    for i in 1:100
        println(file, "This is line number $i")
    end
end