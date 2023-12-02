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

df = DataFrame(CSV.File("E vs e_d V=1 y N=102 y mu=0.5.csv"))
display(df)