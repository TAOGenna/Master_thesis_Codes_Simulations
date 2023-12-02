




# Function to calculate new expectation values
function calculate_expectation_values(old_values)
        espectrum = eigvals(H)
        eigenstates = eigvecs(H)
        expeta = []
        for i in 1:N #corremos a traves de los eigenenergies cogiendo solos los <0
                if espec[i]<0.
                        #calculo de los expectations values para el mean field - self consistency
                        expeta[1]+=conj(eigenstates[3,i])*eigenstates[3,i]# val1 = < c+1 c_1 > 0.5
                        expeta[2]+=conj(eigenstates[1,i])*eigenstates[1,i]# val2 = < c+d c_d > 0.5
                        expeta[3]+=conj(eigenstates[1,i])*eigenstates[3,i]# val3 = < c+d c_1 > 0.25
                        expeta[4]+=conj(eigenstates[1,i])*eigenstates[4,i]# val4 = < c+d c+1 > 0.14
                end
        end        
    # Replace this with your actual calculation.
    # For the sake of simplicity, this function just adds 1 to each value
    return old_values .+ 1
end


# Convergence check function
function check_convergence(old_values, new_values, tolerance)
    return norm(old_values - new_values) < tolerance
end


V=0
expectation_values=[0.5        ,    0.5   ,   0.25      , 0.14]
#              < c+1 c_1 >  < c+d c_d >   < c+d c_1 >    < c+d c+1 >
while V<2.8
        V += 0.01
        complete(i_values,V)
        tolerance = 1e-6
        max_iterations = 100
        # Self-consistent loop
        for i in 1:max_iterations
            new_expectation_values = calculate_expectation_values(expectation_values)
            if check_convergence(expectation_values, new_expectation_values, tolerance)
                println("Converged after $i iterations")
                expectation_values = new_expectation_values
                break
            end
            expectation_values = new_expectation_values
            if i == max_iterations
                println("Did not converge after maximum number of iterations")
            end
        end        
end