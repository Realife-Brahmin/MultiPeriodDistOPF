function count_ways_dp(T::Int)
    # Initialize a (T + 1) x 3 array of zeros
    f = zeros(Int, T + 1, 3)
    f[0+1, 0+1] = 1  # Base case: f[0][0] = 1

    # Print initial state
    println("Initial state:")
    println("f[0][0] = 1")

    # Iterate over positions from 1 to T
    for n in 1:T
        # Initialize f[n][0], f[n][1], f[n][2] to 0
        f_n0 = 0
        f_n1 = 0
        f_n2 = 0

        # Equation 1: f[n][0] = f[n-3][0] + f[n-3][1] + f[n-3][2]
        if n - 3 >= 0
            sum_prev = f[n-3+1, 0+1] + f[n-3+1, 1+1] + f[n-3+1, 2+1]
            f_n0 = sum_prev
            println("f[$n][0] = f[$n-3][0] + f[$n-3][1] + f[$n-3][2] = $sum_prev")
        else
            println("f[$n][0] remains 0 (n - 3 < 0)")
        end

        # Equation 2: f[n][1] = f[n-1][0]
        if n - 1 >= 0
            f_n1 = f[n-1+1, 0+1]
            println("f[$n][1] = f[$n-1][0] = f[$n-1][0] = $(f_n1)")
        else
            println("f[$n][1] remains 0 (n - 1 < 0)")
        end

        # Equation 3: f[n][2] = f[n-2][0]
        if n - 2 >= 0
            f_n2 = f[n-2+1, 0+1]
            println("f[$n][2] = f[$n-2][0] = f[$n-2][0] = $(f_n2)")
        else
            println("f[$n][2] remains 0 (n - 2 < 0)")
        end

        # Update the DP table
        f[n+1, 0+1] = f_n0
        f[n+1, 1+1] = f_n1
        f[n+1, 2+1] = f_n2
    end

    # Calculate the total number of ways
    total_ways = f[T+1, 0+1] + f[T+1, 1+1] + f[T+1, 2+1]
    println("Total ways for T = $T:")
    println("f[$T][0] = $(f[T + 1, 0 + 1])")
    println("f[$T][1] = $(f[T + 1, 1 + 1])")
    println("f[$T][2] = $(f[T + 1, 2 + 1])")
    println("Total ways = f[$T][0] + f[$T][1] + f[$T][2] = $total_ways")

    return total_ways
end

# Example usage:
T = 24
println("Computing the total number of ways for T = $T:")
total = count_ways_dp(T)
println("The total number of ways for T = $T is $total.")
