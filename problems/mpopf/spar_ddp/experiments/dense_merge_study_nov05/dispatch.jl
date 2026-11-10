using LinearAlgebra, InteractiveUtils, Profile
U = UpperTriangular(rand(32,32) + 32I)
B = rand(32,32)
println("DISPATCH ", @which ldiv!(U', B))
println("DISPATCH_UNDERLYING")
println(code_typed(ldiv!, Tuple{typeof(U'),typeof(B)}; optimize=false))
println("BLAS_BACKEND ", BLAS.get_config())
println("BLAS_THREADS ", BLAS.get_num_threads())
println("LOADED_LIBRARIES ", BLAS.lbt_get_config())
