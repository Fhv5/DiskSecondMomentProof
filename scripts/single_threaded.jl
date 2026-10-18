using ArbNumerics
using LinearAlgebra

# Loaded before either CLI entry point by the final evidence launcher.
Threads.nthreads(:default) == 1 && Threads.nthreads(:interactive) == 0 ||
    error("The proof launcher requires --threads=1,0")
Base.JLOptions().nmarkthreads == 1 && Base.JLOptions().nsweepthreads == 0 ||
    error("The proof launcher requires --gcthreads=1")

BLAS.set_num_threads(1)
ccall(ArbNumerics.@libflint(flint_set_num_threads), Cvoid, (Int,), 1)
let flint_threads = ccall(ArbNumerics.@libflint(flint_get_num_threads), Int, ())
    BLAS.get_num_threads() == 1 || error("BLAS did not accept one thread")
    flint_threads == 1 || error("FLINT did not accept one thread")
    println("Thread policy: Julia=1 interactive=0 GC=1 BLAS=1 FLINT=1")
end
