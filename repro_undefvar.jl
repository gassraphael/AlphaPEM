
using AlphaPEM.Fuelcell

try
    # Directly call the function that failed
    params = Fuelcell.eh31_physical_params()
    println("Successfully loaded EH31 physical params")
catch e
    println("Caught expected error: ", e)
    # If it's a stacktrace, print it
    for (exc, bt) in current_exceptions()
        showerror(stdout, exc, bt)
        println()
    end
end
