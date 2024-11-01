module openDSSValidator

using OpenDSSDirect

export validate_opf_against_opendss

function validate_opf_against_opendss(model, data)
    # Extract systemName from data dictionary
    systemName = data[:systemName]
    
    # Build the path to the Master.dss file
    filename = joinpath(dirname(@__DIR__), "rawData", systemName, "Master.dss")
    
    # Set up and run the OpenDSS commands
    dss("""
        clear
        redirect "$filename"
        ! solve
    """)

    # Initialize load aggregation variables
    loadnumber = Loads.First()
    kWsum = 0.0
    kvarsum = 0.0
    
    # Sum up kW and kvar across all loads in OpenDSS
    while loadnumber > 0
        kWsum += Loads.kW()
        kvarsum += Loads.kvar()
        loadnumber = Loads.Next()
    end

    # Return the aggregated results
    return kWsum, kvarsum
end

end # module openDSSValidator
