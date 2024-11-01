using OpenDSSDirect

systemName = "ads10_1ph" # this is something which the user will specify but will get saved into data

filename = joinpath(@__DIR__, "rawData", systemName, "Master.dss")

dss("""
    clear
    redirect "$filename"
    ! solve
""")

function main()

    loadnumber = Loads.First()
    kWsum = 0.0
    kvarsum = 0.0
    
    while loadnumber > 0
        kWsum += Loads.kW()
        kvarsum += Loads.kvar()
        loadnumber = Loads.Next()
    end

    kWsum, kvarsum
end

main()