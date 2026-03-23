using JuMP
using GLPK
using CSV
using XLSX
using DataFrames
using Dates
using Plots
using StatsPlots

dir = @__DIR__

include("$dir/4_inputs.jl")

# Create model
begin
    model = Model(GLPK.Optimizer)

    # Variables
    @variable(model, S[i in nanogrids] >= 0)                          # Solar capacity invested (kWp)
    @variable(model, B[i in nanogrids] >= 0)                          # Battery capacity invested (kWh)
    @variable(model, E_sol[t in timesteps, i in nanogrids] >= 0)      # Solar energy used (kWh)
    @variable(model, E_imp[t in timesteps, i in nanogrids] >= 0)      # Energy import (kWh)
    @variable(model, E_exp[t in timesteps, i in nanogrids] >= 0)      # Energy export (kWh)
    @variable(model, E_batplus[t in timesteps, i in nanogrids] >= 0)  # Charge of battery (kWh)
    @variable(model, E_batminus[t in timesteps, i in nanogrids] >= 0) # Discharge of battery (kWh)
    @variable(model, SoC[t in timesteps, i in nanogrids] >= 0)        # State of Charge (kWh)

    @expression(model, p_nano[t in timesteps, i in nanogrids],
                        2)

    set_start_value.(S, [8.22, 6.91, 4.01])
    set_start_value.(B, [41.99, 32.40, 21.15])

    # Objective function: Maximize total profit
    @objective(model, Min, sum(sum(c_S[i] * S[i] + c_B[i] * B[i] +
                            c_bat * (E_batplus[t,i] + E_batminus[t,i]) +
                            c_ie * (E_imp[t,i] + E_exp[t,i]) 
                            - p_nano[t, i] * (E_sol[t,i] + E_batminus[t,i] + E_imp[t,i] - E_batplus[t,i] - E_exp[t,i])
                            for i in nanogrids)
                            for t in timesteps
                                )
                )

    # Equality constraints
    @constraint(model, [t in timesteps], sum(-E_imp[t,i] + E_exp[t,i] for i in nanogrids) == 0)
    @constraint(model, [i in nanogrids], SoC[1,i] - z_B * B[i] - u_bat * E_batplus[1,i] + (1/u_bat) * E_batminus[1,i] == 0.)
    @constraint(model, [t=2:Num_of_ts, i in nanogrids], SoC[t,i] - (1 - sigma) * SoC[t-1,i] - u_bat * E_batplus[t,i] + (1/u_bat) * E_batminus[t,i] == 0.)

    @constraint(model, [i in nanogrids], c_S[i] * S[i] + c_B[i] * B[i] <= m[i])

    # Energy balance constraint
    @constraint(model, [t in timesteps, i in nanogrids], E_sol[t,i] + E_batminus[t,i] + E_imp[t,i] - E_batplus[t,i] - E_exp[t,i] - e_dem_df[t, i] == 0.)

    # Upper bounds
    @constraint(model, [t in timesteps, i in nanogrids], E_sol[t,i] <= e_solar[t] * S[i])
    @constraint(model, [t in timesteps, i in nanogrids], E_imp[t,i] <= g)
    @constraint(model, [t in timesteps, i in nanogrids], E_exp[t,i] <= g)
    @constraint(model, [t in timesteps, i in nanogrids], E_batplus[t,i] <= v_B * B[i])
    @constraint(model, [t in timesteps, i in nanogrids], E_batminus[t,i] <= v_B * B[i])
    @constraint(model, [t in timesteps, i in nanogrids], SoC[t,i] <= B[i])


    @expression(model, indiprofit[i in nanogrids], 
        c_S[i] * S[i] + c_B[i] * B[i] +
        sum(t -> c_bat * (E_batplus[t, i] + E_batminus[t, i]), timesteps) +
        sum(t -> c_ie * (E_imp[t, i] + E_exp[t, i]), timesteps) -
        sum(t -> p_nano[t, i] *  e_dem_df[t, i], timesteps)
    )

    # Define negative profit as the sum of individual profits
    @expression(model, negprofit, sum(indiprofit[i] for i in nanogrids))

    # Define individual cost for each i
    @expression(model, indicost[i in nanogrids], 
        c_S[i] * S[i] + c_B[i] * B[i] +
        sum(t -> c_bat * (E_batplus[t, i] + E_batminus[t, i]), timesteps) +
        sum(t -> c_ie * (E_imp[t, i] + E_exp[t, i]), timesteps)
    )

    # Define total cost as the sum of individual costs
    @expression(model, totalcost, sum(indicost[i] for i in nanogrids))
end 

# Solve the model
optimize!(model)

message = termination_status(model)

# Check if the problem was solved and feasible
@assert is_solved_and_feasible(model) "Termination status: $message"

for i in nanogrids
    println("Nanogrid $i: Solar capacity: ", value(S[i]), " kWp")
    println("Nanogrid $i: Battery capacity: ", value(B[i]), " kWh")
    println("Nanogrid $i: Individual profit:", value(indiprofit[i]))
end

# Print results
println("Negative Profit: ", value(negprofit))
println("Individual Costs: ", value.(indicost))
println("Total Cost: ", value(totalcost))

# Extract optimal values from the model
p_nano_values = [Array(value.(p_nano[:, i])) for i in nanogrids]
SoC_values = [Array(value.(SoC[:, i])) for i in nanogrids]
solar_generation = [Array(value.(E_sol[:, i])) for i in nanogrids]
energy_import = [Array(value.(E_imp[:, i])) for i in nanogrids]
energy_export = [Array(value.(E_exp[:, i])) for i in nanogrids]
battery_discharge = [Array(value.(E_batminus[:, i])) for i in nanogrids]
battery_charge = [Array(value.(E_batplus[:, i])) for i in nanogrids]

# Create the Excel output path
output_path = "$dir/Results_NanoGrid_LP.xlsx"

# Helper function to get or add a sheet
function get_or_add_sheet(xf, sheet_name)
    # Check if the sheet already exists using `sheetnames`
    if sheet_name in XLSX.sheetnames(xf)
        return XLSX.getsheet(xf, sheet_name)   # Return the existing sheet
    else
        return XLSX.addsheet!(xf, sheet_name)  # Add a new sheet
    end
end

# Write data to Excel with overwriting
XLSX.openxlsx(output_path, mode = "rw") do xf
    # Get existing sheets or add new ones if they don't exist
    sheet_solar = get_or_add_sheet(xf, "solar")
    sheet_demand = get_or_add_sheet(xf, "demand")
    sheet_profit = get_or_add_sheet(xf, "profit")
    sheet_investment = get_or_add_sheet(xf, "investment")
    sheet_solgen = get_or_add_sheet(xf, "solgen")
    sheet_import = get_or_add_sheet(xf, "import")
    sheet_export = get_or_add_sheet(xf, "export")
    sheet_charge = get_or_add_sheet(xf, "charge")
    sheet_discharge = get_or_add_sheet(xf, "discharge")
    sheet_soc = get_or_add_sheet(xf, "soc")
    sheet_p_micro = get_or_add_sheet(xf, "p_micro")
    sheet_p_nano = get_or_add_sheet(xf, "p_nano")
    
    # Parameters (starting at B2, which is column 2, row 2)
    # Delete existing values
    for l in 1:1000
        sheet_solar[l, 2] = nothing
        sheet_solar[l, 3] = nothing
        sheet_demand[l, 2] = nothing
        for i in nanogrids
            sheet_demand[l, 2 + i] = nothing  # demand!B2
        end 
    end 

    for t in timesteps
        # Fill with new values 
        sheet_solar[1 + t, 2] = t
        sheet_demand[2 + t, 2] = t 
        sheet_solar[1 + t, 3] = value(e_solar[t])  # solar!B2
        for i in nanogrids
            sheet_demand[2 + t, 2 + i] = value(e_dem[t, i])   # demand!B2
        end 
    end
    
    # Variables (investment results)
    # Delete existing values
    for l in 1:1000
        sheet_investment[l, 3] = nothing  # investment!B2
        sheet_investment[l, 5] = nothing  # investment!D2
    end 

    for i in nanogrids
        # Fill with new values 
        sheet_investment[1 + i, 3] = value(S[i])   # investment!B2
        sheet_investment[1 + i, 5] = value(B[i])   # investment!D2
    end

    # Profit results
    # Delete existing values
    sheet_profit[2, 2] = nothing
    sheet_profit[2, 6] = nothing

    for l in 1:1000
        sheet_profit[l, 5] = nothing
        sheet_profit[l, 9] = nothing
    end 

    # Fill with new values 
    sheet_profit[2, 2] = value(negprofit)                # profit!B2
    for i in nanogrids
        sheet_profit[1 + i, 5] = value(indiprofit[i])    # profit!D2
    end
    sheet_profit[2, 6] = value(totalcost)                # profit!F2
    for i in nanogrids
        sheet_profit[1 + i, 9] = value(indicost[i])      # profit!H2
    end

    # Energy variables (solgen, import, export)
    # Delete existing values
    for l in 1:1000
        sheet_solgen[l, 2] = nothing       # solgen!B2
        sheet_import[l, 2] = nothing       # import!B2
        sheet_export[l, 2] = nothing       # export!B2
        sheet_charge[l, 2] = nothing       # charge!B2
        sheet_discharge[l, 2] = nothing    # discharge!B2
        sheet_soc[l, 2] = nothing          # soc!B2
        sheet_p_micro[l, 2] = nothing      # p_micro!B2
        sheet_p_nano[l, 2] = nothing       # p_nano!B2
        for i in nanogrids
            # Delete existing values
            sheet_solgen[l, 2 + i] = nothing       # solgen!B2
            sheet_import[l, 2 + i] = nothing       # import!B2
            sheet_export[l, 2 + i] = nothing       # export!B2
            sheet_charge[l, 2 + i] = nothing       # charge!B2
            sheet_discharge[l, 2 + i] = nothing    # discharge!B2
            sheet_soc[l, 2 + i] = nothing          # soc!B2
            sheet_p_micro[l, 3] = nothing          # p_micro!B2
            sheet_p_nano[l, 2 + i] = nothing       # p_nano!B2
        end 
    end 

    # Fill with new values
    for t in timesteps
        sheet_solgen[2 + t, 2] = t       # solgen!B2
        sheet_import[2 + t, 2] = t       # import!B2
        sheet_export[2 + t, 2] = t       # export!B2
        sheet_charge[2 + t, 2] = t       # charge!B2
        sheet_discharge[2 + t, 2] = t    # discharge!B2
        sheet_soc[2 + t, 2] = t          # soc!B2
        sheet_p_micro[1 + t, 2] = t      # p_micro!B2
        sheet_p_nano[1 + t, 2] = t       # p_nano!B2

        for i in nanogrids
            # Fill with new values 
            sheet_solgen[2 + t, 2 + i] = value(E_sol[t, i])         # solgen!B2
            sheet_import[2 + t, 2 + i] = value(E_imp[t, i])         # import!B2
            sheet_export[2 + t, 2 + i] = value(E_exp[t, i])         # export!B2
            sheet_charge[2 + t, 2 + i] = value(E_batplus[t, i])     # charge!B2
            sheet_discharge[2 + t, 2 + i] = value(E_batminus[t, i]) # discharge!B2
            sheet_soc[2 + t, 2 + i] = value(SoC[t, i])              # soc!B2
            sheet_p_nano[1 + t, 2 + i] = value(p_nano[t, i])        # p_nano!B2
        end 
    end
end