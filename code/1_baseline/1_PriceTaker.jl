using JuMP
using CSV
using DataFrames
using XLSX
using Dates
using Plots
using StatsPlots

using PATHSolver
PATHSolver.c_api_License_SetString("1259252040&Courtesy&&&USR&GEN2035&5_1_2026&1000&PATH&GEN&31_12_2035&0_0_0&6000&0_0")

dir = @__DIR__

include("$dir/1_inputs.jl")

begin
    model = Model(PATHSolver.Optimizer)
    set_optimizer_attribute(model, "major_iteration_limit", 600000)
    set_optimizer_attribute(model, "minor_iteration_limit", 20000000)
    set_optimizer_attribute(model, "cumulative_iteration_limit", 10000000)
    set_optimizer_attribute(model, "crash_iteration_limit", 1000)

    # Decision variables
    @variable(model, S[i in nanogrids] >= 0)                            # Solar capacity (kWp)
    @variable(model, B[i in nanogrids] >= 0)                            # Battery capacity 
    @variable(model, E_sol[t in timesteps, i in nanogrids] >= 0)        # Solar energy used
    @variable(model, E_imp[t in timesteps, i in nanogrids] >= 0)        # Energy imported
    @variable(model, E_exp[t in timesteps, i in nanogrids] >= 0)        # Energy exported
    @variable(model, E_batplus[t in timesteps, i in nanogrids] >= 0)    # Battery charging
    @variable(model, E_batminus[t in timesteps, i in nanogrids] >= 0)   # Battery discharging
    @variable(model, SoC[t in timesteps, i in nanogrids] >= 0)          # State of Charge
    @variable(model, E[t in timesteps, i in nanogrids] >= 0)            # Slack component

    # Set start values
    set_start_value.(S, [8.22, 6.91, 4.01])
    set_start_value.(B, [41.99, 32.40, 21.15])

    # Complementarity slackness variables for upper bounds
    @variable(model, μ_sol[t in timesteps, i in nanogrids] >= 0)
    @variable(model, μ_imp[t in timesteps, i in nanogrids] >= 0)
    @variable(model, μ_exp[t in timesteps, i in nanogrids] >= 0)
    @variable(model, μ_batplus[t in timesteps, i in nanogrids] >= 0)
    @variable(model, μ_batminus[t in timesteps, i in nanogrids] >= 0)
    @variable(model, μ_SoC[t in timesteps, i in nanogrids] >= 0)

    @variable(model, λ_nano[t in timesteps, i in nanogrids] >= 0)
    @variable(model, λ_micro[t in timesteps])
    @variable(model, λ_SoC[t in timesteps, i in nanogrids])

    @expression(model, p_nano[t in timesteps, i in nanogrids],
                       λ_nano[t, i])

    @expression(model, p_micro[t in timesteps],
                       λ_micro[t])

    # Complementarity constraints
    @constraint(model, poeq1[t in timesteps, i in nanogrids],  -λ_nano[t, i] + μ_sol[t, i] ⟂ E_sol[t, i])
    @constraint(model, poeq2[t in timesteps, i in nanogrids], c_bat + λ_nano[t, i] + μ_batplus[t, i] - u_bat* λ_SoC[t, i] ⟂ E_batplus[t, i]) 
    @constraint(model, poeq3[t in timesteps, i in nanogrids], c_bat - λ_nano[t, i] + μ_batminus[t, i] +  1/u_bat* λ_SoC[t, i] ⟂ E_batminus[t, i])
    @constraint(model, poeq4[t in timesteps, i in nanogrids], c_ie - λ_nano[t, i] + p_micro[t] + λ_micro[t]  + μ_imp[t, i] ⟂ E_imp[t, i])
    @constraint(model, poeq5[t in timesteps, i in nanogrids], c_ie + λ_nano[t, i] - p_micro[t] - λ_micro[t] + μ_exp[t, i] ⟂ E_exp[t, i])
    @constraint(model, poeq6[i in nanogrids], c_S[i] - sum(t -> e_solar[t] * μ_sol[t, i], timesteps) ⟂ S[i])
    @constraint(model, poeq7[i in nanogrids], c_B[i] - z_B * λ_SoC[1, i] - v_B * sum(t -> μ_batplus[t, i] + μ_batminus[t, i], timesteps) - sum(t -> μ_SoC[t, i], timesteps) ⟂ B[i])
    @constraint(model, poeq8[t in 1:Num_of_ts-1, i in nanogrids], λ_SoC[t, i]
                                            - (1-sigma) * λ_SoC[t+1, i]
                                            + μ_SoC[t, i] ⟂ SoC[t, i])
    @constraint(model, poeq9[i in nanogrids], λ_SoC[Num_of_ts, i] + μ_SoC[Num_of_ts, i] ⟂ SoC[Num_of_ts, i])

    # Dual derivatives and system balance equations
    @constraint(model, doeq1[t in timesteps], sum(i -> - E_imp[t, i] + E_exp[t, i], nanogrids) ⟂ λ_micro[t])
    @constraint(model, doeq5[i in nanogrids], SoC[1, i] - z_B * B[i] + (1/u_bat) * E_batminus[1, i]
                                            - u_bat * E_batplus[1, i] ⟂ λ_SoC[1, i])
    @constraint(model, doeq6[t in 2:Num_of_ts-1, i in nanogrids], SoC[t, i] - (1-sigma) * SoC[t-1, i]
                                                    + (1/u_bat) * E_batminus[t, i] - u_bat * E_batplus[t, i] ⟂ λ_SoC[t, i])
    @constraint(model, doeq8[i in nanogrids], SoC[Num_of_ts, i] - z_B * B[i] ⟂ λ_SoC[Num_of_ts, i])
 
    @constraint(model, doeq7[t in timesteps, i in nanogrids], E_sol[t, i] + E_batminus[t, i] + E_imp[t, i]
                                    - E_batplus[t, i] - E_exp[t, i] - e_dem_df[t, i] ⟂ λ_nano[t, i])
    # Upper bounds
    @constraint(model, IneqCon1[t in timesteps, i in nanogrids], -E_sol[t, i] + e_solar[t] * S[i] ⟂ μ_sol[t, i])
    @constraint(model, IneqCon2[t in timesteps, i in nanogrids], -E_imp[t, i] + g ⟂ μ_imp[t, i])
    @constraint(model, IneqCon3[t in timesteps, i in nanogrids], -E_exp[t, i] + g ⟂ μ_exp[t, i])
    @constraint(model, IneqCon4[t in timesteps, i in nanogrids], -E_batplus[t, i] + v_B * B[i] ⟂ μ_batplus[t, i])
    @constraint(model, IneqCon5[t in timesteps, i in nanogrids], -E_batminus[t, i] + v_B * B[i] ⟂ μ_batminus[t, i])
    @constraint(model, IneqCon6[t in timesteps, i in nanogrids], -SoC[t, i] + B[i] ⟂ μ_SoC[t, i])

    @expression(model, indiprofit[i in nanogrids], 
        c_S[i] * S[i] + c_B[i] * B[i] +
        sum(t -> c_bat * (E_batplus[t, i] + E_batminus[t, i]), timesteps) +
        sum(t -> c_ie * (E_imp[t, i] + E_exp[t, i]), timesteps) -
        sum(t -> p_nano[t, i] *  e_dem_df[t, i], timesteps)
        + sum(t -> p_micro[t]  * (E_imp[t, i] - E_exp[t, i]), timesteps)
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

    optimize!(model)

    message = termination_status(model)

end

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
λ_micro_values = Array(value.(λ_micro))
p_micro_values = Array(value.(p_micro))
λ_nano_values = [Array(value.(λ_nano[:, i])) for i in nanogrids]
p_nano_values = [Array(value.(p_nano[:, i])) for i in nanogrids]
SoC_values = [Array(value.(SoC[:, i])) for i in nanogrids]
solar_generation = [Array(value.(E_sol[:, i])) for i in nanogrids]
energy_import = [Array(value.(E_imp[:, i])) for i in nanogrids]
energy_export = [Array(value.(E_exp[:, i])) for i in nanogrids]
battery_discharge = [Array(value.(E_batminus[:, i])) for i in nanogrids]
battery_charge = [Array(value.(E_batplus[:, i])) for i in nanogrids]

# Plot λ_micro and p_micro over time
plot_micro = plot(timesteps, λ_micro_values, title="\nλ_micro over Time\n", xlabel="\nTime", ylabel="\nλ_micro", label="λ_micro", lw=2)
plot!(p_micro_values, label="p_micro")
savefig(plot_micro, "$dir/p_&_lambda_micro")

# Plot λ_nano over time for each nanogrid
plot_nano = plot()
for i in 1:length(nanogrids)
    plot!(timesteps, λ_nano_values[i], label="λ_nanogrid $i", lw=2)
    plot!(p_nano_values[i], label="p_nanogrid $i")
end
plot!(title="\n λ_nano and p_nano over Time\n", xlabel="\nTime", ylabel="\nλ_nano")
savefig(plot_nano, "$dir/p_&_lambda_nano")

for i in 1:length(nanogrids)
    plot_nano_NG = plot()
    plot!(timesteps, λ_nano_values[i], label="λ_nanogrid $i", lw=2)
    plot!(p_nano_values[i], label="p_nanogrid $i")
    plot!(title="\n λ_nano and p_nano over Time for Nanogrid $i \n", xlabel="\nTime", ylabel="\nλ_nano")
    savefig(plot_nano_NG, "$dir/p_&_lambda_nano_$i")
end

# Create the Excel output path
output_path = "$dir/Results_NanoGrid_MCP.xlsx"

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
        sheet_investment[l, 5] = nothing # investment!D2
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
    sheet_profit[2, 2] = value(negprofit)               # profit!B2
    for i in nanogrids
        sheet_profit[1 + i, 5] = value(indiprofit[i])    # profit!D2
    end
    sheet_profit[2, 6] = value(totalcost)             # profit!F2
    for i in nanogrids
        sheet_profit[1 + i, 9] = value(indicost[i])     # profit!H2
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
            sheet_p_micro[1 + t, 3] = value(p_micro[t])             # p_micro!B2
            sheet_p_nano[1 + t, 2 + i] = value(p_nano[t, i])        # p_nano!B2
        end 
    end
end