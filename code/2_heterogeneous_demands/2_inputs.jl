# Load solar generation and demand data
e_solar = CSV.read("$dir/output_pv_15min_2022.csv", DataFrame)[:, 2] 
e_dem = CSV.read("$dir/Nanogrids_demands_365days_15min.csv", DataFrame)[:, 2:4]

nanogrid_count = 3      # Number of nanogrids

Num_of_ts = 7 * 24 * 4  # Number of timesteps
resample_hour = false
factor = 1/4
annualized = Num_of_ts / (7 * 24* 4) 

timesteps = 1:Num_of_ts
nanogrids = 1:nanogrid_count
timesteps_oneday = 1:24*4
# Make demand data a dataframe
begin
    e_dem_df = []
    for i in nanogrids
        global e_dem_df = [e_dem_df 
        DataFrame(
        timestamp = DateTime(2022, 1, 1, 0, 0):Minute(15):DateTime(2022, 12, 31, 23, 45),
        value = e_dem[:, i]
        )]
    end 
end 

if resample_hour
    # Add an "hour" column to group by
    for i in nanogrids
        e_dem_df[i].hour = floor.(e_dem_df[i].timestamp, Dates.Hour)
        rename!(e_dem_df[i], :hour => Symbol("hour_$i"))
    end 

    # Aggregate to hourly data (e.g., summing the values)
    a1 = combine(groupby(e_dem_df[1], :hour_1), :value => sum => :hourly_sum1)
    a2 = combine(groupby(e_dem_df[2], :hour_2), :value => sum => :hourly_sum2)
    a3 = combine(groupby(e_dem_df[3], :hour_3), :value => sum => :hourly_sum3)

    hourly_dem_data = [a1 a2 a3]
    select!(hourly_dem_data, Not([:hour_1, :hour_2, :hour_3]))

    b1 = hourly_dem_data[:, :hourly_sum1]
    b2 = hourly_dem_data[:, :hourly_sum2]
    b3 = hourly_dem_data[:, :hourly_sum3]

    e_dem_df = [b1 b2 b3] .* 1/4

    e_solar = DataFrame(timestamp = DateTime(2022, 1, 1, 0, 0):Minute(15):DateTime(2022, 12, 31, 23, 45),
                        value = e_solar
                        )
    e_solar.hour = floor.(e_solar.timestamp, Dates.Hour)
    e_solar = combine(groupby(e_solar, :hour), :value => sum => :hourly_sum)
    e_solar = e_solar[:, :hourly_sum] # .* 1/4

    factor = 1
else 
    a1 = e_dem_df[1][:, :value]
    a2 = e_dem_df[2][:, :value]
    a3 = e_dem_df[3][:, :value]
    e_dem_df = [a1 a2 a3] .* factor
    e_solar = e_solar * factor
end

# Parameters
c_S = [0.518, 0.518, 0.518] .* annualized  # Cost per unit capacity of PV system investment (euro/kW)
c_B = [0.384, 0.384, 0.384] .* annualized  # Cost per unit capacity of battery system (euro/kWh)
m = [10., 10., 10.]

c_bat = 0.0005      # Marginal cost for battery operation (charge and discharge) (euro/kWh)
c_ie = 0.0001       # Marginal cost of energy sharing
sigma = 0.00001     # Self-discharge rate of battery
u_bat = 0.9         # Charge and discharge efficiency
v_B = 0.5           # Max charge/discharge fraction per timestep
z_B = 0.7           # Starting SoC fraction of installed capacity
g = 100             # Max export and import limit 
