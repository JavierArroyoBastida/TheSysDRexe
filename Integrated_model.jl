
################################################################################
# INTEGRATED MODEL OF ELECTRICITY GENERATION AND THERMOSTATICALLY CONTROLLED LOADS

# Authors: Dieter Patteeuw. Kenneth Bruninx. Glenn Reynders.
# Full model description: https://lirias.kuleuven.be/handle/123456789/545500.
# Version 07 December 2017.

# Adapted for exercise session of thermal systems by Anke Uytterhoeven.
# Version March 2020.
################################################################################

# clearconsole()
println("\n***** INTEGRATED MODEL *****")

@time begin

#_______________________________________________________________________________
## USER DEFINED SETTINGS
#_______________________________________________________________________________
# working directory
# cd("C:/Users/u0110910/workspace/TheSysDRexe")



#_______________________________________________________________________________
## LOADING PACKAGES AND GENERAL FUNCTIONALITIES
#_______________________________________________________________________________

println("----- LOADING PACKAGES AND GENERAL FUNCTIONALITIES -----")

current_dir = pwd()
mkpath(string(current_dir,"/RESULTS"))                                          # create a folder to save the results

# global GUROBI_ENV = Gurobi.Env()

# using PyCall
# using PyPlot
# pygui(true)
# close("all")



#_______________________________________________________________________________
## TIME INVARIANT SETTINGS - GENERAL
#_______________________________________________________________________________

println("----- PROCESSING GENERAL TIME INVARIANT SETTINGS -----")

#...............................................................................
# Read input data from Excel file
mainNames    = XLSX.readdata("InputData.xlsx","MainParameters","A2:A20")        # Read the nametags
mainValues   = XLSX.readdata("InputData.xlsx","MainParameters","B2:B20")        # Read the values

supplyNames  = XLSX.readdata("InputData.xlsx","SupplyParameters","A2:A6")       # Read the nametags
supplyValues = XLSX.readdata("InputData.xlsx","SupplyParameters","B2:B6")       # Read the values

demandNames  = XLSX.readdata("InputData.xlsx","DemandParameters","A2:A41")      # Read the nametags
demandValues = XLSX.readdata("InputData.xlsx","DemandParameters","B2:B41")      # Read the values

timeNames    = XLSX.readdata("InputData.xlsx","YearData","A2:Y2")               # Read the nametags
timeValues   = XLSX.readdata("InputData.xlsx","YearData","A3:Y8954")            # Read the values

# Dictionary for reading the main parameters, supply parameters and demand parameters
mainParams = Dict(zip(mainNames,mainValues))
delete!(mainParams,missing)
supplyParams = Dict(zip(supplyNames,supplyValues))
delete!(supplyParams,missing)
demandParams = Dict(zip(demandNames,demandValues))
delete!(demandParams,missing)
mainDict = merge(mainParams,supplyParams,demandParams)

# DataFrame for the time series data
timeSeries = DataFrame(timeValues)
rename!(timeSeries,[Symbol(timeNames[i]) for i=1:length(timeNames)])

#...............................................................................
# Simulation time settings
heatingSeasonVec = [mainDict["StartHeatingDay"]:mainDict["OptDays"]:364 ; 1:mainDict["OptDays"]:mainDict["EndHeatingDay"]]
if mainDict["SimFullYear"]
    simulTimeVec = [mainDict["StartHeatingDay"]:mainDict["OptDays"]:364 ; 1:mainDict["OptDays"]:mainDict["StartHeatingDay"]-mainDict["OptDays"]]
else
    simulTimeVec = mainDict["SimStartDay"]
end

#...............................................................................
# Tracking results
if mainDict["SimFullYear"]
    resHours = mainDict["fullyearEffective"]                                    # 8760
else
    resHours = mainDict["TimeSteps"]
end
results_format = DataFrame(cost_tot    = zeros(resHours),
                           gen_ccgt    = zeros(resHours),
                           gen_ocgt    = zeros(resHours),
                           gen_solar   = zeros(resHours),
                           gen_wind    = zeros(resHours),
                           curtail     = zeros(resHours),
                           load_shed   = zeros(resHours),
                           dem_fxd     = zeros(resHours),
                           dem_var     = zeros(resHours),
                           cost_ccgt   = zeros(resHours),
                           cost_ocgt   = zeros(resHours),
                           cost_LL     = zeros(resHours),
                           dem_var_h1  = zeros(resHours),
                           dem_var_h2  = zeros(resHours),
                           Php_sh_h1   = zeros(resHours),
                           Php_sh_h2   = zeros(resHours),
                           Php_dhw_h1  = zeros(resHours),
                           Php_dhw_h2  = zeros(resHours),
                           Paux_sh_h1  = zeros(resHours),
                           Paux_sh_h2  = zeros(resHours),
                           Paux_dhw_h1 = zeros(resHours),
                           Paux_dhw_h2 = zeros(resHours),
                           QemD_h1     = zeros(resHours),
                           QemD_h2     = zeros(resHours),
                           QemN_h1     = zeros(resHours),
                           QemN_h2     = zeros(resHours),
                           Qtank_h1    = zeros(resHours),
                           Qtank_h2    = zeros(resHours),
                           Tday_h1     = zeros(resHours),
                           Tday_h2     = zeros(resHours),
                           Tnight_h1   = zeros(resHours),
                           Tnight_h2   = zeros(resHours),
                           Xw_h1       = zeros(resHours),
                           Xw_h2       = zeros(resHours),
                           Tmin_h1     = zeros(resHours),
                           Tmin_h2     = zeros(resHours),
                           Tmax_h1     = zeros(resHours),
                           Tmax_h2     = zeros(resHours),
                           Xmin_h1     = zeros(resHours),
                           Xmin_h2     = zeros(resHours),
                           Xmax_h1     = zeros(resHours),
                           Xmax_h2     = zeros(resHours),
                           dem_noSH_fxd = zeros(resHours),
                           dem_HP_fxd  = zeros(resHours))

global results_tot = deepcopy(results_format)



#_______________________________________________________________________________
## TIME INVARIANT SETTINGS - SUPPLY SIDE MODEL
#_______________________________________________________________________________

println("----- PROCESSING TIME INVARIANT SETTINGS OF SUPPLY SIDE MODEL -----")

#...............................................................................
# Read out time invariant characteristics of considered power plants (only minimum and maximum operating point and associated CO2 and fuel cost
gen_max_ccgt = mainDict["CapCCGT"]                                               # [MW] maximum power of the combined cycle gas turbine power plant, put equal to installed capacity
gen_max_ocgt = mainDict["CapOCGT"]                                               # [MW] maximum power of the open cycle gas turbine power plant, put equal to installed capacity
MC_ccgt    = ( mainDict["fuelcost"] + mainDict["CO2_price"]*mainDict["fuelco2"] ) /mainDict["effCCGT"]   # marginal cost of producing one MWh_electrical with a combined cycle gas turbine (taking into account fuel cost and CO2 cost)
MC_ocgt    = ( mainDict["fuelcost"] + mainDict["CO2_price"]*mainDict["fuelco2"] ) /mainDict["effOCGT"]   # marginal cost of producing one MWh_electrical with an open cycle gas turbine (taking into account fuel cost and CO2 cost)

VOLL = mainDict["VOLL"]                                                         # value of lost load, if power plant is not able to deliver requested electricity



#_______________________________________________________________________________
## TIME INVARIANT SETTINGS - DEMAND SIDE MODEL
#_______________________________________________________________________________

println("----- PROCESSING TIME INVARIANT SETTINGS OF DEMAND SIDE MODEL -----")

#...............................................................................
# Resizing of the demand from the modeled houses (10^(-6) because the demand side model is in W and the supply model in MW)
factorBuildingsDR = mainDict["share_HP_withDR"] * 10^(-6) * mainDict["nbuildsStudy"]   # take fraction of buildings into account
ResizeVecBuildings = factorBuildingsDR*[0.5 ; 0.5]                                     # equal division of building between first and second type

#...............................................................................
# System dynamics space heating

# 9 state State Space Model (SSM)
# AbuildDisc, BbuildDisc and EbuildDisc are in discrete timedomain with timestep 1 hour
# model: X(i+1) = A*X(i) + B*U(i) + E*D(i)  with X states, U inputs and D disturbances
# 1   X = [ TiD  ]    [°C] day zone temperature
# 2       [ TwD  ]    [°C] external walls day temperature
# 3       [ TwiD ]    [°C] internal walls  day temperature
# 4       [ TflD ]    [°C] floor day temperature
# 5       [ TfDN ]    [°C] floor day to night temperature dayside
# 6       [ TiN  ]    [°C] night-zone temperture
# 7       [ TwN  ]    [°C] external walls night temperature
# 8       [ TwiN ]    [°C] internal walls night temperature
# 9       [ TfND ]    [°C] floor day to night temperature nightside
#
#   and inputs:
# 1   U = [ QemD]     [W]  Heat added to the heat emission system of the day zone
# 2       [ QemN]     [W]  Heat added to the heat emission system of the night zone
#
#   and disturbances:
# 1   E = [ Te    ]   [°C] Ambient air temperature
# 2       [ Tgr   ]   [°C] Ground temperature
# 3       [ QsolN ]   [W]  Solar gain from north
# 4       [ QsolE ]   [W]  Solar gain from east
# 5       [ QsolS ]   [W]  Solar gain from south
# 6       [ QsolW ]   [W]  Solar gain from west
# 7       [ QintD ]   [W]  Internal heat gains (user and appliances) to the day zone
# 8       [ QintN ]   [W]  Internal heat gains (user and appliances) to night zone

AshS         = zeros(mainDict["SH_numberStates"],mainDict["SH_numberStates"],mainDict["numberBuildings"])
BshS         = zeros(mainDict["SH_numberStates"],mainDict["SH_numberInputs"],mainDict["numberBuildings"])
EshS         = zeros(mainDict["SH_numberStates"],mainDict["SH_numberDisturbances"],mainDict["numberBuildings"])

rangeAss  = XLSX.readdata("InputData.xlsx","BuildingModel","B2:B4")             # Range of the A,B,E matrices data of the state space model of the building

# Read out the building model
AbuildCont      = convert(Array{Float64},XLSX.readdata("InputData.xlsx","BuildingModel",rangeAss[1]))   # A matrix building, continuous model
BbuildCont      = convert(Array{Float64},XLSX.readdata("InputData.xlsx","BuildingModel",rangeAss[2]))   # B matrix building, continuous model
EbuildCont      = convert(Array{Float64},XLSX.readdata("InputData.xlsx","BuildingModel",rangeAss[3]))   # E matrix building, continuous model

# Discretize the state space model
Ccont         = one(AbuildCont)
DcontB        = zeros(size(BbuildCont))
DcontE        = zeros(size(EbuildCont))
systemContB   = ss(AbuildCont,BbuildCont,Ccont,DcontB)                          # construct system with A and B
systemContE   = ss(AbuildCont,EbuildCont,Ccont,DcontE)                          # construct system with A and E
systemDiscB   = c2d(systemContB,mainDict["TimeStepLength"])
systemDiscB   = systemDiscB[1]
systemDiscE   = c2d(systemContE,mainDict["TimeStepLength"])
systemDiscE   = systemDiscE[1]
AbuildDisc    = systemDiscB.A
BbuildDisc    = systemDiscB.B
EbuildDisc    = systemDiscE.B

# The two archetype buildings only differ from each other by their domestic hot water tank; hence, they are both characterized by the same space heating model
AshS[:,:,1] = AshS[:,:,2] = AbuildDisc
BshS[:,:,1] = BshS[:,:,2] = BbuildDisc
EshS[:,:,1] = EshS[:,:,2] = EbuildDisc

#...............................................................................
# System dynamics domestic hot water (DHW) tank

# Model of a perfectly stratified, heat pump heated TES (thermal energy storage) for DHW (domestic hot water)
#
# Assumptions
# - perfectly mixed thermal storage
# - surroundings temperature around the storage constant
# - hot water demand at constant temperature
#
# Adhw and Bdhw are in discrete timedomain with timestep 1 hour
# model: X(i+1) = A*X(i) + B*U(i) + E*D(i)  with X states, U inputs and D disturbances
# 1   X = [ Tdhw  ]   [°C] Domestic hot water temperature inside tank
#
#   and inputs:
# 1   U = [ Qtank]    [W]  Heat added to the tank
#
#   and disturbances:
# 1   E = [ Qdhw  ]   [W]  Domestic hot water demand
# 2       [ Tsurr ]   [°C] Temperature of tank surroundings
#
# rows of models: if 2 DWH tanks, then 2 rows

AdhwS        = zeros(mainDict["DHW_numberStates"],mainDict["DHW_numberStates"],mainDict["numberBuildings"])
BdhwS        = zeros(mainDict["DHW_numberStates"],mainDict["DHW_numberInputs"],mainDict["numberBuildings"])
EdhwS        = zeros(mainDict["DHW_numberStates"],mainDict["DHW_numberDisturbances"],mainDict["numberBuildings"])

dhwData     = convert(Matrix{Float64},XLSX.readdata("InputData.xlsx","DHWtank","B3:E4")) # Range of the A,B,E matrices data of the state space model of the DHW tank
Cdhw        = 1                                                                 # Dummy matrix for state space model operation
Ddhw        = 0                                                                 # Dummy matrix for state space model operation

# First DHW tank model
sys            = ss(dhwData[1,1],dhwData[1,2],Cdhw,Ddhw)                        # Continuous DHW tank model for inputs
sys_d          = c2d(sys,mainDict["TimeStepLength"])                            # Discretize with Ts
sys_d          = sys_d[1]
AdhwS[:,:,1]   = sys_d.A
BdhwS[:,:,1]   = sys_d.B
sys            = ss(dhwData[1,1],transpose(dhwData[1,3:4]),Cdhw,Ddhw)           # Continuous DHW tank model for disturbances
sys_d          = c2d(sys,mainDict["TimeStepLength"])                            # Discretize with Ts
sys_d          = sys_d[1]
EdhwS[:,:,1]   = sys_d.B

# Second DHW tank model
sys            = ss(dhwData[2,1],dhwData[2,2],Cdhw,Ddhw)                        # Continuous DHW tank model for inputs
sys_d          = c2d(sys,mainDict["TimeStepLength"])                            # Discretize with Ts
sys_d          = sys_d[1]
AdhwS[:,:,2]   = sys_d.A
BdhwS[:,:,2]   = sys_d.B
sys            = ss(dhwData[2,1],transpose(dhwData[2,3:4]),Cdhw,Ddhw)           # Continuous DHW tank model for disturbances
sys_d          = c2d(sys,mainDict["TimeStepLength"])                            # Discretize with Ts
sys_d          = sys_d[1]
EdhwS[:,:,2]   = sys_d.B


#_______________________________________________________________________________
## PROCESSING THE TIME VARIANT SETTINGS OF THE MODELS AND RUNNING THE MODEL
#_______________________________________________________________________________

println("----- RUNNING THE MODEL -----")

for day_start = simulTimeVec
    println("Starting the simulation for day $day_start")

    #...............................................................................
    # Time series indices to read out
    time_indices_here = collect((day_start*24-23):1:((day_start*24-23)+mainDict["TimeSteps"]-1))


    #...............................................................................
    # [SUPPLY SIDE] Unit commitment and economic dispatch model

    ## Read time series data
    gen_wind_on    = mainDict["CapWindOn"]*timeSeries[time_indices_here,:g_wind_on]
    gen_wind_off   = mainDict["CapWindOff"]*timeSeries[time_indices_here,:g_wind_off]
    gen_wind       = gen_wind_on + gen_wind_off
    gen_solar      = mainDict["CapSol"]*timeSeries[time_indices_here,:g_solar]

    #...............................................................................
    # [DEMAND SIDE] Demand side model

    # [1] Fixed demand

    ## Read time series data
    # electricity demand for all purposes except for space heating
    demand_noSH_fxd = timeSeries[time_indices_here,:dem_fix]
    # electricity demand for space heating for the non flexible heat pumps without demand response
    demand_HP_fxd = mainDict["share_HP_noDR"] * 10^(-6) * mainDict["nbuildsStudy"] * timeSeries[time_indices_here,:Ptot1house] # 10^(-6) because the demand side model is in W and the supply model in MW
    # total fixed electricity demand
    demand_fxd = demand_noSH_fxd + demand_HP_fxd


    # [2] Flexible demand => calculated by solving integrated model

    ## Constructs the demand side model, as part of the integrated model.
    ## Iterates over the number of dwellings considered.
    ## Constructs the constraints from the dwelling type, occupancy profile and DHW consumption profile.

    #-------------------------------------------------------------------------------
    # Read time series data (boundary conditions/disturbances) used for both house types
    ## Meteo data
    Tamb         =   transpose(timeSeries[time_indices_here,:Tamb])
    QsolN        =   transpose(timeSeries[time_indices_here,:QsolN])
    QsolE        =   transpose(timeSeries[time_indices_here,:QsolE])
    QsolS        =   transpose(timeSeries[time_indices_here,:QsolS])
    QsolW        =   transpose(timeSeries[time_indices_here,:QsolW])
    Tground      =   mainDict["Tground"]*ones(size(Tamb))

    #-------------------------------------------------------------------------------
    # Starting conditions
    if ( mainDict["SimFullYear"] == 0 ) | ( day_start == mainDict["StartHeatingDay"] )
        # If simulation is only 1 period or if the start of the heating season: cyclic boundary condtion: T_start = T_end
        global cyclic_here = true
        # dummy values for start states:
        TshStartS  = DataFrame(zeros(2,11))
        rename!(TshStartS,[:TiD, :TwD, :TwiD, :TflD, :TfDN, :TiN, :TwN, :TwiN, :TfND, :TemD, :TemN])
        XdhwStartS = DataFrame(zeros(2,1))
        rename!(XdhwStartS,[:Xdhw])
    else
        # other cases: starting conditions are taken from the previous simulation
        global cyclic_here = false
        startValues = CSV.read("RESULTS/startValues.csv")

        TshStartS = startValues[:,[:TiD, :TwD, :TwiD, :TflD, :TfDN, :TiN, :TwN, :TwiN, :TfND]]
        XdhwStartS = startValues[:,[:Xdhw]]
    end

    TshStartS = Matrix(TshStartS)
    XdhwStartS = Matrix(XdhwStartS)

    #-------------------------------------------------------------------------------
    # Space heating model (initialization; to be determined for every distinct house type)
    global TminS        = zeros(mainDict["SH_numberStates"],mainDict["TimeSteps"],mainDict["numberBuildings"])
    global TmaxS        = zeros(mainDict["SH_numberStates"],mainDict["TimeSteps"],mainDict["numberBuildings"])
    global COPshS       = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global PhpmaxvecS   = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global PauxmaxvecS  = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global P_sh_100S    = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global P_dhw_100S   = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global DisturbSHS   = zeros(mainDict["SH_numberDisturbances"],mainDict["TimeSteps"],mainDict["numberBuildings"])

    #-------------------------------------------------------------------------------
    # Domestic hot water model (initialization; to be determined for every distinct house type)
    global XminS           = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global XmaxS           = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global COPdhwS         = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global QdhwS           = zeros(mainDict["TimeSteps"],mainDict["numberBuildings"])
    global DisturbDHWS     = zeros(mainDict["DHW_numberDisturbances"],mainDict["TimeSteps"],mainDict["numberBuildings"])

    #-------------------------------------------------------------------------------
    # Read out data for this specific dwelling s
    for s = 1:mainDict["numberBuildings"]

            #-------------------------------------------------------------------------------
            # SPACE HEATING MODEL

            # Temperature constraints
            Tmin      = mainDict["TextremeLow"]*ones(mainDict["SH_numberStates"],mainDict["TimeSteps"])
            Tmax      = mainDict["TextremeHigh"]*ones(mainDict["SH_numberStates"],mainDict["TimeSteps"])
            ## overwrite for comfort limits
            Tmin[1,:] = timeSeries[time_indices_here,Symbol(string("TminDZ",s))]
            Tmin[6,:] = timeSeries[time_indices_here,Symbol(string("TminNZ",s))]
            Tmax[1,:] = timeSeries[time_indices_here,Symbol(string("TmaxDZ",s))]
            Tmax[6,:] = timeSeries[time_indices_here,Symbol(string("TmaxNZ",s))]

            # Heating season?
            if Bool(sum(heatingSeasonVec .== day_start))
                global summerBoolean = false  # heating season
            else
                global summerBoolean = true  # summer season
            end

            # Heat pump model
            ## Determine COP based on the nominal heat demand of that day
            QradDthisCase = timeSeries[time_indices_here,:QradMean]             # [W] Radiator heat demand from the data. Used to determine the supply temperature.
            QradthisCase_nonZero = QradDthisCase[QradDthisCase.>0]
            if ( isempty(QradDthisCase[QradDthisCase.>0]) | (mean(mean(QradDthisCase))==0) )
                QmeanRadD = 100
            else
                QradDnonZero = mean(QradDthisCase[QradDthisCase.>0])
                QmeanRadD = mean(QradDnonZero)
            end

            ## Floor heating model
            Tsupply = 20 + mainDict["deltaTsh"]/2 + QmeanRadD/mainDict["UAfloorheatingD"]

            ## Air coupled heat pump model
            aboveBivalence = Tamb .>= mainDict["ecofysBivalenceStart"]
            aboveBivalence = 1aboveBivalence
            if mean( aboveBivalence) >= 23/24
                # this case, almost all temperatures above Tbivalence
                COP_shFirst   = mainDict["ecofysACHP_A"]*( (Tsupply+273.15) ./ (Tsupply.-Tamb.+mainDict["ecofysACHP_B"]) )
                # for DHW
                COP_dhwFirst  = mainDict["ecofysACHP_A"]*( (mainDict["TsupDhwCop"]+273.15) ./ (mainDict["TsupDhwCop"].-Tamb.+mainDict["ecofysACHP_B"]) )
            else
                COP_shFirst   = zeros(1,mainDict["TimeSteps"])
                COP_dhwFirst  = zeros(1,mainDict["TimeSteps"])
                COP_bivalence = mainDict["ecofysACHP_A"]*( (Tsupply+273.15) / (Tsupply-mainDict["ecofysBivalenceStart"]+mainDict["ecofysACHP_B"]) )
                COP_bivdhw = mainDict["ecofysACHP_A"]*( (mainDict["TsupDhwCop"]+273.15) / (mainDict["TsupDhwCop"]-mainDict["ecofysBivalenceStart"]+mainDict["ecofysACHP_B"]) )
                for iii = 1:mainDict["TimeSteps"]
                    if Tamb[1,iii] > mainDict["ecofysBivalenceStart"]
                        COP_shFirst[1,iii] = mainDict["ecofysACHP_A"]*( (Tsupply+273.15) / (Tsupply-Tamb[1,iii]+mainDict["ecofysACHP_B"]) )
                        # for DHW
                        COP_dhwFirst[1,iii]= mainDict["ecofysACHP_A"]*( (mainDict["TsupDhwCop"]+273.15) / (mainDict["TsupDhwCop"]-Tamb[1,iii]+mainDict["ecofysACHP_B"]) )
                    else
                        # linear interpolation of COP
                        COP_shFirst[1,iii] = mainDict["ecofysCOPatBivalenceStop"] + (COP_bivalence-mainDict["ecofysCOPatBivalenceStop"])*
                            (Tamb[1,iii] -(mainDict["ecofysBivalenceStop"]))/(-mainDict["ecofysBivalenceStart"] -(mainDict["ecofysBivalenceStop"]))
                        # for DHW
                        COP_dhwFirst[1,iii] = mainDict["ecofysCOPatBivalenceStop"] + (COP_bivdhw-mainDict["ecofysCOPatBivalenceStop"])*
                            (Tamb[1,iii] -(mainDict["ecofysBivalenceStop"]))/(-mainDict["ecofysBivalenceStart"] -(mainDict["ecofysBivalenceStop"]))
                    end
                end
            end

            ## For design purposes
            Tamb_DOT    = -10                                                   # [°C] lowest Tamb for which the heat pump is designed
            COP_bivDOT  = mainDict["ecofysACHP_A"]*( (mainDict["TsupShCop"]+273.15)/(mainDict["TsupShCop"]-mainDict["ecofysBivalenceStart"]+mainDict["ecofysACHP_B"]) )
            COP_DOT     = mainDict["ecofysCOPatBivalenceStop"] + (COP_bivDOT-mainDict["ecofysCOPatBivalenceStop"])*
                          (Tamb_DOT-(mainDict["ecofysBivalenceStop"]))/(-mainDict["ecofysBivalenceStart"]-(mainDict["ecofysBivalenceStop"]))

            COP_shClipped   = min.(COP_shFirst,mainDict["ecofysCOPmax"]*ones(size(COP_shFirst)) )
            COP_dhwClipped  = min.(COP_dhwFirst,mainDict["ecofysCOPmax"]*ones(size(COP_dhwFirst)) )

            ## Put constant?
            if mainDict["constantCOP"]
                COP_sh = mean(COP_shClipped)*ones(size(COP_shClipped))
                COPdhw = mean(COP_dhwClipped)*ones(size(COP_dhwClipped))
            else
                COP_sh = COP_shClipped
                COPdhw = COP_dhwClipped
            end

            ## Maximal powers, from fitting, Tsupply = 45°C  NEEDS TO BE RESIZED ACCORDING TO DESIGN HEAT DEMAND
            P_sh_100data        = 24.8 .*Tamb .+ 2745.2                         # Power at 100% modulation. Correlation to Tamb, R² = 0,9536
            P_sh_100DOTdata     = 24.8 .* mainDict["TambDOT"] .+ 2745.2	        # Heat pump power from data at DOT
            P_dhw_100data       = 39.0 .* Tamb .+ 3406.5000                     # Power at 100% modulation Fit of extrapolated data at Tsup of 60°C with R²= 0,9698 wrt Tamb

            ## System constraints
            Qhpmax          = mainDict["heatpumpSizing"]*mainDict["QDOT"]
            Qauxmax        = mainDict["auxHeatingSizing"]*mainDict["QDOT"]
            Phpmax          = Qhpmax/COP_DOT
            Phpmaxvec       = Phpmax*ones(size(Tamb))
            Pauxmaxvec      = Qauxmax*ones(size(Tamb))                         # electrical heating efficiency 1

            ## More detailled Pmax constraint: still need to resize it to design conditions
            resizeFactorWrtData = Phpmax/P_sh_100DOTdata
            P_sh_100    = resizeFactorWrtData*P_sh_100data
            P_dhw_100   = resizeFactorWrtData*P_dhw_100data

            # Occupancy and internal gains
            QintDay   = mainDict[string("Qintdayzone",s)]*ones(1,mainDict["TimeSteps"])      # [W] internal heat gains due to appliances
            QintNight = mainDict["constantHeatGainNightZone"]*ones(1,mainDict["TimeSteps"])  # [W] internal heat gains due to appliances

            # Disturbances
            DisturbSH = [Tamb; Tground; QsolN; QsolE; QsolS; QsolW; QintDay; QintNight]      # Disturbances vector

            # Put the data per house
            global TminS[:,:,s] = Tmin
            global TmaxS[:,:,s] = Tmax
            global COPshS[:,s]  = COP_sh
            global PhpmaxvecS[:,s]  = Phpmaxvec
            global PauxmaxvecS[:,s] = Pauxmaxvec
            global P_sh_100S[:,s]   = P_sh_100
            global P_dhw_100S[:,s]  = P_dhw_100
            global DisturbSHS[:,:,s] = DisturbSH

            #-------------------------------------------------------------------------------
            # DOMESTIC HOT WATER TANK MODEL

            # Disturbances

            # Minimum tank temperature
            Xmin = timeSeries[time_indices_here,Symbol(string("TminDHW",s))]                   # (°C) minimum temperature of the DHW tank in order to have temperature for comfort (aggregated)
            Xmax = mainDict["TmaxDHW"]*ones(mainDict["TimeSteps"],1)                           # (°C) maximum temperature of the DHW tank

            # Energy demand for domestic hot water
            DHW_demand_lmin = timeSeries[time_indices_here,Symbol(string("mDHW",s))]                # (l/min) domestic hot water demand at 38°C in liter per minute (aggregated)
            DHW_demand_lhour =  DHW_demand_lmin*60                                                  # [liter/hour] at 38°C
            DHW_demand = DHW_demand_lhour * (mainDict["dhwReferenceTemperature"] - mainDict["coldTapWaterTemperature"]) /
                         (mainDict["dhwTappingTemperature"] - mainDict["coldTapWaterTemperature"])  # [liter/hour] at the demanded temperature (50°C)
            demand_factor = mainDict["densityWaterKgLiter"]*mainDict["heatCapacityWater"]*(mainDict["dhwTappingTemperature"] - mainDict["coldTapWaterTemperature"])  # [J/liter]
            Qdhw = transpose(1/mainDict["TimeStepLength"] * demand_factor * DHW_demand)             # in W per family

            # Temperature of environment surrounding tank
            Tsurr = mainDict["dhwTankSurroundings"]*transpose(ones(size(DHW_demand)))

            # Disturbances
            DisturbDHW = [Qdhw; Tsurr]                                                          # Disturbances vector

            # Put the data per house
            global XminS[:,s]      = Xmin
            global XmaxS[:,s]      = Xmax
            global COPdhwS[:,s]    = COPdhw
            global QdhwS[:,s]      = Qdhw
            global DisturbDHWS[:,:,s] = DisturbDHW

    end

    #...............................................................................
    # [INTEGRATED MODEL] OVERALL OPTIMIZATION PROBLEM

    # Initialize optimization problem
    # integratedModel = Model(with_optimizer(Ipopt.Optimizer,OutputFlag=0))
    integratedModel = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))


    # Indexing parameters
    timeStep = mainDict["TimeStep"]

    nbStatesSH = mainDict["SH_numberStates"]
    nbInputsSH = mainDict["SH_numberInputs"]
    nbDisturbSH = mainDict["SH_numberDisturbances"]

    nbStatesDhw = mainDict["DHW_numberStates"]
    nbInputsDhw = mainDict["DHW_numberInputs"]
    nbDisturbDhw = mainDict["DHW_numberDisturbances"]

    nbHouses = mainDict["numberBuildings"]
    nbTimesteps = mainDict["TimeSteps"]


    # Optimization variables
    @variables(integratedModel,begin
        gen_ccgt[t=1:nbTimesteps] >= 0                                          # electricity generatad by combined cycle gas turbine           [MW]
        gen_ocgt[t=1:nbTimesteps] >= 0                                          # electricity generatad by open cycle gas turbine               [MW]
        curtail[t=1:nbTimesteps] >= 0                                           # curtailment                                                   [MW]
        load_shed[t=1:nbTimesteps] >= 0                                         # load shedding                                                 [MW]
        demand_var[t=1:nbTimesteps] >= 0                                        # resized variable demand for all households altogether         [MW]
        demand_var_onehouse[t=1:nbTimesteps,s=1:nbHouses] >= 0                  # electricity consumption per house                             [W]
        Tsh[q=1:nbStatesSH,t=1:nbTimesteps,s=1:nbHouses] >= 0                   # state variables of space heating model                        [°C]
        Xdhw[t=1:nbTimesteps,s=1:nbHouses] >= 0                                 # state variable of DHW tank model                              [°C]
        Php_sh[t=1:nbTimesteps,s=1:nbHouses] >= 0                               # power of the heat pump for space heating (sh)                 [W]
        Php_dhw[t=1:nbTimesteps,s=1:nbHouses] >= 0                              # power of the heat pump for domestic hot water (dhw)           [W]
        Paux_sh[t=1:nbTimesteps,s=1:nbHouses] >= 0                              # power of the auxiliary heater for space heating (sh)          [W]
        Paux_dhw[t=1:nbTimesteps,s=1:nbHouses] >= 0                             # power of the auxiliary heater for domestic hot water (dhw)    [W]
        QemD[t=1:nbTimesteps,s=1:nbHouses] >= 0                                 # heat for heat emission system day zone                        [W]
        QemN[t=1:nbTimesteps,s=1:nbHouses] >= 0                                 # heat for heat emission system night zone                      [W]
        Qtank[t=1:nbTimesteps,s=1:nbHouses] >= 0                                # heat for domestic hot water tank                              [W]
    end)


    # Objective function
    ## Total cost for the system
    @objective(integratedModel,Min,
               timeStep*sum(MC_ccgt*gen_ccgt[t] + MC_ocgt*gen_ocgt[t] + VOLL*load_shed[t] for t=1:nbTimesteps))


    # Contraints
    ## Market clearing (balance between supply and demand)
    @constraint(integratedModel,[t=1:nbTimesteps],
               (gen_ccgt[t] + gen_ocgt[t]) +                                    # electricity generation by conventional power plants
               (gen_wind[t] + gen_solar[t] - curtail[t]) +                      # electricity generation by renewable power plants (minus curtailment)
               load_shed[t] -                                                   # shedded load
               (demand_fxd[t] + demand_var[t]) == 0)                            # fixed and variable demand (the latter can be altered via demand reponse)

   ## Curtailment RES
   @constraint(integratedModel,[t=1:nbTimesteps],
               curtail[t] <= gen_wind[t] + gen_solar[t])

   ## Conventional power plants
   @constraint(integratedModel,[t=1:nbTimesteps],
               gen_ccgt[t] <= gen_max_ccgt)
   @constraint(integratedModel,[t=1:nbTimesteps],
               gen_ocgt[t] <= gen_max_ocgt)

   ## Demand
   @constraint(integratedModel,[t=1:nbTimesteps],                               # coupling supply side model (in MW, all houses) and demand side model (in W, two representative archetypes)
               demand_var[t] == sum(demand_var_onehouse[t,s]*ResizeVecBuildings[s] for s=1:nbHouses))

   # IN CASE OF FLEXIBLE DEMAND
   if ~(mainDict["share_HP_withDR"] == 0)

       # IN WINTER
       if ~summerBoolean

            ## Power consumption (in winter: space heating and domestic hot water production)
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        demand_var_onehouse[t,s] == Php_sh[t,s] + Php_dhw[t,s] + Paux_sh[t,s] + Paux_dhw[t,s])
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        Php_sh[t,s] <= P_sh_100S[t,s])
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        Php_dhw[t,s] <= P_dhw_100S[t,s])
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        Php_sh[t,s] + Php_dhw[t,s] <= 0.5*P_sh_100S[t,s] + 0.5*P_dhw_100S[t,s])
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        Paux_sh[t,s] + Paux_dhw[t,s] <= PauxmaxvecS[t,s])

            ## Heating system model
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        QemD[t,s] + QemN[t,s] == COPshS[t,s]*Php_sh[t,s] + Paux_sh[t,s])
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        Qtank[t,s] == COPdhwS[t,s]*Php_dhw[t,s] + Paux_dhw[t,s])


            ## Building model
            if cyclic_here
                @constraint(integratedModel,[q=1:nbStatesSH,s=1:nbHouses],
                            Tsh[q,1,s] == Tsh[q,end,s])                         # Cyclic boundary conditions
            else
                @constraint(integratedModel,[q=1:nbStatesSH,t=1,s=1:nbHouses],  # Results of initialization with cyclic boundary conditions or results of previous optimization problem (previous optimization horizon)
                            Tsh[q,t,s] == sum(AshS[q,qi,s]*TshStartS[s,qi] for qi=1:nbStatesSH) +
                                          (BshS[q,1,s]*QemD[t,s] + BshS[q,2,s]*QemN[t,s]) +
                                          sum(EshS[q,di,s]*DisturbSHS[di,t,s] for di=1:nbDisturbSH))
            end
            @constraint(integratedModel,[q=1:nbStatesSH,t=2:nbTimesteps,s=1:nbHouses],         # Results of previous time step as starting point
                        Tsh[q,t,s] == sum(AshS[q,qi,s]*Tsh[qi,t-1,s] for qi=1:nbStatesSH) +
                                      (BshS[q,1,s]*QemD[t,s] + BshS[q,2,s]*QemN[t,s]) +
                                      sum(EshS[q,di,s]*DisturbSHS[di,t,s] for di=1:nbDisturbSH))

            @constraint(integratedModel,[q=1:nbStatesSH,t=1:nbTimesteps,s=1:nbHouses],
                         TminS[q,t,s] <= Tsh[q,t,s] <= TmaxS[q,t,s])

            ## Domestic hot water tank model
            if cyclic_here
                @constraint(integratedModel,[s=1:nbHouses],
                            Xdhw[1,s] == Xdhw[end,s])                           # Cyclic boundary conditions
            else
                @constraint(integratedModel,[t=1,s=1:nbHouses],                 # Results of initialization with cyclic boundary conditions or results of previous optimization problem (previous optimization horizon)
                            Xdhw[t,s] == AdhwS[1,1,s].*XdhwStartS[s,1] +
                                         BdhwS[1,1,s].*Qtank[t,s] +
                                         sum(EdhwS[1,di,s]*DisturbDHWS[di,t,s] for di=1:nbDisturbDhw))
            end
            @constraint(integratedModel,[t=2:nbTimesteps,s=1:nbHouses],         # Results of previous time step as starting point
                        Xdhw[t,s] == AdhwS[1,1,s].*Xdhw[t-1,s] +
                                     BdhwS[1,1,s].*Qtank[t,s] +
                                     sum(EdhwS[1,di,s]*DisturbDHWS[di,t,s] for di=1:nbDisturbDhw))
            @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                         XminS[t,s] <= Xdhw[t,s] <= XmaxS[t,s])

      # IN SUMMER
      else

           ## Power consumption (in summer: only domestic hot water production)
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       demand_var_onehouse[t,s] == Php_dhw[t,s] + Paux_dhw[t,s])
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       Php_dhw[t,s] <= P_dhw_100S[t,s])
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       Php_dhw[t,s] <= 0.5*P_sh_100S[t,s] + 0.5*P_dhw_100S[t,s])
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       Paux_dhw[t,s] <= PauxmaxvecS[t,s])

           ## Heating system model
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       QemD[t,s] == 0)
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       QemN[t,s] == 0)
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                       Qtank[t,s] == COPdhwS[t,s]*Php_dhw[t,s] + Paux_dhw[t,s])

           ## Building model
           @constraint(integratedModel,[q=1:nbStatesSH,t=1:nbTimesteps,s=1:nbHouses],
                       Tsh[q,t,s] == 0)

           ## Domestic hot water tank model
           if cyclic_here
               @constraint(integratedModel,[s=1:nbHouses],
                           Xdhw[1,s] == Xdhw[end,s])                           # Cyclic boundary conditions
           else
               @constraint(integratedModel,[t=1,s=1:nbHouses],                 # Results of initialization with cyclic boundary conditions or results of previous optimization problem (previous optimization horizon)
                           Xdhw[t,s] == AdhwS[1,1,s].*XdhwStartS[s,1] +
                                        BdhwS[1,1,s].*Qtank[t,s] +
                                        sum(EdhwS[1,di,s]*DisturbDHWS[di,t,s] for di=1:nbDisturbDhw))
           end
           @constraint(integratedModel,[t=2:nbTimesteps,s=1:nbHouses],         # Results of previous time step as starting point
                       Xdhw[t,s] == AdhwS[1,1,s].*Xdhw[t-1,s] +
                                    BdhwS[1,1,s].*Qtank[t,s] +
                                    sum(EdhwS[1,di,s]*DisturbDHWS[di,t,s] for di=1:nbDisturbDhw))
           @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                        XminS[t,s] <= Xdhw[t,s] <= XmaxS[t,s])

       end

   # IN CASE OF NO FLEXIBLE DEMAND
   else
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   demand_var_onehouse[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   Php_sh[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   Php_dhw[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   Paux_sh[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   Paux_dhw[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   QemD[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   QemN[t,s] == 0)
       @constraint(integratedModel,[q=1:nbStatesSH,t=1:nbTimesteps,s=1:nbHouses],
                   Tsh[q,t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   Qtank[t,s] == 0)
       @constraint(integratedModel,[t=1:nbTimesteps,s=1:nbHouses],
                   Xdhw[t,s] == 0)
   end

   #-------------------------------------------------------------------------------
   # Solve optimization problem
   JuMP.optimize!(integratedModel)
   status = JuMP.termination_status(integratedModel)

   if status == MOI.LOCALLY_SOLVED
           println("     ",status)
   else
           error("ERROR: Integrated model for day $day_start was not solved to optimality")
   end

   #-------------------------------------------------------------------------------
   # Save optimization results
   ## Only for the currently considered optimization horizon
   results_index = 1:mainDict["OptDays"]*24
   results = deepcopy(results_format[results_index,:])                          # set back to zero

   results[results_index,:cost_tot] = JuMP.objective_value(integratedModel)*ones(mainDict["OptDays"]*24)

   results[results_index,:gen_ccgt] = collect(JuMP.value.(gen_ccgt))[results_index]
   results[results_index,:gen_ocgt] = collect(JuMP.value.(gen_ocgt))[results_index]
   results[results_index,:gen_solar] = gen_solar[results_index]
   results[results_index,:gen_wind] = gen_wind[results_index]
   results[results_index,:curtail] = collect(JuMP.value.(curtail))[results_index]
   results[results_index,:load_shed] = collect(JuMP.value.(load_shed))[results_index]
   results[results_index,:dem_fxd] = demand_fxd[results_index]
   results[results_index,:dem_var] = collect(JuMP.value.(demand_var))[results_index]

   results[results_index,:cost_ccgt] = MC_ccgt.*collect(JuMP.value.(gen_ccgt))[results_index]
   results[results_index,:cost_ocgt] = MC_ocgt.*collect(JuMP.value.(gen_ocgt))[results_index]
   results[results_index,:cost_LL] = VOLL.*collect(JuMP.value.(load_shed))[results_index]

   results[results_index,:dem_var_h1] = collect(JuMP.value.(demand_var_onehouse[:,1]))[results_index]
   results[results_index,:dem_var_h2] = collect(JuMP.value.(demand_var_onehouse[:,2]))[results_index]
   results[results_index,:Php_sh_h1] = collect(JuMP.value.(Php_sh[:,1]))[results_index]
   results[results_index,:Php_sh_h2] = collect(JuMP.value.(Php_sh[:,2]))[results_index]
   results[results_index,:Php_dhw_h1] = collect(JuMP.value.(Php_dhw[:,1]))[results_index]
   results[results_index,:Php_dhw_h2] = collect(JuMP.value.(Php_dhw[:,2]))[results_index]
   results[results_index,:Paux_sh_h1] = collect(JuMP.value.(Paux_sh[:,1]))[results_index]
   results[results_index,:Paux_sh_h2] = collect(JuMP.value.(Paux_sh[:,2]))[results_index]
   results[results_index,:Paux_dhw_h1] = collect(JuMP.value.(Paux_dhw[:,1]))[results_index]
   results[results_index,:Paux_dhw_h2] = collect(JuMP.value.(Paux_dhw[:,2]))[results_index]
   results[results_index,:QemD_h1] = collect(JuMP.value.(QemD[:,1]))[results_index]
   results[results_index,:QemD_h2] = collect(JuMP.value.(QemD[:,2]))[results_index]
   results[results_index,:QemN_h1] = collect(JuMP.value.(QemN[:,1]))[results_index]
   results[results_index,:QemN_h2] = collect(JuMP.value.(QemN[:,2]))[results_index]
   results[results_index,:Qtank_h1] = collect(JuMP.value.(Qtank[:,1]))[results_index]
   results[results_index,:Qtank_h2] = collect(JuMP.value.(Qtank[:,2]))[results_index]
   results[results_index,:Tday_h1] = collect(JuMP.value.(Tsh[1,:,1]))[results_index]
   results[results_index,:Tday_h2] = collect(JuMP.value.(Tsh[1,:,2]))[results_index]
   results[results_index,:Tnight_h1] = collect(JuMP.value.(Tsh[6,:,1]))[results_index]
   results[results_index,:Tnight_h2] = collect(JuMP.value.(Tsh[6,:,2]))[results_index]
   results[results_index,:Xw_h1] = collect(JuMP.value.(Xdhw[:,1]))[results_index]
   results[results_index,:Xw_h2] = collect(JuMP.value.(Xdhw[:,2]))[results_index]

   results[results_index,:Tmin_h1] = TminS[1,1:mainDict["OptDays"]*24,1]
   results[results_index,:Tmin_h2] = TminS[1,1:mainDict["OptDays"]*24,2]
   results[results_index,:Tmax_h1] = TmaxS[1,1:mainDict["OptDays"]*24,1]
   results[results_index,:Tmax_h2] = TmaxS[1,1:mainDict["OptDays"]*24,2]
   results[results_index,:Xmin_h1] = XminS[1:mainDict["OptDays"]*24,1]
   results[results_index,:Xmin_h2] = XminS[1:mainDict["OptDays"]*24,2]
   results[results_index,:Xmax_h1] = XmaxS[1:mainDict["OptDays"]*24,1]
   results[results_index,:Xmax_h2] = XmaxS[1:mainDict["OptDays"]*24,2]

   results[results_index,:dem_noSH_fxd] = demand_noSH_fxd[results_index]
   results[results_index,:dem_HP_fxd] = demand_HP_fxd[results_index]

   CSV.write(string("RESULTS/results_d",day_start,".csv"),results;delim=';',decimal='.')

   ## for the overall considered horizon (= multiple optimization horizons)
   results_tot_startIndex = day_start*24-23
   global results_tot_index = results_tot_startIndex:results_tot_startIndex+mainDict["OptDays"]*24-1
   global results_tot[results_tot_index,:] .= results

   CSV.write(string("RESULTS/results_tot.csv"),results_tot;delim=';',decimal='.')   # Constantly overwritten


   # Save starting values for next optimization
   startValues = DataFrame(TiD =    collect(JuMP.value.(Tsh[1,results_index[end],:])),
                           TwD =    collect(JuMP.value.(Tsh[2,results_index[end],:])),
                           TwiD =   collect(JuMP.value.(Tsh[3,results_index[end],:])),
                           TflD =   collect(JuMP.value.(Tsh[4,results_index[end],:])),
                           TfDN =   collect(JuMP.value.(Tsh[5,results_index[end],:])),
                           TiN =    collect(JuMP.value.(Tsh[6,results_index[end],:])),
                           TwN =    collect(JuMP.value.(Tsh[7,results_index[end],:])),
                           TwiN =   collect(JuMP.value.(Tsh[8,results_index[end],:])),
                           TfND =   collect(JuMP.value.(Tsh[9,results_index[end],:])),
                           Xdhw =   collect(JuMP.value.(Xdhw[results_index[end],:])))

   CSV.write("RESULTS/startValues.csv",startValues;delim=';',decimal='.')       # continuously overwritten

end

#_______________________________________________________________________________
## SUMMARIZING RESULTS
#_______________________________________________________________________________
println("----- RESULTS SUMMARY -----")
print("Electricity generation CCGT:  ", sum(results_tot[:,:gen_ccgt]), " MWhel")
print(" \nElectricity generation OCGT:  ", sum(results_tot[:,:gen_ocgt]), " MWhel")
print(" \nElectricity generation PV:    ", sum(results_tot[:,:gen_solar]), " MWhel")
print(" \nElectricity generation WIND:  ", sum(results_tot[:,:gen_wind]), " MWhel")
print(" \nRES curtailment:              ", sum(results_tot[:,:curtail]), " MWhel")
print(" \nLoad shedding:                ", sum(results_tot[:,:load_shed]), " MWhel")
print(" \nNG consumption :              ", sum(results_tot[:,:gen_ccgt])/mainDict["effCCGT"] + sum(results_tot[:,:gen_ocgt])/mainDict["effOCGT"], " MWhprim")
println("")

#_______________________________________________________________________________
## CREATING PLOTS
#_______________________________________________________________________________

println("----- CREATING PLOTS -----")

if mainDict["DoPlots"]

    mkpath(string(current_dir,"/RESULTS/PLOTS"))

    plot_index = 1:mainDict["OptDays"]*24
    plot_tot_index = 1:mainDict["fullyearEffective"]
    Plots.scalefontsizes(0.9)

    for day_start in simulTimeVec

        results_index = day_start*24-23:day_start*24-23+mainDict["OptDays"]*24-1
        results = results_tot[results_index,:]

        ## Market clearing output
        #     Supply side
        x = plot_index
        data = [results.gen_ccgt results.gen_ocgt (results.gen_solar.+results.gen_wind) results.load_shed]
        p1 = areaplot(x,data,
                      label=["CCGT" "OCGT" "RES" "LL"],
                      legend=:outerbottomright,
                      seriescolor=[:steelblue :green :orange :purple],
                      linecolor=[:steelblue :green :orange :purple],
                      xlabel="Time \n [h]",
                      ylabel="Electric power \n [MW]",
                      title=" \n Week started at day $day_start \n($(Dates.monthday(Date(2013,1,1)+Dates.Day(day_start-1))[2]) $(Dates.monthname(Dates.monthday(Date(2013,1,1)+Dates.Day(day_start-1))[1])))",
                      show=false,
                      fmt=:svg)
                      # size=(1000,700),xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,titlefontsize=12,legendfontvalign=:center,legendfonthalign=:left,top_margin=10Plots.mm)

        #     Demand side
        x = plot_index
        data = [(results.dem_noSH_fxd+results.dem_HP_fxd+results.dem_var) (results.dem_noSH_fxd+results.dem_HP_fxd) results.dem_noSH_fxd]
        p1 = plot!(x,data,
                   label=["w HP all" "w HP no DR      " "w/o HP"],
                   color=[:red :black :darkgrey],
                   linewidth=1.5,
                   show=false,
                   fmt=:svg)
                   # size=(1000,700),xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,titlefontsize=12,legendfontvalign=:center,legendfonthalign=:left)

        ## Detailed demand side results (demand response)
        if ~(mainDict["share_HP_withDR"] == 0)
             x = plot_index
             data = [results.Tmin_h1 results.Tday_h1]
             p2 = plot(x,data,
                      label=["Lower bound      " "Actual T"],
                      legend=:outerbottomright,
                      color=[:red :red],
                      linestyle = [:dot :solid],
                      linewidth=2,
                      xlabel="Time \n [h]",
                      ylabel="Indoor temperature \n [°C]",
                      show=false,
                      fmt=:svg)
                      # size=(1000,700),xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,titlefontsize=12,legendfontvalign=:center,legendfonthalign=:left)

              x = plot_index
              data = [results.Xmin_h1 results.Xw_h1]
              p3 = plot(x,data,
                        label=["Lower bound      " "Actual T"],
                        legend=:outerbottomright,
                        color=[:red :red],
                        linestyle = [:dot :solid],
                        linewidth=2,
                        xlabel="Time \n [h]",
                        ylabel="DHW temperature \n [°C]",
                        show=false,
                        fmt=:svg)
                        # size=(1000,700),xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,titlefontsize=12,legendfontvalign=:center,legendfonthalign=:left)

             p4=plot(p1,p2,p3,
                     layout=(3,1),
                     size=(1000,800),
                     left_margin=15Plots.mm,right_margin=15Plots.mm,top_margin=1Plots.mm,bottom_margin=5Plots.mm,
                     titlefontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,
                     show=true,
                     fmt=:svg)
         else
             p4=plot(p1,
                     size=(1000,300),
                     left_margin=15Plots.mm,right_margin=15Plots.mm,top_margin=2Plots.mm,bottom_margin=10Plots.mm,
                     titlefontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,
                     show=true,
                     fmt=:svg)

        end
        savefig(p4,"RESULTS/PLOTS/fig_MarketClearing_d$day_start.pdf")

        ## Boundary conditions
        Tamb         = Float64.(timeSeries[results_index,:Tamb])
        QsolN        = timeSeries[results_index,:QsolN]
        QsolE        = timeSeries[results_index,:QsolE]
        QsolS        = timeSeries[results_index,:QsolS]
        QsolW        = timeSeries[results_index,:QsolW]
        Tground      = Float64(mainDict["Tground"])*ones(size(Tamb))
        QintDay      = mainDict[string("Qintdayzone",1)]*ones(size(Tamb))      # [W] internal heat gains due to appliances
        QintNight    = mainDict["constantHeatGainNightZone"]*ones(size(Tamb))

        x = plot_index
        p5 = plot(
                  plot(x,[Tamb Tground],
                       label=["Ambient T                        " "Ground T"],
                       legend=:outerbottomright,
                       color=[:red :black],
                       linestyle = [:solid :solid],
                       linewidth=2,
                       xlabel="Time \n [h]",
                       ylabel="Temperature \n [°C]",
                       title=" \n Week started at day $day_start \n($(Dates.monthday(Date(2013,1,1)+Dates.Day(day_start-1))[2]) $(Dates.monthname(Dates.monthday(Date(2013,1,1)+Dates.Day(day_start-1))[1])))",
                       show=false,
                       fmt=:svg),
                  plot(x,[QsolN QsolE QsolS QsolW],
                       label=["Solar heat gains N               " "Solar heat gains E" "Solar heat gains S" "Solar heat gains W"],
                       legend=:outerbottomright,
                       color=[:red :black :green :blue],
                       linestyle = [:solid :solid :solid :solid],
                       linewidth=2,
                       xlabel="Time \n [h]",
                       ylabel="Solar heat gains \n [W/m²]",
                       show=false,
                       fmt=:svg),
                  plot(x,[QintDay QintNight],
                       label=["Internal heat gains DZ           " "Internal heat gains NZ"],
                       legend=:outerbottomright,
                       color=[:red :black],
                       linestyle = [:solid :solid],
                       linewidth=2,
                       xlabel="Time \n [h]",
                       ylabel="Internal heat gains \n [W]",
                       show=false,
                       fmt=:svg),
                  link=:x,
                  layout=(3,1),
                  size=(1000,800),
                  left_margin=15Plots.mm,right_margin=15Plots.mm,top_margin=1Plots.mm,bottom_margin=5Plots.mm,
                  titlefontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,
                  show=true,
                  fmt=:svg)
          savefig(p5,"RESULTS/PLOTS/fig_BoundaryConditions_d$day_start.pdf")


    end

    if mainDict["SimFullYear"]

        ## Residual load curves
        residual_load_year = results_tot.gen_ccgt + results_tot.gen_ocgt - results_tot.curtail + results_tot.load_shed  # full electricity demand minus renewables (based on balance between supply and demand)
        residual_load_year_noDR = residual_load_year - results_tot.dem_var                                              # subtract the buildings with heat pumps and demand response
        residual_load_year_noHP = residual_load_year_noDR - results_tot.dem_HP_fxd                                      # subtract the heat pump demand without demand response

        residual_load = DataFrame(residual_load_year = sort(residual_load_year,rev=true),
                                  residual_load_year_noDR = sort(residual_load_year_noDR,rev=true),
                                  residual_load_year_noHP = sort(residual_load_year_noHP,rev=true)
                                  )
        CSV.write(string("RESULTS/residual_load_year.csv"),residual_load;delim=';',decimal='.')

        x = plot_tot_index
        data = [sort(residual_load_year,rev=true) sort(residual_load_year_noDR,rev=true) sort(residual_load_year_noHP,rev=true)]

        p6 = plot(x,data,
                  label=["with DR" "with HP no DR        " "w/o HP"],
                  legend=:outerbottomright,
                  color=[:red :black :darkgrey],
                  linewidth=2,
                  xlabel="Time \n [h]",
                  ylabel="Electric power \n [MW]",
                  title=" \n Residual load curves electricity",
                  show=false,
                  fmt=:svg)

         ## Gas demand of buildings without a heat pump
         demand_gas_noHP = mainDict["share_noHP"] * 10^(-6) * mainDict["nbuildsStudy"] * timeSeries[1:mainDict["fullyearEffective"],:Qtot1house]

         x = plot_tot_index
         data = sort(demand_gas_noHP,rev=true)

         p7 = plot(x,data,
                   label="GB        ",
                   legend=:outerbottomright,
                   color=:black,
                   linewidth=2,
                   xlabel="Time \n [h]",
                   ylabel="Thermal power \n [MW]",
                   title=" \n Residual load curves gas",
                   show=false,
                   fmt=:svg)

         p8= plot(p6,p7,layout=(2,1),
                  size=(1000,700),
                  left_margin=15Plots.mm,right_margin=15Plots.mm,top_margin=2Plots.mm,bottom_margin=5Plots.mm,
                  titlefontsize=12,xtickfontsize=12,ytickfontsize=12,xguidefontsize=12,yguidefontsize=12,legendfontsize=8,
                  show=true,
                  fmt=:svg)
         savefig(p8,"RESULTS/PLOTS/fig_ResidualLoad_fullYear.pdf")

    end
end


#_______________________________________________________________________________
println("")
print("RUNNING INTEGRATED MODEL FINISHED! \nYou can start analyzing the results.")
println("")

end
