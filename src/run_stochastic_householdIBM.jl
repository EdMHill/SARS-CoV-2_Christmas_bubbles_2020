#=============================
Purpose:
Run the stochastic IBM household bubble toy model
=============================#

#=============================
Set paths & load environment
=============================#

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../")

#=============================
Load packages
=============================#
#Required packages
using Distributions
using LinearAlgebra, Random, DelimitedFiles, Parameters
using MAT
# using Profile  # For use with memory allocation checks

#=============================
Load supporting files
=============================#
# Files containing other functions needed to run the model
include("include_files_stochastic_householdIBM/main_function.jl")
include("include_files_stochastic_householdIBM/stochastic_householdIBM_parametertypes.jl")
include("include_files_stochastic_householdIBM/stochastic_householdIBM_support_fns.jl")
include("include_files_stochastic_householdIBM/seed_initial_states_fns.jl")
include("include_files_stochastic_householdIBM/populate_hh_visit_plan_fns.jl")
include("include_files_stochastic_householdIBM/contact_tracing_fns.jl")

#=============================
Set variables from ARGS
=============================#
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
if length(ARGS)==0
    args = ["1","1234","2","100000","15"]
end

# Set identifier for job
job_ID = parse(Int64, args[1])

# Use job_ID to:
#   - Choose the function to construct the household visit plan
#   - Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
if mod(job_ID,10) == 1
    # Scenario A, main analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = support_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 0
elseif mod(job_ID,10) == 2
    # Scenario B, main analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_faithful_bubble_shorter_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 0
elseif mod(job_ID,10) == 3
    # Scenario C, main analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_faithful_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 0
elseif mod(job_ID,10) == 4
    # Scenario D, main analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_fixed_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 0
elseif mod(job_ID,10) == 5
    # Scenario E, main analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_unfaithful_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 0
elseif mod(job_ID,10) == 6
    # Scenario A, alternative analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = support_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 10
elseif mod(job_ID,10) == 7
    # Scenario B, alternative analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_faithful_bubble_shorter_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 10
elseif mod(job_ID,10) == 8
    # Scenario C, alternative analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_faithful_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 10
elseif mod(job_ID,10) == 9
    # Scenario D, alternative analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_fixed_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 10
elseif mod(job_ID,10) == 0
    # Scenario E, alternative analysis

    # Set the function to construct the household visit plan
    populate_hh_visit_plan_fn = three_household_unfaithful_bubble_visit_plan!

    # Set the number of days the simulation should begin before the christmas bubble window (pre_christmas_bubble_time)
    pre_christmas_bubble_time = 10
end

# Set adherence level to be used
if job_ID <= 10
    adherence_prob = 0.7
elseif job_ID <= 20
    adherence_prob = 0.3
else
    error("No adherence level specified for the requested job_ID")
end

#=============================
Specify global variables
=============================#
# Set the RNG
RNGseed = parse(Int64, args[2])

# Specify the number of simulation replicates to be performed
n_simns = parse(Int64, args[3])

# Specify the household associated inputs
n_households = parse(Int64, args[4])
support_bubble_flag = true

# State the number of timesteps per simulation
time_duration_from_christmas_bubble_start = parse(Int64, args[5])
endtime = time_duration_from_christmas_bubble_start + pre_christmas_bubble_time

# Choose the seed initial infected fn
seed_initial_states_fn = seed_states_using_age #set_five_symp_five_asymp_initial

# Specify if contact tracing is active
contact_tracing_active = true

#=============================
Call the main function
=============================#
# Assign outputs to variable
output = run_stochastic_household_IBM_fn(RNGseed,
                                          n_households,
                                          support_bubble_flag,
                                          seed_initial_states_fn,
                                          populate_hh_visit_plan_fn,
                                          n_simns,
                                          endtime,
                                          pre_christmas_bubble_time,
                                          contact_tracing_active,
                                          adherence_prob
                                          )

#=============================
Save outputs to file
=============================#
# Set filename
if job_ID == 1
    output_filename = "../results/MAT_files/support_bubble_visit_plan_main_analysis.mat"
elseif job_ID == 2
    output_filename = "../results/MAT_files/three_household_faithful_bubble_shorter_visit_plan_main_analysis.mat"
elseif job_ID == 3
    output_filename = "../results/MAT_files/three_household_faithful_bubble_visit_plan_main_analysis.mat"
elseif job_ID == 4
    output_filename = "../results/MAT_files/three_household_fixed_bubble_visit_plan_main_analysis.mat"
elseif job_ID == 5
    output_filename = "../results/MAT_files/three_household_unfaithful_bubble_visit_plan_main_analysis.mat"
elseif job_ID == 6
    output_filename = "../results/MAT_files/support_bubble_visit_plan_alternative_analysis.mat"
elseif job_ID == 7
    output_filename = "../results/MAT_files/three_household_faithful_bubble_shorter_visit_plan_alternative_analysis.mat"
elseif job_ID == 8
    output_filename = "../results/MAT_files/three_household_faithful_bubble_visit_plan_alternative_analysis.mat"
elseif job_ID == 9
    output_filename = "../results/MAT_files/three_household_fixed_bubble_visit_plan_alternative_analysis.mat"
elseif job_ID == 10
    output_filename = "../results/MAT_files/three_household_unfaithful_bubble_visit_plan_alternative_analysis.mat"
elseif job_ID == 11
    output_filename = "../results/MAT_files/support_bubble_visit_plan_main_analysis_lower_adherence.mat"
elseif job_ID == 12
    output_filename = "../results/MAT_files/three_household_faithful_bubble_shorter_visit_plan_main_analysis_lower_adherence.mat"
elseif job_ID == 13
    output_filename = "../results/MAT_files/three_household_faithful_bubble_visit_plan_main_analysis_lower_adherence.mat"
elseif job_ID == 14
    output_filename = "../results/MAT_files/three_household_fixed_bubble_visit_plan_main_analysis_lower_adherence.mat"
elseif job_ID == 15
    output_filename = "../results/MAT_files/three_household_unfaithful_bubble_visit_plan_main_analysis_lower_adherence.mat"
elseif job_ID == 16
    output_filename = "../results/MAT_files/support_bubble_visit_plan_alternative_analysis_lower_adherence.mat"
elseif job_ID == 17
    output_filename = "../results/MAT_files/three_household_faithful_bubble_shorter_visit_plan_alternative_analysis_lower_adherence.mat"
elseif job_ID == 18
    output_filename = "../results/MAT_files/three_household_faithful_bubble_visit_plan_alternative_analysis_lower_adherence.mat"
elseif job_ID == 19
    output_filename = "../results/MAT_files/three_household_fixed_bubble_visit_plan_alternative_analysis_lower_adherence.mat"
elseif job_ID == 20
    output_filename = "../results/MAT_files/three_household_unfaithful_bubble_visit_plan_alternative_analysis_lower_adherence.mat"
else
    error("Invalid populate_hh_visit_plan_fn name provided.")
end

file = matopen(output_filename, "w")
write(file, "numlat", output.numlat)
write(file, "numinf", output.numinf)
write(file, "numrep", output.numrep)
write(file, "prevlat", output.prevlat)
write(file, "prevsymp", output.prevsymp)
write(file, "prevasymp", output.prevasymp)
write(file, "prevrec", output.prevrec)
write(file, "newlat", output.newlat)
write(file, "newinf", output.newinf)
write(file, "newasymp", output.newasymp)
write(file, "numlat_by_age_grp", output.numlat_by_age_grp)
write(file, "numinf_by_age_grp", output.numinf_by_age_grp)
write(file, "numrep_by_age_grp", output.numrep_by_age_grp)
write(file, "prevlat_by_age_grp", output.prevlat_by_age_grp)
write(file, "prevsymp_by_age_grp", output.prevsymp_by_age_grp)
write(file, "prevasymp_by_age_grp", output.prevasymp_by_age_grp)
write(file, "prevrec_by_age_grp", output.prevrec_by_age_grp)
write(file, "newlat_by_age_grp", output.newlat_by_age_grp)
write(file, "newinf_by_age_grp", output.newinf_by_age_grp)
write(file, "newasymp_by_age_grp", output.newasymp_by_age_grp)
write(file, "n_ppl_by_age_grp", output.n_ppl_by_age_grp)
write(file, "num_isolating", output.num_isolating)
write(file, "num_household_isolating", output.num_household_isolating)
write(file, "num_symp_isolating", output.num_symp_isolating)
write(file, "num_isolating_CTcause", output.num_isolating_CTcause)
write(file, "n_isol_adhering", output.n_isol_adhering)
close(file)
