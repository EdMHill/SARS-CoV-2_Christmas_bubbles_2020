#=
Purpose:
Parameter types to be used with the stochastic IBM household bubble toy model

- household_info: Household parameter type to have information on the household each individual is in
- person_info: Person parameter type to have information on each individual in the system
- generate_household_params: Parameters used in the generation of the households in the system
- indiv_states: Used for the status of each node
- infection_params: Variables associated with the epidemiological process
- contacts_struct: Used for contact structures
- CTParams: Contact tracing related
- contact_tracing_vars: Variables used for contact tracing
- sim_outputs: Outputs to be saved from the simulations
=#

# Household parameter type to have information on the household each individual is in
@with_kw mutable struct household_info

  # Household the individual has been assigned to. ID for the global system.
  household_ID::Int64 = 0

  # State the age composition of the household
  # Entry defns: [0-19yrs,20-64yrs,65+yrs]
  household_age_composition::Array{Int64,1} = [0,0,0]

  # A log of person IDs that reside in the household
  person_IDs_in_household::Array{Int64,1} = Int64[]

  # Store whether household follows isolation guidance
  hh_will_isolate_flag::Bool = false

  # Track the number of household members that are in isolation/not in isolation
  n_hh_members_in_isol::Int64 = 0
  n_hh_members_not_in_isol::Int64 = 0

  # Specify if in a support bubble
  in_support_bubble::Bool = false

  # ID of other household in the support bubble
  support_bubble_household_ID::Int64 = 0 # Empty if not support bubble household assigned
end

# Parameters used in the generation of the households in the system
@with_kw struct generate_household_params

  # Proportion of support bubble eligible households that are in a support bubble
  support_bubble_propn_dist::Distribution = Uniform(0.5,0.75)
end


# Person parameter type to have information on each individual in the system
@with_kw mutable struct person_info

  # Age category
  age_group::Int64 = 0  # Defns - 1: 0-19yrs; 2: 20-64yrs; 3: 65+yrs.

  # household assigned to
  household_ID::Int64 = 0

  # Secondary attack rate within households of various sizes were that individual the index case
  transrisk_household::Array{Float64,1} = [0.,0.,0.,0.] # Household sizes: 2, 3, 4, 5+

  # If student is symptomatically infected and it is reported, the time the case report occurs.
  # Asymptomatic infection also logged on day of positive test.
  # If never reporting infection, keeps value 0.
  time_of_reporting_infection::Int64 = 0
end

# Used for the status of each node
@with_kw mutable struct indiv_states
   n_persons::Int64 = 0
   timelat::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent in latent state
   timeinf::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent in infectious state
   timesymp::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent in symptomatic state
   lattime::Array{Int64,1} = zeros(Int64,n_persons) # time to spend in latent state
   inftime::Int64 = 4 # time to spend in infectious, pre-symptomatic state
   symptime::Int64 = 10 # time to spend in symptomatic state
   asymp::Array{Int64,1} = zeros(Int64,n_persons) # whether the node will be asymptomatic
   timeisol::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent in isolation due to housemate symptoms
   symp_timeisol::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent isolating due to symptoms
   asymp_timeisol::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent isolating due to positive result in mass testing (finding asymptomatic infection)
   timeisol_CTcause::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent isolating due to being contact traced
   visit_hh_timeisol::Array{Int64,1} = zeros(Int64,n_persons) # time currently spent isolating due to a household visited having a positive case
   rep_inf_this_timestep::Array{Int64,1} = zeros(Int64,n_persons) # whether the node reports symptoms this timestep
   hh_isolation::Array{Int64,1} = zeros(Int64,n_persons) # Whether individual adheres to isolation guidance. (1) Yes. (0) No.
   delay_adherence::Array{Int64,1} = zeros(Int64,n_persons) # Individual may not report symptoms immediately.
   acquired_infection::Array{Int64,1} = zeros(Int64,n_persons) # time node acquired infection
end

# Used for parameters relating to the infection process
@with_kw mutable struct infection_params

  # Distribution of incubation period
  d_incub::Distribution = Erlang(6,0.88)

  # Distribution of infectiousness
  dist_infectivity::Array{Float64,1} = [0.0369, 0.0491, 0.0835, 0.1190, 0.1439, 0.1497, 0.1354, 0.1076, 0.0757, 0.0476, 0.0269, 0.0138, 0.0064, 0.0044] # Corrected He (4 day pre-symptomatic period) & 10 day symptomatic period

  # probability of transmission within a household, based on secondary attack rate
  # Can differ per household group. Values for household size of 2, 3, 4 or more.
  transrisk_household_group_mean::Array{Float64,1} = [0.48,0.40,0.33,0.22]
  transrisk_household_group_sd::Array{Float64,1} = [0.06,0.06,0.05,0.05]

  # Scale the force of infection
  scale_FOI::Float64 = 1.

  # Young people age group relative susceptibility
  young_people_relative_suscep_dist::Uniform{Float64} = Uniform(0.4,0.6)

  # Initialise susceptiblity by age group vector.
  # Entries correspond to: 0-19yrs, 20-64yrs, 65+yrs
  # First entry will be updated in each simulation.
  suscep_by_age_grp = [0.5,1.0,1.0]

  # Also add in the possibility that individuals are asymptomatic.
  # In that case, their infection potential is lower but they dont isolate.
  # asymp_trans_scaling::Float64 = 0.165
  asymp_trans_scaling_dist::Uniform{Float64} = Uniform(0.3,0.7)

  # probability of being asymptomatic
  probasymp_adults_dist::Uniform{Float64} = Uniform(0.05,0.2)
  probasymp_children_dist::Uniform{Float64} = Uniform(0.2,0.35)

  # Flag variable indicating if self-isolation is an active measure
  isolation::Int64 = 1
  symp_isoltime::Int64 = 10
  asymp_isoltime::Int64 = 10
  household_isoltime::Int64 = 14

  # Set proportion of population who will adhere
  adherence::Float64 = 0.7

  # Distribution of delay in reporting symptomns (for those that do/eventually adhere)
  delay_adherence_pmf::Array{Float64,1} = [1.,0.,0.] # Note first entry is 0 day delay,
                                                        # second entry 1 day delay etc

  # environmental transmission risk
  environmental_infection_prob::Float64 = 0.
end

# Used for contact structures
@with_kw mutable struct contacts_struct

  # Per person, household contacts
  household_contacts::Array{Array{Int64,1},1} = Array{Int64,1}[]

  # Record of household visits that occur
  # Per person and per timestep, a vector of vectors (with each vector being a list of household IDs contacted)
  hh_visit_plan::Array{Array{Array{Int64,1},1},2} = Array{Array{Array{Int64,1},1},2}(undef,0,0)

  # Record of non-immediate household contacts made each day
  other_hh_contact_by_timestep::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)

  # Flag arrays to keep a log of whether an individual is in isolation each timestep
  # Used in each replicate
  daily_record_inisol::Array{Int64,2} = zeros(Int64,0,0)
end

@with_kw mutable struct CT_params
# used for parameters relating to contact tracing

   # For those adhering to self-isolating, engagement with contact tracing aspect
   # If 1, as well as self-isolating, all adhering individuals do give contacts
   CT_engagement::Float64 = 1.

   # Set time delay. 0 corresponds to result on day of reporting.
   CT_delay_until_test_result_pmf::Array{Float64,1} = [1., 0., 0.,]

   # Set number of days before symptoms CT will attempt to capture
   CT_days_before_symptom_included::Int64 = 2

   # Propotions of tests returning a false negative outcome
   # Entry per day since infected
   test_false_negative_vec::Array{Float64,1} = 0.13*ones(20)

   # Amount of time spent in isolation if contact traced
   CT_caused_isol_limit::Int64 = 14

   # proportion of people that can identify their infector (set to 0 for no backwards CT)
   prob_backwards_CT::Float64 = 0.

   # Parameters for performing forward contact tracing from identified infectors
   perform_CT_from_infector::Bool = false
   infector_engage_with_CT_prob::Float64 = 1.0  # For those not complying with other isolation guidance,
                                       # probability they do if idenfitied as possible infector
                                       # Note, those that are set to adhere to isolation guidance
                                       # are also assuemd fully compliant with CT from infector measures
end

@with_kw mutable struct contact_tracing_vars
# only needed for contact tracing

   # sizes to initialise arrays
   n_persons::Int64 = 0
   n_households::Int64 = 0

   # The number of days prior to symptoms that each node remembers
   relevant_prev_days_for_CT::Array{Int64,1} = zeros(Int64,n_persons)

   # Vector of vectors for storing IDs of those to be contacted in CT
   Inds_to_be_contacted::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,n_persons)

   # vector tracking symptomatic cases (positive confirmed or untested)
   Symp_cases_per_household_pos_or_unknown::Array{Int64,1} = zeros(Int64,n_households)

   # Variables for waiting for test results, set at -1 until activated
   Time_to_test_result::Array{Int64,1} = -1*ones(Int64,n_persons)

   # Boolean vector to store whether a false negative test result would be returned
   # Populated in each replicate
   Test_result_false_negative::Array{Bool,1} = Array{Bool,1}(undef,n_persons)

   # Boolean vector to store if individual will engage with CT or not
   # Populated in each replicate
   Engage_with_CT::Array{Bool,1} = Array{Bool,1}(undef,n_persons)

   # Array to keep track of whether an infected recalls their infector
   Recall_infector::Array{Int64,1} = zeros(Int64,n_persons)

   # Delay before test result is returned
   CT_delay_until_test_result::Array{Int64,1} = zeros(Int64,n_persons)
end

# Outputs to be saved from the simulations
@with_kw mutable struct sim_outputs
   endtime::Int64 = 0
   countfinal::Int64 = 0
   n_persons::Int64 = 0

   # 2D outputs
   numlat::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number latently infected
   numinf::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number infectious
   numrep::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number reporting infection
   prevlat::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for latently infected
   prevsymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for symptomatic infectious
   prevasymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for asymptomatic infectious
   prevpresymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # Prevalence for recovereds
   newlat::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of new latent infecteds
   newinf::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of new infectious
   newasymp::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of newly infected asymptomatics
   num_CT::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # total number of recallable contacts
   num_infected::Array{Array{Int64,1}} = Array{Array{Int64,1}}(undef,countfinal) # number of infections caused by that individual
   num_init_infected::Array{Array{Int64,1}} = Array{Array{Int64,1}}(undef,countfinal) # number of infections caused by the initially infected nodes
   Rt::Array{Float64,2} = zeros(Float64,endtime+1,countfinal) # real time R value (number of secondary infections caused by nodes that were newly infected on that day)

   # 2D outputs. Age group summary stats
   n_ppl_by_age_grp::Array{Int64,2} = zeros(Int64,countfinal,3)     # Number of people of age: [0-19yrs,20-64yrs,65+yrs]

   # 2D outputs. Isolation related
   num_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating at each timepoint
   num_household_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating due to household members having symptoms at each timepoint
   num_symp_isolating::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating due to symptoms at each timepoint
   num_isolating_CTcause::Array{Int64,2} = zeros(Int64,endtime+1,countfinal) # number of individuals isolating as a contact of a positive test

   #2D outputs. Testing related
   tests_performed::Array{Int64,2} = zeros(Int64,endtime+1,countfinal)

   # 3D outputs - age group disease status
   numlat_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # number latently infected
   numinf_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # number infectious
   numrep_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # number reporting infection
   prevlat_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # Prevalence for latently infected
   prevsymp_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # Prevalence for symptomatic infectious
   prevasymp_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # Prevalence for asymptomatic infectious
   prevpresymp_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # Prevalence for recovereds
   newlat_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # number of new latent infecteds
   newinf_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # number of new infectious
   newasymp_by_age_grp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,3) # number of newly infected asymptomatics

   #3D outputs. Testing related
    # Slices: 1 - True positive; 2 - false negative; 3 - True negative; 4 - false positive.
   test_outcomes::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,4)

   # 1D outputs - Isolation
   n_isol_adhering::Array{Int64,1} = zeros(Float64,countfinal) # Per replicate, number of individuals that are set to adhere to isolation guidance

   # 1D outputs
   infected_by::Array{Int64,1} = zeros(Int64,n_persons)  # who each node was infected by
   var_num_infected::Array{Int64,1} = zeros(Int64,countfinal) # variance in the number of infections caused by each node
   mean_init_generation_time::Array{Float64,1} = zeros(Float64,countfinal) # mean initial generation time
end
