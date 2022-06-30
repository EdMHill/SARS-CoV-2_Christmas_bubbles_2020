#=============================
Purpose:
Supporting functions for running the stochastic IBM household bubble toy model

Functions to generate households & person info:
- sample_households (Sample the requested number of households)
- generate_persons (Assign required number of people to each household)

Functions to increment state tracking variables
- increment_counters! (Used each timestep to increase the value of time in state variables)

Transmission functions contained within this file include:
- transmit_over!  (check infection to possible contacts)

Functions to reinitialise states at start of each run
- reinitialise_node_states! (multiply time series vectors by 0)
- reinitialise_CT_vars! (reset contact tracing associated variables)

Miscellaneous functions contained within this file include:
- draw_sample_from_pmf!
- set_infection_related_times! (set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence)
==============================#

#=============================
Household generation functions
==============================#

# Sample the requested number of households
function sample_households(n_households::Int64,
                            support_bubble_flag::Bool,
                            support_bubble_propn::Float64,
                            rng::MersenneTwister)
# Inputs:
# n_households::Int64 - Number of households in the system
# support_bubble_flag::Bool - Specify if household composition caters for support bubbles
# support_bubble_propn::Float64 - Proportion of support bubble eligible households that are in a support bubble
# rng::MersenneTwister - The random numer generator

# Outputs:
# household_composition::Array{Int64,2} - Per household, the number of people aged 0-19yrs, 20-64yrs, 65+yrs

    # Load the required data
    household_compositon_data = readdlm("../data/Census2011_household_size_age_strat.csv",',',Int64)
    household_composition_vals = household_compositon_data[:,1:3]
    household_dist = household_compositon_data[:,4]

    # Initialise the output household_composition array
    # Row per household. Column per age bracket (0-19yrs, 20-64yrs, 65+yrs)
    household_composition_idxs = zeros(Int64,n_households)

    # Create household info paramter type for each household to store information
    info_per_household = Array{household_info,1}(undef,n_households)

    # If applicable, specify the household composition indexes that are eligible to make a support bubble
    if support_bubble_flag == true
        elig_support_bubble_flag::BitArray{1} = (sum(household_composition_vals[:,2:3],dims=2) .== 1)[:] # Ensure output is a vector
        elig_support_bubble_flag[1] = true # Include household with single 0-19yrs individual

        # Get the vector indexes corresponding to compositions eligible to make a support bubble
        c = collect(1:length(elig_support_bubble_flag))
        elig_support_bubble_idxs::Array{Int64,1} = c[elig_support_bubble_flag]

        # Initialise indicator variable for household that are elgibile & want to form a support bubble
        # Default value is to be NOT eligible
        households_elig_and_forms_support_bubble = zeros(Int64,n_households)

        # Initialise indicator variable for household being selectable as part of a support bubble
        # Default value is YES
        selectable_household_for_support_bubble = ones(Int64,n_households)
    end

    # Sample from the available household compositions
    csum_household_dist = cumsum(household_dist)/sum(household_dist)
    println("csum_household_dist[1:10]: $(csum_household_dist[1:10])")
    println("csum_household_dist[1:10]: $(csum_household_dist[1:10])")
    for household_itr = 1:n_households
        household_composition_idx = draw_sample_from_pmf(csum_household_dist,rng)
        household_composition_idxs[household_itr] = household_composition_idx
        info_per_household[household_itr] = household_info(household_ID = household_itr,
                        household_age_composition = household_composition_vals[household_composition_idx,:])

        # If applicable, update support bubble eligibility status
        if support_bubble_flag == true
            if household_composition_idx ∈ elig_support_bubble_idxs
                # For those allowed to form a support bubble, determine if they do or not
                if rand(rng) < support_bubble_propn
                    households_elig_and_forms_support_bubble[household_itr] = true
                    selectable_household_for_support_bubble[household_itr] = true
                else
                    households_elig_and_forms_support_bubble[household_itr] = false
                    selectable_household_for_support_bubble[household_itr] = false
                end
            end
            # Note, for households not eligible to form a support bubble, we assume
            # they are all available to be selected to be in a support bubble
        end
    end

    # If support bubbles permitted, check for households that have single individual or
    # a single adult with children
    # If support bubbles active, assign support bubbles if option to use one is chosen
    if (support_bubble_flag == true) && (sum(households_elig_and_forms_support_bubble) > 0)
        # Get indexes of households eligible to be in a support bubble
        temp_vec = collect(1:length(households_elig_and_forms_support_bubble))
        household_forms_support_bubble_idxs::Array{Int64,1} = temp_vec[households_elig_and_forms_support_bubble.==true]

        # Shuffle order
        shuffle!(rng,household_forms_support_bubble_idxs)

        # Iterate over those households & check if assigned to support bubble
        n_support_bubble_households = length(household_forms_support_bubble_idxs)
        for support_bubble_itr = 1:n_support_bubble_households
            # Get relevant household index
            selected_household_ID::Int64 = household_forms_support_bubble_idxs[support_bubble_itr]

            # Check not already been allocated to a support bubble
            if info_per_household[selected_household_ID].in_support_bubble == false
                # If not, choose a selectable household
                reselect_flag = true
                bubble_household_ID = 0 # Initialise bubble_household_ID
                while reselect_flag == true
                    # Draw a household
                    bubble_household_ID = rand(rng,1:n_households)

                    # Check not already allocated to a support bubble
                    # AND it is a selectable household (recall some individual households are not)
                    # AND bubble_household_ID is not the same as selected_household_ID
                    if (info_per_household[bubble_household_ID].in_support_bubble == false) &&
                        (selectable_household_for_support_bubble[bubble_household_ID] == true) &&
                        (bubble_household_ID != selected_household_ID)
                        # If conditions are satisfied, leave the loop
                        reselect_flag = false
                    end
                end

                # Update the suport bubble statuses of the two households
                info_per_household[selected_household_ID].in_support_bubble = true
                info_per_household[selected_household_ID].support_bubble_household_ID = bubble_household_ID
                info_per_household[bubble_household_ID].in_support_bubble = true
                info_per_household[bubble_household_ID].support_bubble_household_ID = selected_household_ID
            end
        end
    end

    # Return the sampled household composition
    return info_per_household::Array{household_info,1}
end

#=============================
Assign individuals to households
=============================#
# Use household compositions to generate required number of individuals per household
# Also assign household transmission risks to each person
function generate_persons!(info_per_household::Array{household_info,1},
                            inf_params::infection_params,
                            rng::MersenneTwister)
# Inputs:
# info_per_household::Array{household_info,1} - Structure containing parameters associatd with each household
# inf_parameters::infection_params - Parameter structure for infection params
# rng::MersenneTwister - The random numer generator

# Outputs:
# info_per_person::Array{Int64,2} - Per person, a parameter structure containing info on relevant attributes
# household_contacts_per_person::Array{Int64,1} - Per person, number of others in same home

    # Get number of households in use
    n_households = length(info_per_household)

    # First need to get total number of individuals in the system.
    # Can then initialise an array of type person_info of the required size
    total_n_inds = 0
    for household_itr = 1:n_households
        age_composition_data = info_per_household[household_itr].household_age_composition
        total_n_inds += sum(age_composition_data) # Increment the cumulative number of individuals
    end

    # Initialise a parameter structure for each individual
    # & for number of household contacts
    info_per_person = Array{person_info,1}(undef,total_n_inds)
    household_contacts_per_person = zeros(Int64,total_n_inds)

    # Set up distributions for sampling transmission risk
    mean_val = inf_params.transrisk_household_group_mean
    sd_val = inf_params.transrisk_household_group_sd
    household_trans_risk_dist = Normal.(mean_val,sd_val)

    # Initialise person_ID
    person_ID::Int64 = 0

    # Iterate over each household.
    # Check age composition. Iterate over each category. Generate required number of individuals
    for household_itr = 1:n_households
        # Get age composition info
        age_composition_data = info_per_household[household_itr].household_age_composition
        n_people_in_household = sum(age_composition_data)

        # Initialise the vector to store IDs of people in the household
        info_per_household[household_itr].person_IDs_in_household = zeros(Int64,n_people_in_household)
        person_in_household_idx = 0

        for age_group_itr = 1:length(age_composition_data)
            n_ind_in_age_group = age_composition_data[age_group_itr]
            for person_itr = 1:n_ind_in_age_group
                # Increment vector access indexes
                person_ID += 1
                person_in_household_idx += 1

                # Create the person_info parameter type
                info_per_person[person_ID] = person_info(age_group = age_group_itr,
                                                            household_ID = household_itr)

                # Assign number of household contacts to vector
                household_contacts_per_person[person_ID] = n_people_in_household - 1

                # Sample household transmission risk for the individual
                # Secondary attack rate within households of various sizes were that individual the index case
                # Household sizes: 2, 3, 4, 5+
                info_per_person[person_ID].transrisk_household = rand.(rng,household_trans_risk_dist)

                # Populate IDs of household contacts
                info_per_household[household_itr].person_IDs_in_household[person_in_household_idx] = person_ID
            end
        end
    end

    # Return the peron information array
    return info_per_person::Array{person_info,1},
            household_contacts_per_person::Array{Int64,1}
end

#=============================
Functions to increment state tracking variables
=============================#
# Used each timestep to increase the value of time in state variables
function increment_counters!(states::indiv_states,
    household_isoltime::Int64,
    symp_isoltime::Int64,
    asymp_isoltime::Int64,
    contact_tracing_active::Bool;
    timeisol_CTcause::Array{Int64,1}=zeros(Int64,1),
    CT_caused_isol_limit::Int64=0)

    @unpack timelat, timeinf, timesymp, timeisol, symp_timeisol, asymp_timeisol, timeisol_CTcause,
            lattime, inftime, symptime, asymp, timeisol, symp_timeisol = states

    # Increments time counters
    n_persons = length(timelat)
    for person_itr = 1:n_persons
        if timelat[person_itr]>0
            timelat[person_itr] += 1
        end

        if timeinf[person_itr]>0
            timeinf[person_itr] += 1
        end

        if timesymp[person_itr]>0
            timesymp[person_itr] += 1
        end

        if timeisol[person_itr]>0
            timeisol[person_itr] += 1
        end

        if timeisol[person_itr]>household_isoltime
            timeisol[person_itr] = 0
        end

        # Isolation due to symptomatic infection
        if symp_timeisol[person_itr]>0
            symp_timeisol[person_itr] += 1
        end

        # Exit isolation at end of symptomatic infection
        # Also resets other isolation counters to zero
        if symp_timeisol[person_itr]>symp_isoltime
            symp_timeisol[person_itr] = 0
            timeisol[person_itr] = 0
            if contact_tracing_active == true
                timeisol_CTcause[person_itr] = 0
            end
        end

        # Isolation due to discovering asymptomatic infection
        if asymp_timeisol[person_itr]>0
            asymp_timeisol[person_itr] += 1
        end

        if asymp_timeisol[person_itr]>asymp_isoltime
            asymp_timeisol[person_itr] = 0
            timeisol[person_itr] = 0
            if contact_tracing_active == true
                timeisol_CTcause[person_itr] = 0
            end
        end

        if contact_tracing_active == true

            # Increment time in self-isolation, caused by Contact Tracing
            if timeisol_CTcause[person_itr]>0
                timeisol_CTcause[person_itr] += 1
            end

            # Reset contact tracing counter if limit is exceeded
            if timeisol_CTcause[person_itr]>CT_caused_isol_limit
                timeisol_CTcause[person_itr] = 0
            end
        end
    end
    return nothing
end

#=============================
Transmission related fns
=============================#

function transmit_over!(info_per_person::Array{person_info,1},
                    infec_params::infection_params,
                    transmission_risk::Float64,
                    infected_by::Array{Int64,1},
                    output::sim_outputs,
                    states::indiv_states,
                    probasymp_adult::Float64,
                    probasymp_children::Float64,
                    rng::MersenneTwister,
                    time::Int64,
                    count::Int64;
                    infecting_by::Int64=0,
                    contacts_to_check::Array{Int64,1}=Array{Int64,1}(),
                    inisol::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    attendance_record::Array{Int64,2} = Array{Int64,2}(undef,0,0))
    # Inputs:
    # info_per_person - Fields with individual level information
    # infec_params::infection_params - Variables associated with the transmission dynamics
    # transmission_risk::Float64 - Probability of transmission given contact made
    # infected_by::Array{Int64,1} - Record of ID of infectors of each node
    # output::sim_outputs - record of outputs from the simulations
    # states::indiv_states - status of each node
    # probasymp_adult::Float64,probasymp_children::Float64 - Asymptomatic probability
    # rng::MersenneTwister - The random number generator
    # time::Int64 - Current timestep
    # count::Int64 - Replicate ID
    # infecting_by::Int64 - Node ID of node transmitting infection
    # contacts_to_check::Int64 - list of nodes IDs to receive infection
    # inisol::Array{Int64,2} - Flag if nodes are in isolation
    # attendance_record::Array{Int64,2} - Flag if individuals are at visiting others

    # Iterate over each contact
    for contact_itr = 1:length(contacts_to_check)
        infecting_to = contacts_to_check[contact_itr]

        # Check if infecting_to is susceptible, not isolating and at work
        # only checks isolation or at work if given the arguments inisol / attendance_record
        if (states.timelat[infecting_to]==0) &&
            (isassigned(inisol) == false || inisol[time,infecting_to] == false) &&
            (isassigned(attendance_record) == false || attendance_record[time,infecting_to] == true)

            # Get susceptibility of contact
            age_grp_of_contact = info_per_person[infecting_to].age_group
            suscep_of_contact = infec_params.suscep_by_age_grp[age_grp_of_contact]

            # Get infection probability
            infec_prob = min(transmission_risk*suscep_of_contact,1)

            # Bernoulli trial to see if infection occurs
            if rand(rng) < infec_prob
                states.timelat[infecting_to] = 1
                infected_by[infecting_to] = infecting_by
                output.num_infected[count][infecting_by] += 1
                states.acquired_infection[infecting_to] = time

                # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                if states.acquired_infection[infecting_by]>0
                    # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                    output.Rt[(states.acquired_infection[infecting_by]+1),count] += 1
                end

                # if this was from an initial infection, modify the generation time
                # this sum will be divided by the total number of secondary initial infections
                if infected_by[infecting_by]==-1
                    output.mean_init_generation_time[count] += time + states.lattime[infecting_by]
                end

                # Check if infection will be asymptomatic
                # Dependent on age group the individual is in
                age_group_of_new_infected = info_per_person[infecting_to].age_group
                if age_group_of_new_infected == 1
                    probasymp = probasymp_children
                else
                    probasymp = probasymp_adult
                end
                if rand(rng) < probasymp
                    states.asymp[infecting_to] = 1
                end

                # Update latent event counter
                output.numlat[time+1,count] += 1
                output.newlat[time+1,count] += 1

                # Update age group latent event counter
                output.numlat_by_age_grp[time+1,count,age_group_of_new_infected] += 1
                output.newlat_by_age_grp[time+1,count,age_group_of_new_infected] += 1
            end
        end
    end
end

# Function to obtain transmission/SAR within group based on size of gathering
function get_transmission_risk(n_people_in_gathering::Int64,
                                transmission_risk_vec::Array{Float64,1})
# Inputs:
# n_people_in_gathering::Int64 - Group size for gathering
# transmission_risk_vec::Array{Float64,1} - Based on group size, the SAR

# Outputs:
# transtemp_gathering::Float64 - Output the relevant entry from transmission risk vector based on household size

    # Assign transmission risk based on the group size
    if n_people_in_gathering == 1
        transtemp_gathering = 0.
    elseif n_people_in_gathering == 2
        transtemp_gathering = transmission_risk_vec[1]
    elseif n_people_in_gathering == 3
        transtemp_gathering = transmission_risk_vec[2]
    elseif n_people_in_gathering == 4
        transtemp_gathering = transmission_risk_vec[3]
    else
        transtemp_gathering = transmission_risk_vec[4]
    end

    return transtemp_gathering::Float64
end

#=============================
Functions to reinitialise states at start of each run
=============================#
# Epi & isolation variables reset
function reinitialise_node_states!(states::indiv_states)
    lmul!(0,states.timelat)
    lmul!(0,states.timeinf)
    lmul!(0,states.timesymp)
    lmul!(0,states.asymp)
    lmul!(0,states.lattime)
    lmul!(0,states.timeisol)
    lmul!(0,states.symp_timeisol)
    lmul!(0,states.hh_isolation)
    lmul!(0,states.delay_adherence)
    lmul!(0,states.acquired_infection)
end

# Contact tracing variables
function reinitialise_CT_vars!(CT_vars::contact_tracing_vars,n_persons::Int64, rng::MersenneTwister,
    CT_parameters::CT_params, delay_adherence::Array{Int64,1},
    csum_test_result_delay::Array{Float64,1},max_test_result_delay::Int64)

@unpack CT_days_before_symptom_included, CT_engagement = CT_parameters

    # Reset vector tracking symptomatic cases (positive confirmed or untested)
    lmul!(0,CT_vars.Symp_cases_per_household_pos_or_unknown)

    # Variables for waiting for test results
    lmul!(0,CT_vars.Time_to_test_result)
    CT_vars.Time_to_test_result .-= 1 # Reset so all values are -1

    # Repopulate Boolean vector stating whether a false negative test result would be returned
    # and the number of days relevant for contact tracing
    lmul!(0,CT_vars.relevant_prev_days_for_CT)
    for ii = 1:n_persons

      # For each worker, initialise CT_vars.Test_result_false_negative as false
      CT_vars.Test_result_false_negative[ii] = false

      # For each worker, check if they engage with contact tracing
      engage_with_CT_rand = rand(rng)
      if engage_with_CT_rand < CT_engagement # engage with contact tracing
          CT_vars.Engage_with_CT[ii] = true
      else # do not engage with contact tracing
          CT_vars.Engage_with_CT[ii] = false
      end

      # Get amount of days to be looked back over
      # Have upper bound of 7 days post symp
      # if we put in reporting delay, needs to be above this
      CT_vars.relevant_prev_days_for_CT[ii] = min(CT_days_before_symptom_included + delay_adherence[ii],
                                          CT_days_before_symptom_included + 7)
    end

    # Repopulate time until test result received for each individual
    lmul!(0,CT_vars.CT_delay_until_test_result)
    for student_itr = 1:n_persons
                      CT_vars.CT_delay_until_test_result[student_itr] = draw_sample_from_pmf(csum_test_result_delay,
                                                                                rng;
                                                                                idx_offset = 1)
    end

    # Set up vector of vectors for storing IDs of those to be contacted in CT
    CT_vars.Inds_to_be_contacted = Array{Array{Int64,1},1}(undef,n_persons)

    # Initialise array to keep track of whether an infected recalls their infector
    lmul!(0,CT_vars.Recall_infector)
end

#=============================
Misc. fns
=============================#

# Sample from a probabiliy mass function
function draw_sample_from_pmf(csum_pmf::Array{Float64,1},
                                rng::MersenneTwister;
                                idx_offset::Int64 = 0)
# Inputs:
# csum_pmf::Array{Float64,1} - Cumulative summed probability mass function. Used to draw value from.
# rng::MersenneTwister - The random number generator
# idx_offset::Int64 = 0 - Links bin index to the quantity value

# Outputs:
# val_to_update::Int64 - Entry sampled value will be assigned to

    # Get number of elements in the pmf
    n_bins = length(csum_pmf)

    # Initialise output value
    val_to_update = 0

    # Draw random number
    # Set delay in adherence/symptoms becoming known to household
    # Find interval random number resides in
    r = rand(rng)
    allocated_flag = false # Intialise allocation flag. Switch to true when allocation done.
    bin_idx = 1   # Current interval being checked
    while (allocated_flag == false)
        if r <= csum_pmf[bin_idx]
            # Assign selected value
            # Subtract idx_offset
            val_to_update = bin_idx - idx_offset

            # Update allocation flag
            allocated_flag = true
        else
            # r does not reside in this interval. Update bin index.
            bin_idx += 1

            # Error check, if not assigned value after checked final bin value
            if bin_idx > n_bins
                error("bin_idx is now $bin_idx. The pmf only has $n_bins bins. Terminating programme.")
            end
        end
    end

    return val_to_update::Int64
end


function set_infection_related_times!(states::indiv_states,
    info_per_household::Array{household_info,1},info_per_person::Array{person_info,1},
    isolation::Int64,adherence::Float64,csum_delay_adherence::Array{Float64,1},
    d_incub::Distribution,n_persons::Int64,rng::MersenneTwister)

    time_to_symps = ceil.(rand(rng,d_incub,n_persons)) # time to symptoms
    # (for asymptomatics, the same from a silent start of "symptoms")

    # iterate over nodes to set lattime and hh_isolation
    for person_itr = 1:n_persons
        # lattime is the time from infection to infectiousness
        if time_to_symps[person_itr]-states.inftime<1
            states.lattime[person_itr] = 1  # Infectiousness can begin the day after becoming infected
        else
            states.lattime[person_itr] = time_to_symps[person_itr]-states.inftime
        end

        # Base isolation on whether the household will isolate
        if isolation==1
            # Get household ID the individual is in
            # Check intention to isolate status
            hh_ID = info_per_person[person_itr].household_ID
            isol_intention = info_per_household[hh_ID].hh_will_isolate_flag
            if isol_intention == true
                states.hh_isolation[person_itr] = 1 # adherence to household isolation = 1 if adherent, 0 if not.
            end

            # p1 = rand(rng)
            # if p1 < adherence # those who adhere will isolate when they get symptoms
            #     states.hh_isolation[person_itr] = 1 # adherence to household isolation = 1 if adherent, 0 if not.
            # end

            # Draw random number
            # Set delay in adherence/symptoms becoming known to household
            # Find interval random number resides in
            states.delay_adherence[person_itr] = draw_sample_from_pmf(csum_delay_adherence,
                                                                    rng;
                                                                    idx_offset = 1)
        end

    end
end
