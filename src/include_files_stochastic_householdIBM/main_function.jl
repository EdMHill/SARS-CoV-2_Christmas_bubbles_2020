#=============================
Main function
=============================#

# Run outbreak on system of households
function run_stochastic_household_IBM_fn(RNGseed::Int64,
                                          n_households::Int64,
                                          support_bubble_flag::Bool,
                                          seed_initial_states_fn::Function,
                                          populate_hh_visit_plan_fn::Function,
                                          n_simns::Int64,
                                          endtime::Int64,
                                          pre_christmas_bubble_time::Int64,
                                          contact_tracing_active::Bool,
                                          adherence_prob::Float64
                                          )

# Inputs:
# RNGseed - Sets the random number generator
# n_households::Int64 - Number of people in the system
# support_bubble_flag::Bool - Declares if support bubbles may be formed
# seed_initial_states_fn::Function - Sets amount of nodes to be seeded as latent/asymp/symp/recovered.
# populate_hh_visit_plan_fn::Function - Function to construct schedule of household visits
# n_simns::Int64, endtime::Int64 - Simulation replicates & timesteps for each replicate
# endtime::Int64 - Number of timesteps each simulation replicate to be performed for
# pre_christmas_bubble_time::Int64 - Number of timesteps to be run before christmas bubbles begin
# contact_tracing_active::Bool - Set if contact tracing is active or not (Bool type variable)
# adherence_prob::Float64 - Probability of each household adhering to test, trace and isolate guidance.

# For parameters within the parameter structures, see "include_files_stochastic_householdIBM/parametertypes.jl"

##  OUTLINE OF THE CODE STRUCTURE
#         - Unpack required variables
#         - Set the random number generator
#         - Assign households (supporting functions in stochastic_householdIBM_support_fns.jl)
#         - Generate household specific transmission (supporting functions in stochastic_householdIBM_support_fns.jl)
#         - Initialise variables
#         - Iterate over replicates
#             - Reinitialisation phase
#             - Set course of infection times (supporting functions in stochastic_householdIBM_support_fns.jl)
#             - Update output time series with initial conditions
#             - Reset contact tracing variables (supporting functions in stochastic_householdIBM_support_fns.jl)
#             - Iterate over time
#                 - Reinitialise variables at start of timestep
#                 - Assign outputs
#                 - Increment counters (supporting functions in stochastic_householdIBM_support_fns.jl)
#                 - Increment infection process (if come to the end of latent time, move to infectious time etc.). Includes household infection loop. (supporting functions in stochastic_householdIBM_support_fns.jl)
#                 - Record whether individuals are in isolation.
#                 - Transmit infections (supporting functions in stochastic_householdIBM_support_fns.jl)
#                 - Perform contact tracing (supporting functions in contact_tracing_fns.jl)
#                 - Assign prevalence & isolation outputs
#         - Return outputs

    #=============================
    Set the random number generator
    =============================#
    rng = MersenneTwister(RNGseed)

    #=============================
    Initialise global parameter structures
    =============================#
    # Initialise parameter types
    gen_household_params = generate_household_params()
    infec_params = infection_params()
    CT_parameters = CT_params()

    # Individial may not report symptoms immediately.
    csum_delay_adherence = cumsum(infec_params.delay_adherence_pmf)

    # Update adherence parameter
    infec_params.adherence = adherence_prob

    #=============================
    Unpack required variables
    =============================#
    if contact_tracing_active==true
        @unpack CT_engagement, CT_delay_until_test_result_pmf, CT_days_before_symptom_included, test_false_negative_vec,
        CT_caused_isol_limit, prob_backwards_CT,
        perform_CT_from_infector, infector_engage_with_CT_prob = CT_parameters
    end

    #=============================
    Initialise output structure
    =============================#
    # Initialise output structure
    output = sim_outputs(endtime=endtime,countfinal=n_simns)

    #=============================
    Run the outbreak
    =============================#
    # Perform the requested amount of simulations
    for simn_itr = 1:n_simns

        #=============================
        Set the RNG
        =============================#
        rng = MersenneTwister(RNGseed+simn_itr)

        #=============================
        Initialisation phase
        =============================#

        # Sample global variables for this simn replicate from distributions
        support_bubble_propn = rand(rng,gen_household_params.support_bubble_propn_dist)
        probasymp_adults = rand(rng,infec_params.probasymp_adults_dist)
        probasymp_children = rand(rng,infec_params.probasymp_children_dist)
        asymp_trans_scaling = rand(rng,infec_params.asymp_trans_scaling_dist)

        # Sample relative susceptibility for young people for this simn replicate from distribution
        young_people_relative_suscep = rand(rng,infec_params.young_people_relative_suscep_dist)
        infec_params.suscep_by_age_grp[1] = young_people_relative_suscep

        # Reset contact structure
        contacts = contacts_struct()

        #=============================
        Generate households and persons
        =============================#
        # Get specified number of households
        # Indicate if support bubbles in use
        info_per_household = sample_households(n_households,
                                                    support_bubble_flag,
                                                    support_bubble_propn,
                                                    rng)

        # Set up individuals according to generated household compositions
        info_per_person,
        household_contacts_per_node = generate_persons!(info_per_household,
                                                        infec_params,
                                                        rng)
        n_persons = length(info_per_person)

        # Initialise num_infected vector for current simulation replicate
        output.num_infected[simn_itr] = zeros(Int64,n_persons)

        # Initialise infected_by tracking vector
        output.infected_by = zeros(Int64,n_persons)

        # Construct household contact vectors
        contacts.household_contacts = Array{Array{Int64,1},1}(undef,n_persons)
        for household_itr = 1:n_households
            n_people_in_household = length(info_per_household[household_itr].person_IDs_in_household)
            for index_person_itr = 1:n_people_in_household
                # Get ID of index person
                person_ID = info_per_household[household_itr].person_IDs_in_household[index_person_itr]

                # Initialise vector to store IDs of household contacts
                # Vector of length one less than total number of people in household
                contacts.household_contacts[person_ID] = Array{Int64,1}(undef,n_people_in_household-1)

                # Iterate through other household contacts and assign their IDs
                access_contact_vec_idx = 1
                for contact_person_itr = 1:n_people_in_household
                    if index_person_itr != contact_person_itr
                        contact_ID = info_per_household[household_itr].person_IDs_in_household[contact_person_itr]
                        contacts.household_contacts[person_ID][access_contact_vec_idx] = contact_ID
                        access_contact_vec_idx += 1
                    end
                end
            end
        end

        # Initialise tracking array to store contacts each timestep with those outside immediate household
        contacts.other_hh_contact_by_timestep = Array{Array{Int64,1},2}(undef,endtime,n_persons)

        #=============================
        Determine which households would follow isolation guidance
        =============================#
        for hh_itr = 1:n_households
            if info_per_household[hh_itr].hh_will_isolate_flag == false
            # Check household not already set to isolate
            # (may have been previously assigned as part of a support bubble)
                if rand(rng) < infec_params.adherence
                    info_per_household[hh_itr].hh_will_isolate_flag = true

                    # Check for support bubble household
                    if info_per_household[hh_itr].in_support_bubble == true
                        support_bubble_hh_ID = info_per_household[hh_itr].support_bubble_household_ID
                        info_per_household[support_bubble_hh_ID].hh_will_isolate_flag = true
                    end
                end
            end
        end

        #=============================
        Generate household visit plan array
        =============================#
        # Initialise an empty vector for each timestep and each person
        contacts.hh_visit_plan = Array{Array{Array{Int64,1},1},2}(undef,endtime,n_persons)
        populate_hh_visit_plan_fn(contacts.hh_visit_plan,
                                    endtime,
                                    pre_christmas_bubble_time,
                                    n_persons,
                                    info_per_person,
                                    info_per_household,
                                    rng)

        #=============================
        Initialise individual states, contact tracing variables &
        record of isolation per timestep
        =============================#
        states = indiv_states(n_persons=n_persons)

        contacts.daily_record_inisol = zeros(Int64,endtime,n_persons)

        #=============================
        Initialise contact tracing variables
        =============================#
         # If required, set up contact tracing related variables
         CT_vars = contact_tracing_vars(n_persons=n_persons,
                                         n_households=n_households)
         if contact_tracing_active == true
             max_test_result_delay = length(CT_delay_until_test_result_pmf) # Vars to be used when allocating test delay times
             csum_test_result_delay = cumsum(CT_delay_until_test_result_pmf)
             reinitialise_CT_vars!(CT_vars, n_persons, rng, CT_parameters, states.delay_adherence,csum_test_result_delay,max_test_result_delay)
         end

        # Initialise counter for being able to identify the infector of an infectee
        recall_infector_count = 0

        # Counter for identified infector engaging in contact tracing
        # (having not participated in CT before)
        infector_trace_count = [0]

        #=============================
        Set course of infection times
        =============================#
        # set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence
        set_infection_related_times!(states,
                                    info_per_household,info_per_person,
                                    infec_params.isolation,infec_params.adherence,
                                    csum_delay_adherence,
                                    infec_params.d_incub,
                                    n_persons,
                                    rng)

        #=============================
        Seed initial infecteds
        =============================#

        # Get number of people in each age group
        # Vector entries correspond to: [0-19yrs,20-64yrs,65+yrs]
        cumul_total_per_age_group = zeros(Int64,3)
        for household_itr = 1:n_households
            cumul_total_per_age_group += info_per_household[household_itr].household_age_composition
        end
        output.n_ppl_by_age_grp[simn_itr,:] = cumul_total_per_age_group

        # Sets latent, asymptomatic, symptomatic, recovered nodes
        n_initial_latent_by_age_group::Array{Int64,1},
        n_initial_asymp_by_age_group::Array{Int64,1},
        n_initial_symp_by_age_group::Array{Int64,1},
        n_initial_recovereds_by_age_group::Array{Int64,1}, = seed_initial_states_fn(rng,
                                                        n_persons,
                                                        cumul_total_per_age_group,
                                                        simn_itr,
                                                        info_per_person,
                                                        states,
                                                        probasymp_adults,
                                                        probasymp_children,
                                                        output.infected_by,
                                                        output)

        # Initialise amount of initial infecteds for use in tracking output vector
        simn_itr_n_initial_inf = sum(n_initial_latent_by_age_group) + sum(n_initial_asymp_by_age_group) + sum(n_initial_symp_by_age_group)
        output.num_init_infected[simn_itr] = zeros(Int64,simn_itr_n_initial_inf)

        #=============================
        Update output time series with initial conditions
        =============================#
        # Update time series for latent & infecteds after assigning initial
        # infecteds
        output.numlat[1,simn_itr] = sum(n_initial_latent_by_age_group)
        output.numinf[1,simn_itr] = sum(n_initial_asymp_by_age_group) + sum(n_initial_symp_by_age_group)
        output.numlat_by_age_grp[1,simn_itr,:] = n_initial_latent_by_age_group
        output.numinf_by_age_grp[1,simn_itr,:] = n_initial_asymp_by_age_group + n_initial_symp_by_age_group

        # Update prevalences
        output.prevlat[1,simn_itr] = sum(n_initial_latent_by_age_group)
        output.prevasymp[1,simn_itr] = sum(n_initial_asymp_by_age_group)
        output.prevpresymp[1,simn_itr] = sum(n_initial_symp_by_age_group)
        output.prevrec[1,simn_itr] = sum(n_initial_recovereds_by_age_group)

        output.prevlat_by_age_grp[1,simn_itr,:] = n_initial_latent_by_age_group
        output.prevasymp_by_age_grp[1,simn_itr,:] = n_initial_asymp_by_age_group
        output.prevpresymp_by_age_grp[1,simn_itr,:] = n_initial_symp_by_age_group
        output.prevrec_by_age_grp[1,simn_itr,:] = n_initial_recovereds_by_age_group

        #=============================
        Run single replicate
        =============================#
        for time=1:endtime

            # Initial timepoint is for initial conditions
            # Set row to access in output arrays for this timestep
            output_time_idx = time + 1

            #=============================
            Reinitialise variables at start of timestep
            =============================#
            # Reinitialise timestep specific values
            lmul!(0,states.rep_inf_this_timestep)

            #=============================
            Assign outputs
            =============================#
            # Assign counts in each disease state to array
            output.numlat[output_time_idx,simn_itr] = output.numlat[output_time_idx-1,simn_itr]
            output.numinf[output_time_idx,simn_itr] = output.numinf[output_time_idx-1,simn_itr]
            output.numrep[output_time_idx,simn_itr] = output.numrep[output_time_idx-1,simn_itr]

            for age_grp_itr = 1:3
                output.numlat_by_age_grp[output_time_idx,simn_itr,age_grp_itr] = output.numlat_by_age_grp[output_time_idx-1,simn_itr,age_grp_itr]
                output.numinf_by_age_grp[output_time_idx,simn_itr,age_grp_itr] = output.numinf_by_age_grp[output_time_idx-1,simn_itr,age_grp_itr]
                output.numrep_by_age_grp[output_time_idx,simn_itr,age_grp_itr] = output.numrep_by_age_grp[output_time_idx-1,simn_itr,age_grp_itr]
            end

            #=============================
            Increment counters
            =============================#
            # Increment counters if person is currently in that state.
            # increment_counters! in stochastic_householdIBM_support_fns.jl
            if contact_tracing_active==true
                increment_counters!(states,
                                    infec_params.household_isoltime,
                                    infec_params.symp_isoltime,
                                    infec_params.asymp_isoltime,
                                    contact_tracing_active,
                                    timeisol_CTcause=states.timeisol_CTcause,
                                    CT_caused_isol_limit=CT_caused_isol_limit)
            else
                increment_counters!(states,
                                    infec_params.household_isoltime,
                                    infec_params.symp_isoltime,
                                    infec_params.asymp_isoltime,
                                    contact_tracing_active)
            end

            #=============================
            Increment infection process
            =============================#
            # If come to the end of latent time, move to infectious time etc
            for person_itr = 1:n_persons
                # if the node has reached the end of latent infection
                if states.timelat[person_itr]>states.lattime[person_itr]
                    # move to being infectious
                    states.timelat[person_itr] = -1
                    states.timeinf[person_itr] = 1

                    # Increment time series counts
                    output.numinf[output_time_idx,simn_itr] += 1
                    output.newinf[output_time_idx,simn_itr] += 1

                    # Update age group latent event counter
                    age_group_of_indiv = info_per_person[person_itr].age_group
                    output.numinf_by_age_grp[time+1,simn_itr,age_group_of_indiv] += 1
                    output.newinf_by_age_grp[time+1,simn_itr,age_group_of_indiv] += 1

                    # check if new infected will be asymptomatic
                    # Also set household transmission risk
                    if states.asymp[person_itr] > 0
                        output.newasymp[output_time_idx,simn_itr] += 1
                        output.newasymp_by_age_grp[output_time_idx,simn_itr,age_group_of_indiv] += 1
                    end
                end

                # Update node disease state time vectors
                if states.timeinf[person_itr]>states.inftime
                    # the node becomes symptomatic (if they develop symptoms)
                    states.timeinf[person_itr] = -1
                    states.timesymp[person_itr] = 1

                    # Increment time series counts
                    output.numrep[output_time_idx,simn_itr] += 1
                    output.numrep_by_age_grp[output_time_idx,simn_itr,info_per_person[person_itr].age_group] += 1

                    # Check if index case are symptomatic & would have zero adherence delay
                    if (states.asymp[person_itr] == 0) && (states.delay_adherence[person_itr]==0)
                        # Check if infected will isolate
                        if (states.hh_isolation[person_itr]==1)
                            states.symp_timeisol[person_itr] = 1

                            # Set that the unit has reported infection this timestep
                            states.rep_inf_this_timestep[person_itr] = 1

                            # Assign time of reporting to field in person info parameter type
                            info_per_person[person_itr].time_of_reporting_infection = time

                            # Household members will also isolate
                            for hh = 1:household_contacts_per_node[person_itr]
                                contact_ID = contacts.household_contacts[person_itr][hh]
                                if (states.symp_timeisol[contact_ID]==0) # Individual not already symptomatic themselves
                                    states.timeisol[contact_ID] = 1
                                end
                            end

                            # Support bubble household would also isolate
                            hh_ID = info_per_person[person_itr].household_ID
                            if info_per_household[hh_ID].in_support_bubble == true
                                support_bubble_hh_ID = info_per_household[hh_ID].support_bubble_household_ID
                                persons_in_support_bubble_hh = info_per_household[support_bubble_hh_ID].person_IDs_in_household
                                for support_bubble_persons_itr = 1:length(persons_in_support_bubble_hh)
                                    contact_ID = info_per_household[support_bubble_hh_ID].person_IDs_in_household[support_bubble_persons_itr]
                                    if (states.symp_timeisol[contact_ID]==0) # Individual not already symptomatic themselves
                                        states.timeisol[contact_ID] = 1
                                    end
                                end
                            end
                        end

                        # If contact tracing active, increase number of symptomatic infections
                        # in household by one
                        if contact_tracing_active == true
                            if (states.asymp[person_itr] == 0) # Check case is symptomatic
                                current_person_household_ID = info_per_person[person_itr].household_ID
                                CT_vars.Symp_cases_per_household_pos_or_unknown[current_person_household_ID] += 1
                            end
                        end
                    end
                end

                # Check if node, if having a delayed adherence, begins adherence on current day
                if (states.timesymp[person_itr] > 1)&&((states.timesymp[person_itr]-1)==states.delay_adherence[person_itr]) # Condition for node beginning adherence on current day & has been symptomatic for at least one day
                    if states.asymp[person_itr] == 0 # Check node is symptomatic and will adhere
                        if states.hh_isolation[person_itr]==1 # Check node will adhere
                            states.symp_timeisol[person_itr] = 1 + states.delay_adherence[person_itr]

                            # Set that the unit has reported infection this timestep
                            states.rep_inf_this_timestep[person_itr] = 1

                            # Assign time of reporting to field in student parameter type
                            info_per_person[person_itr].time_of_reporting_infection = time

                            # Household members will also now isolate
                            for hh = 1:household_contacts_per_node[person_itr]
                                contact_ID = contacts.household_contacts[person_itr][hh]
                                if (states.symp_timeisol[contact_ID]==0) # Individual not already symptomatic themselves
                                    states.timeisol[contact_ID] = 1 + states.delay_adherence[person_itr]
                                        # Individual shortens isolation by length of time since
                                        # unwell individual began displaying symptoms
                                end
                            end

                            # Support bubble household would also isolate
                            hh_ID = info_per_person[person_itr].household_ID
                            if info_per_household[hh_ID].in_support_bubble == true
                                support_bubble_hh_ID = info_per_household[hh_ID].support_bubble_household_ID
                                persons_in_support_bubble_hh = info_per_household[support_bubble_hh_ID].person_IDs_in_household
                                for support_bubble_persons_itr = 1:length(persons_in_support_bubble_hh)
                                    contact_ID = info_per_household[support_bubble_hh_ID].person_IDs_in_household[support_bubble_persons_itr]
                                    if (states.symp_timeisol[contact_ID]==0) # Individual not already symptomatic themselves
                                        states.timeisol[contact_ID] = 1
                                    end
                                end
                            end
                        end
                    end
                end

                # Check if node has reached end of symptom period
                if states.timesymp[person_itr]>states.symptime
                    states.timesymp[person_itr] = -1

                    # If contact tracing active and case was symptomatic,
                    # decrease number of symptomatic infections in household by one
                    if contact_tracing_active == true
                        # Check case is symptomatic & not returned a false negative (if false negative, has already been subtracted)
                        if (states.asymp[person_itr] == 0) && (CT_vars.Test_result_false_negative[person_itr] == false)
                            current_node_household_ID = info_per_person[person_itr].household_ID
                            CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] -= 1

                            # Error check
                            if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] < 0
                                error("CT_vars.Symp_cases_per_household_pos_or_unknown contains a negative entry. Terminate programme.")
                            end

                        end
                    end
                end
            end

            #=============================
            Record whether individuals are in isolation.
            Update household record of number of people in isolation
            =============================#
            # record whether individuals are in isolation on this timestep
            for person_itr = 1:n_persons
                # Record whether the individual is in isolation on current day
                if (infec_params.isolation>0) &&
                    ((states.timeisol[person_itr]>0) || (states.symp_timeisol[person_itr]>0) ||
                      (states.timeisol_CTcause[person_itr]>0))
                    # Default value is 0. So only update if individual is in isolation for any reason
                    contacts.daily_record_inisol[time,person_itr] = 1
                end
            end

            # Iterate over each household. Count amount in isolation/not in isolation
            for household_itr = 1:n_households
                # Reset isolation number tracking values
                info_per_household[household_itr].n_hh_members_not_in_isol = 0
                info_per_household[household_itr].n_hh_members_in_isol = 0

                # Iterate over each indiviual in the household. Record isolation status
                hh_size = length(info_per_household[household_itr].person_IDs_in_household)
                for within_hh_itr = 1:hh_size
                    within_hh_person_ID = info_per_household[household_itr].person_IDs_in_household[within_hh_itr]
                    if contacts.daily_record_inisol[time,within_hh_person_ID] == false
                        # Household member NOT in isolation
                        info_per_household[household_itr].n_hh_members_not_in_isol += 1
                    else # Household member IS in isolation
                        info_per_household[household_itr].n_hh_members_in_isol += 1
                    end
                end
            end

            #=============================
            For each person, add to record those visited from other households on
            the current timestep
            =============================#
            for person_itr = 1:n_persons
                # Initialise record of contacts vector
                contacts.other_hh_contact_by_timestep[time,person_itr] = []

                # Check if individual is in isolation on current day
                # If in isolation, no contacts outside immediate household
                # If individual not in isolation, check other housholds visited
                if contacts.daily_record_inisol[time,person_itr] == false
                    # Find individuals not in isolation. Add to contact record.
                    index_hh_ID = info_per_person[person_itr].household_ID

                    # Check daily visit record (Array of 1D vectors. Column per person, row per timestep)
                    # For each visit, check if happens (i.e. check households not isolating) & get size of gathering
                    # Check for infection against those in visit (outside immediate household)
                    if !isempty(contacts.hh_visit_plan[time,person_itr][1])
                        n_visits_today = length(contacts.hh_visit_plan[time,person_itr])
                    else # Empty vector entry, so no visits between households
                        n_visits_today = 0
                    end
                    for hh_visit_itr = 1:n_visits_today
                        # List of household IDs person scheduled to see in this visit
                        current_visit_hh_list = contacts.hh_visit_plan[time,person_itr][hh_visit_itr]

                        # Iterate over households & add those not isolating to contact tracing record
                        for current_visit_itr = 1:length(current_visit_hh_list)
                            # Get ID of household involved in the visit
                            visited_hh_ID = current_visit_hh_list[current_visit_itr]

                            # Check if individuals in visited_hh_ID are isolating
                            n_inds_to_check = length(info_per_household[visited_hh_ID].person_IDs_in_household)
                            for contact_check_itr = 1:n_inds_to_check
                                contact_to_check_ID = info_per_household[visited_hh_ID].person_IDs_in_household[contact_check_itr]
                                if contacts.daily_record_inisol[time,contact_to_check_ID] == false
                                    append!(contacts.other_hh_contact_by_timestep[time,person_itr],contact_to_check_ID)
                                end
                            end
                        end
                    end
                end
            end

            #=============================
            Transmit infections

            Structure:
            - Check infection of immediate household members
            - Check daily visit record
            - For each visit, check if happens & get size of gathering
            - Check for infection against those in visit (outside immediate household)
            =============================#
            # Iterate over nodes that may be able to transmit infection
            for person_itr = 1:n_persons
                if ((states.timeinf[person_itr]>0) | (states.timesymp[person_itr]>0))
                        # Only enter loop if node is capable of transmitting infection

                    # find the total time infectious
                    if states.timeinf[person_itr]>0
                        tot_time_inf = states.timeinf[person_itr]
                    else
                        tot_time_inf = states.timesymp[person_itr]+states.inftime
                    end

                    # find the infectiousness
                    infectiousness = infec_params.dist_infectivity[tot_time_inf]
                    current_person = info_per_person[person_itr]
                    if states.asymp[person_itr]>0 # Asymptomatic
                        transtemp_household_vec = current_person.transrisk_household.*infectiousness*asymp_trans_scaling*infec_params.scale_FOI
                    else
                        transtemp_household_vec = current_person.transrisk_household.*infectiousness*infec_params.scale_FOI
                    end

                    # Infection check for other household members
                    # Transmit over household_contacts[person_itr]
                    # checking that contacts are susceptible
                    n_hh_contacts = length(contacts.household_contacts[person_itr])
                    if n_hh_contacts > 0
                        hh_size = n_hh_contacts + 1 # Total number of people in the household
                        transtemp_household::Float64 = get_transmission_risk(hh_size,
                                                                            transtemp_household_vec)

                        transmit_over!(info_per_person,
                                                infec_params,
                                                transtemp_household,
                                                    output.infected_by,
                                                    output,
                                                    states,
                                                    probasymp_adults,
                                                    probasymp_children,
                                                    rng,
                                                    time,
                                                    simn_itr,
                                                    infecting_by = person_itr,
                                                    contacts_to_check = contacts.household_contacts[person_itr])
                    end

                    # Check if index person is isolating or not.
                    # If not, check other households that may be visited on this timestep
                    if contacts.daily_record_inisol[time,person_itr] == false
                        # Get number of individuals from index person household that are not isolating
                        # and so are also part of the gathering
                        index_hh_ID = info_per_person[person_itr].household_ID
                        index_hh_ppl_in_gathering = info_per_household[index_hh_ID].n_hh_members_not_in_isol

                        # Check daily visit record (Array of 1D vectors. Column per person, row per timestep)
                        # For each visit, check if happens (i.e. check households not isolating) & get size of gathering
                        # Check for infection against those in visit (outside immediate household)
                        if !isempty(contacts.hh_visit_plan[time,person_itr][1])
                            n_visits_today = length(contacts.hh_visit_plan[time,person_itr])
                        else # Empty vector entry, so no visits between households
                            n_visits_today = 0
                        end
                        for hh_visit_itr = 1:n_visits_today
                            # List of household IDs person scheduled to see in this visit
                            current_visit_hh_list = contacts.hh_visit_plan[time,person_itr][hh_visit_itr]

                            # Get number of individuals collectively involved in the gathering
                            # Iterate over households & get number not isolating
                            ppl_in_gathering_total = 0
                            ppl_in_gathering_total += index_hh_ppl_in_gathering
                            for current_visit_itr = 1:length(current_visit_hh_list)
                                # Get ID of household involved in the visit
                                visited_hh_ID = current_visit_hh_list[current_visit_itr]

                                # Get number of non-isolating individuals in the visiting household
                                n_ppl_not_isol_in_visited_hh = info_per_household[visited_hh_ID].n_hh_members_not_in_isol

                                # Add to the tracking total of people involved in the gathering
                                ppl_in_gathering_total += n_ppl_not_isol_in_visited_hh
                            end

                            # Set overall transmission risk (based on numbers involved in the gathering)
                            gathering_transmission_risk = get_transmission_risk(ppl_in_gathering_total,
                                                                                transtemp_household_vec)

                            # Iterate over households in current list.
                            # Check if household takes part in visit. If so, run infection check.
                            for current_visit_itr = 1:length(current_visit_hh_list)
                                # Visit check.
                                visited_hh_ID = current_visit_hh_list[current_visit_itr]

                                # Infection check.
                                transmit_over!(info_per_person,
                                                        infec_params,
                                                        gathering_transmission_risk,
                                                            output.infected_by,
                                                            output,
                                                            states,
                                                            probasymp_adults,
                                                            probasymp_children,
                                                            rng,
                                                            time,
                                                            simn_itr,
                                                            infecting_by = person_itr,
                                                            contacts_to_check = info_per_household[visited_hh_ID].person_IDs_in_household,
                                                            inisol = contacts.daily_record_inisol)
                            end
                        end
                    end
                end
            end

            #=============================
            Perform contact tracing
            =============================#
            # If in use, enact contact tracing from index cases reported today
            if contact_tracing_active == true

                # Store contacts made during day.
                # For those reporting symptoms, start delay to test result (if needed)
                for person_itr = 1:n_persons

                    # Increment time to test result if currently waiting for that to occur
                    if CT_vars.Time_to_test_result[person_itr]>=0
                        CT_vars.Time_to_test_result[person_itr] += 1
                    end

                    # For current person, check if would be leaving pre-symptomatic phase
                    # and displaying symptoms
                    # If so, and will not return a false negative test result, gather traceable contacts
                    if (states.rep_inf_this_timestep[person_itr] == 1)

                        # Initialise CT_vars.Time_to_test_result value
                        CT_vars.Time_to_test_result[person_itr] = 0

                        # Increment test counter
                        output.tests_performed[output_time_idx,simn_itr] += 1

                        # Determine whether test result will return a negative outcome
                        # - Get time since person_itr became infected
                        # - Given time since infected, look up probability case will return negative test result
                        if states.timeinf[person_itr]>0
                            tot_time_inf = states.timeinf[person_itr]
                        else
                            tot_time_inf = states.timesymp[person_itr]+states.inftime
                        end
                        test_false_negative_prob = test_false_negative_vec[tot_time_inf]

                        # Bernoulli trial to determine if false negative returned
                        if rand(rng) < test_false_negative_prob
                            CT_vars.Test_result_false_negative[person_itr] = true
                        end

                        # Check if the individual is destined to return a false negative result
                        # If false negative to be returned, do not need to work out who traceable contacts are
                        # Otherwise, gather traceable contacts
                        # Also, if no isolation in use, no need to gather contacts.
                        CT_vars.Inds_to_be_contacted[person_itr] = Int64[] # Initialise vector to store contacts
                        if (CT_vars.Test_result_false_negative[person_itr] == false) &&
                            (infec_params.isolation>0) && (CT_vars.Engage_with_CT[person_itr] == true)

                            trace_node!(person_itr,time,CT_vars,CT_parameters,contacts.other_hh_contact_by_timestep,rng)

                            # if we are doing "backward contact tracing"
                            # some small chance the infector is included in this
                            # don't try to backwards trace the initial infections
                            if ( (rand(rng)<prob_backwards_CT) && (infected_by[person_itr]!=-1) &&
                                (infected_by[person_itr]!=-2) && (infected_by[person_itr]!=-3) )
                                append!(CT_vars.Inds_to_be_contacted[person_itr],infected_by[person_itr])
                                CT_vars.Recall_infector[person_itr] = 1
                            end
                        end

                        # Remove duplicates in CT_vars.Inds_to_be_contacted[person_itr]
                        unique!(CT_vars.Inds_to_be_contacted[person_itr])
                    end

                    # Check if delay to test result reached
                    if CT_vars.Time_to_test_result[person_itr] >= CT_vars.CT_delay_until_test_result[person_itr]

                        # Reset the CT_vars.Time_to_test_result counter
                        CT_vars.Time_to_test_result[person_itr] = -1

                        # If delay time passed, check if test returned a false negative.
                        if CT_vars.Test_result_false_negative[person_itr] == true

                            # Increment false negative counter
                            output.test_outcomes[output_time_idx,simn_itr,2] += 1

                            # Release index case from symptomatic caused isolation
                            # If yes, they cannot be released
                            # If no, they can be released.
                            states.symp_timeisol[person_itr] = 0

                            # Amend tracker of symptomatic cases, unknown test result
                            current_node_household_ID = info_per_person[person_itr].household_ID
                            CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] -= 1

                            # Error check
                            if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] < 0
                                error("CT_vars.Symp_cases_per_household_pos_or_unknown contains a negative entry. Terminate programme.")
                            end

                            # If returning a false negative, can release other household members from isolation
                            # Only release if no other household members are symptomatic
                            if CT_vars.Symp_cases_per_household_pos_or_unknown[current_node_household_ID] == 0
                                n_household_contacts = length(contacts.household_contacts[person_itr])
                                for household_contact_itr = 1:n_household_contacts
                                    household_contact_ID = contacts.household_contacts[person_itr][household_contact_itr]
                                    if states.timeisol_CTcause[household_contact_ID] == 0 # Check not also self isolating due to contact tracing
                                        states.timeisol[household_contact_ID] = 0
                                    end
                                end
                            end

                        else

                            # Increment true positive counter
                            output.test_outcomes[output_time_idx,simn_itr,1] += 1

                            # If test is positive, contacts told to isolate
                            # Get number of recallable contacts
                            n_recallable_contacts = length(CT_vars.Inds_to_be_contacted[person_itr])
                            output.num_CT[output_time_idx,simn_itr] += n_recallable_contacts

                            # Do contacts isolate or not?
                            # Isolation based on adherence to household isolation of that individual
                            if infec_params.isolation>0
                                for recallable_contact_idx = 1:n_recallable_contacts
                                    recallable_contact_ID = CT_vars.Inds_to_be_contacted[person_itr][recallable_contact_idx]

                                    # Check if individual will adhere to guidance
                                    # If so, they self-isolate
                                    if states.hh_isolation[recallable_contact_ID] == 1
                                        # Time needed to spend in isolation reduced by test result delay time for index case
                                        states.timeisol_CTcause[recallable_contact_ID] = 1 + CT_vars.CT_delay_until_test_result[person_itr]
                                    end
                                end

                                # Perform forwards CT from infector, if infector has been identified
                                # and such a policy is active
                                if perform_CT_from_infector == true
                                    if CT_vars.Recall_infector[person_itr] == 1
                                        recall_infector_count += 1

                                        forwardCT_from_infector!(infected_by[person_itr],
                                        CT_vars,
                                        CT_parameters,
                                        CT_vars.Engage_with_CT,
                                        time,
                                        simn_itr,
                                        output.num_CT,
                                        states.hh_isolation,
                                        states.timeisol_CTcause,
                                        contacts.other_hh_contact_by_timestep,
                                        rng,
                                        infector_trace_count)
                                    end
                                end
                            end
                        end
                    end
                end
            end


            #=============================
            Assign prevalence & isolation outputs
            =============================#
            # For this timestep, get number isolating
            # and if they are latent infected or infectious on current timestep
            for person_itr=1:n_persons

                # Get age group individual is in
                age_grp_idx = info_per_person[person_itr].age_group

                # Isolating due to housemate symptoms
                isolating_for_any_reason = false
                if states.timeisol[person_itr]>0
                    output.num_household_isolating[output_time_idx,simn_itr] += 1
                    isolating_for_any_reason = true
                end

                # Isolating due to symptoms
                if states.symp_timeisol[person_itr]>0
                    output.num_symp_isolating[output_time_idx,simn_itr] += 1
                    isolating_for_any_reason = true
                end

                # Isolating due to finding asymptomatic infection (due to testing)
                if states.asymp_timeisol[person_itr]>0
                    output.num_asymp_isolating[output_time_idx,simn_itr] += 1
                    isolating_for_any_reason = true
                end

                # Isolating as close contact of positive test symptomtic
                if contact_tracing_active == true
                    if states.timeisol_CTcause[person_itr]>0
                        output.num_isolating_CTcause[output_time_idx,simn_itr] += 1
                        isolating_for_any_reason = true
                    end
                end

                # Isolating for any reason
                if isolating_for_any_reason == true
                    output.num_isolating[output_time_idx,simn_itr] += 1
                end

                # Check if latently infected
                if states.timelat[person_itr]>0
                    output.prevlat[output_time_idx,simn_itr] += 1
                    output.prevlat_by_age_grp[output_time_idx,simn_itr,age_grp_idx] += 1
                end

                # In presymptomatic infectious period.
                if (states.timeinf[person_itr]>0)
                    if states.asymp[person_itr] > 0 # asymptomatic
                        output.prevasymp[output_time_idx,simn_itr] += 1
                        output.prevasymp_by_age_grp[output_time_idx,simn_itr,age_grp_idx] += 1
                    else # will be symptomatic
                        output.prevpresymp[output_time_idx,simn_itr] += 1
                        output.prevpresymp_by_age_grp[output_time_idx,simn_itr,age_grp_idx] += 1
                    end
                end

                # After presymp period, check if symptomatic or asymptomatic
                if (states.timesymp[person_itr]>0)
                    if states.asymp[person_itr] > 0 # asymptomatic
                        output.prevasymp[output_time_idx,simn_itr] += 1
                        output.prevasymp_by_age_grp[output_time_idx,simn_itr,age_grp_idx] += 1
                    else # symptomatic
                        output.prevsymp[output_time_idx,simn_itr] += 1
                        output.prevsymp_by_age_grp[output_time_idx,simn_itr,age_grp_idx] += 1
                    end
                end

                # Check if recovered
                if states.timesymp[person_itr] == -1
                    output.prevrec[output_time_idx,simn_itr] += 1
                    output.prevrec_by_age_grp[output_time_idx,simn_itr,age_grp_idx] += 1
                end
            end
        end


        #=============================
        End of simulation replicate summary stats
        =============================#
        # Find how many individuals were infected by the initial infecteds
        initial_infected_indiv = findall(output.infected_by.==-1)
        #println("initial_infected_indiv: $(length(initial_infected_indiv))")
        sum_infections = 0
        output.num_init_infected[simn_itr] = zeros(Int64,length(initial_infected_indiv)) # Initialise output vector
        for initial_infected_it = 1:length(initial_infected_indiv)
            output.num_init_infected[simn_itr][initial_infected_it] = output.num_infected[simn_itr][initial_infected_indiv[initial_infected_it]]
            sum_infections+=output.num_init_infected[simn_itr][initial_infected_it]
        end

        # Find how many nodes were set to adhere to isolation guidance for this replicate
        output.n_isol_adhering[simn_itr] = sum(states.hh_isolation)

        # find mean generation time
        if sum_infections>0
            output.mean_init_generation_time[simn_itr] = output.mean_init_generation_time[simn_itr]/sum_infections
        end

        # divide number of infections by number of infectors to find Rt
        for time=1:(endtime+1)
            # divide by the number of nodes that were infected (entered the latent state)
            # at time
            if time == 1
                output.Rt[time,simn_itr] = output.Rt[time,simn_itr] / output.numlat[time,simn_itr]
            else
                output.Rt[time,simn_itr] = output.Rt[time,simn_itr] / (output.numlat[time,simn_itr]-output.numlat[time-1,simn_itr])
            end
        end

        # Print to screen info on run just completed
        println("Run $simn_itr complete.")
    end

  # Specify what is output from the function
  return output::sim_outputs
end
