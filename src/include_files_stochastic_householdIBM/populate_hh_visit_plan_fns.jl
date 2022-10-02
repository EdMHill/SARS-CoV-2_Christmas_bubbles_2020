#=============================
Purpose:
Functions to populate household visit plan structure.
Gives household IDs contacted in each visit. Stratified per person and per timestep
=============================#

# Generic function outline

# Inputs:
# hh_visit_plan::Array{Array{Array{Int64,1},1},2}
#     - Per person and per timestep, a vector of vectors (with each vector being a list of household IDs contacted)
# endtime::Int64 - Runtime of the simulation
# pre_christmas_bubble_time::Int64 - Number of timesteps to be run before christmas bubbles begin
# n_persons::Int64 - Number of individuals in the system
# info_per_person::Array{person_info,1}
# info_per_household::Array{household_info,1}

# Outputs: Directly modify the visit log


#=============================
Supporting functions
=============================#
function update_visitation_schedule!(christmas_bubble_hh_IDs::Array{Int64,1},
                                        hh_mixed_with_by_day::Array{Array{Int64,1},1})

    # Iterate over pairs of households within the household triplet.
    # Record on visitation schedule for period of time where social restrictions
    # are eased
    for index_itr = 1:length(christmas_bubble_hh_IDs)
        index_hh_ID = christmas_bubble_hh_IDs[index_itr]
        for contact_itr = 1:length(christmas_bubble_hh_IDs)
            if index_itr != contact_itr
                contact_hh_ID = christmas_bubble_hh_IDs[contact_itr]
                append!(hh_mixed_with_by_day[index_hh_ID],contact_hh_ID)
            end
        end

        # Remove any duplicate entries
        sort!(hh_mixed_with_by_day[index_hh_ID])
        unique!(hh_mixed_with_by_day[index_hh_ID])

        # Error check. If three extended households form a Christmas bubble,
        # that would give a maximum of six individual households.
        # Therefore, cannot visit more than five other households
        if length(hh_mixed_with_by_day[index_hh_ID]) > 5
            error("More than five households visited by hh $index_hh_ID: $(hh_mixed_with_by_day[index_hh_ID])")
        end
    end

    return nothing
end

function update_visitation_schedule_scenario_D!(christmas_bubble_hh_IDs::Array{Int64,1},
                                                hh_mixed_with_by_day::Array{Array{Int64,1},1})

    # Iterate over pairs of households within the household triplet.
    # Record on visitation schedule for period of time where social restrictions
    # are eased
    for index_itr = 1:length(christmas_bubble_hh_IDs)
        index_hh_ID = christmas_bubble_hh_IDs[index_itr]
        for contact_itr = 1:length(christmas_bubble_hh_IDs)
            if index_itr != contact_itr
                contact_hh_ID = christmas_bubble_hh_IDs[contact_itr]
                append!(hh_mixed_with_by_day[index_hh_ID],contact_hh_ID)
            end
        end

        # Remove any duplicate entries
        sort!(hh_mixed_with_by_day[index_hh_ID])
        unique!(hh_mixed_with_by_day[index_hh_ID])

        # Error check. If three extended households form a Christmas bubble,
        # that would give a maximum of six individual households.
        # Therefore, cannot visit more than five other households
        if length(hh_mixed_with_by_day[index_hh_ID]) > 5
            error("More than five households visited by hh $index_hh_ID: $(hh_mixed_with_by_day[index_hh_ID])")
        end
    end

    return nothing
end

function update_hh_visit_plan!(hh_visit_plan::Array{Array{Array{Int64,1},1},2},
                                    hh_mixed_with_by_day::Array{Array{Int64,1},1},
                                    endtime::Int64,
                                    n_persons::Int64,
                                    info_per_person::Array{person_info,1},
                                    info_per_household::Array{household_info,1},
                                    mixing_start_time::Int64,
                                    mixing_end_time::Int64)

    for person_itr = 1:n_persons
        hh_ID::Int64 = info_per_person[person_itr].household_ID

        for timestep_itr = 1:endtime
            if mixing_start_time <= timestep_itr <= mixing_end_time
                # A day where larger household bubbles allowed

                # If household not meeting with any others, intiailise entry in
                # hh_mixed_with_by_day
                if isassigned(hh_mixed_with_by_day,hh_ID) == false
                    println("Warning: hh_mixed_with_by_day[$hh_ID] not assigned. No more than one household should have unassigned visit schedule.")
                    hh_mixed_with_by_day[hh_ID] = []
                end

                # Update visit plan
                hh_visit_plan[timestep_itr,person_itr] = [hh_mixed_with_by_day[hh_ID]]
            else
                # A day where only meeting support bubble is allowed
                if info_per_household[hh_ID].in_support_bubble == true
                    # ID of other household in the support bubble (if index household is in one)
                    support_bubble_hh_ID = info_per_household[hh_ID].support_bubble_household_ID

                    # Set up as a vector of vectors
                    # (multiple vector for multiple different visits on the same timestep)
                    hh_visit_plan[timestep_itr,person_itr] = [[support_bubble_hh_ID]]
                else
                    # If not in a support bubble, then no other households met
                    hh_visit_plan[timestep_itr,person_itr] = [[]]
                end
            end
        end
    end

    return nothing
end

function construct_christmas_bubbles_and_update_visit_schedule!(rng::MersenneTwister,
                                                                n_households::Int64,
                                                                info_per_household::Array{household_info,1},
                                                                hh_mixed_with_by_day::Array{Array{Int64,1},1})

    # Get a list of households, with just one of pair of households included if in support bubble
    # Iterate over each household.
    # Update flag variable: Initially true. Set to false if the support bubble household of another household
    hh_to_include_in_bubble_assignment_flag_vec = ones(Bool,n_households)
    for hh_ID = 1:n_households
        if hh_to_include_in_bubble_assignment_flag_vec[hh_ID] == true

            # Check if in a support bubble
            # If true, for other household in support bubble set
            # hh_to_include_in_bubble_assignment_flag_vec to false
            if info_per_household[hh_ID].in_support_bubble == true
                support_bubble_hh_ID = info_per_household[hh_ID].support_bubble_household_ID
                hh_to_include_in_bubble_assignment_flag_vec[support_bubble_hh_ID] = false
            end
        end
    end

    # Construct lookup vector of houehold IDs to use for constructing household
    # triplets
    hh_ID_vec = collect(1:1:n_households)
    hh_bubble_assignment_IDs = hh_ID_vec[hh_to_include_in_bubble_assignment_flag_vec .== true]

    # Shuffle order of hh_bubble_assignment_IDs
    shuffle!(rng,hh_bubble_assignment_IDs)

    # Contact assignment for household triplets
    n_hh_bubble_triplets = ceil(Int64,length(hh_bubble_assignment_IDs)/3)
    for hh_bubble_itr = 1:n_hh_bubble_triplets
        start_access_idx = ((hh_bubble_itr-1)*3) + 1
        end_access_idx = (hh_bubble_itr*3)

        # Have check if end index exceeds length of hh_bubble_assignment_IDs
        if end_access_idx > length(hh_bubble_assignment_IDs)
            end_access_idx = length(hh_bubble_assignment_IDs)
        end

        # Construct list of households in Christmas bubble
        christmas_bubble_hh_IDs = Int64[]
        for hh_bubble_assignment_ID_itr = start_access_idx:1:end_access_idx
            current_itr_hh_ID = hh_bubble_assignment_IDs[hh_bubble_assignment_ID_itr]
            append!(christmas_bubble_hh_IDs,current_itr_hh_ID)

            # Check if household is in a support bubble
            # If true, add support bubble household to christmas_bubble_hh_IDs
            if info_per_household[current_itr_hh_ID].in_support_bubble == true
                support_bubble_hh_ID = info_per_household[current_itr_hh_ID].support_bubble_household_ID

                # Error check
                if info_per_household[support_bubble_hh_ID].household_ID != support_bubble_hh_ID
                    error("support_bubble_hh_ID, $support_bubble_hh_ID, and info_per_household[support_bubble_hh_ID].household_ID $(info_per_household[support_bubble_hh_ID].household_ID) do not match.")
                end

                # Add support bubble household ID to list of households in Christmas bubble
                append!(christmas_bubble_hh_IDs,support_bubble_hh_ID)
            end
        end

        # Update hh_mixed_with_by_day
        update_visitation_schedule!(christmas_bubble_hh_IDs,
                                    hh_mixed_with_by_day)
    end

end

#=============================
Main bubble formation functions
=============================#
function support_bubble_visit_plan!(hh_visit_plan::Array{Array{Array{Int64,1},1},2},
                                        endtime::Int64,
                                        pre_christmas_bubble_time::Int64,
                                        n_persons::Int64,
                                        info_per_person::Array{person_info,1},
                                        info_per_household::Array{household_info,1},
                                        rng::MersenneTwister)

    # Check if in support bubble.
    # Populate each timestep with support bubble household ID
    for person_itr = 1:n_persons
        hh_ID::Int64 = info_per_person[person_itr].household_ID

        if info_per_household[hh_ID].in_support_bubble == true
            # ID of other household in the support bubble (if index household is in one)
            support_bubble_hh_ID::Int64 = info_per_household[hh_ID].support_bubble_household_ID

            for timestep_itr = 1:endtime
                # Set up as a vector of vectors
                # (multiple vector for multiple different visits on the same timestep)
                hh_visit_plan[timestep_itr,person_itr] = [[support_bubble_hh_ID]]
            end
        else
            for timestep_itr = 1:endtime
                hh_visit_plan[timestep_itr,person_itr] = [[]]
            end
        end
    end

    return nothing
end

function three_household_faithful_bubble_shorter_visit_plan!(hh_visit_plan::Array{Array{Array{Int64,1},1},2},
                                        endtime::Int64,
                                        pre_christmas_bubble_time::Int64,
                                        n_persons::Int64,
                                        info_per_person::Array{person_info,1},
                                        info_per_household::Array{household_info,1},
                                        rng::MersenneTwister)

    # Assign timesteps Christmas bubbles occur to variables
    mixing_start_time = 3 + pre_christmas_bubble_time
    mixing_end_time = 4 + pre_christmas_bubble_time

    #Initialise households visited per timestep for each person
    n_households = length(info_per_household)
    hh_mixed_with_by_day = Array{Array{Int64,1},1}(undef,n_households)
    for hh_itr = 1:n_households
        hh_mixed_with_by_day[hh_itr] = []
    end

    # Construct Christmas bubbles and update visit schedule
    construct_christmas_bubbles_and_update_visit_schedule!(rng,
                                                            n_households,
                                                            info_per_household,
                                                            hh_mixed_with_by_day)

    # Populate each timestep with relevant bubble household ID
    update_hh_visit_plan!(hh_visit_plan,
                            hh_mixed_with_by_day,
                            endtime,
                            n_persons,
                            info_per_person,
                            info_per_household,
                            mixing_start_time,
                            mixing_end_time)

    return nothing
end

function three_household_faithful_bubble_visit_plan!(hh_visit_plan::Array{Array{Array{Int64,1},1},2},
                                        endtime::Int64,
                                        pre_christmas_bubble_time::Int64,
                                        n_persons::Int64,
                                        info_per_person::Array{person_info,1},
                                        info_per_household::Array{household_info,1},
                                        rng::MersenneTwister)

    # Assign timesteps Christmas bubbles occur to variables
    mixing_start_time = 1 + pre_christmas_bubble_time
    mixing_end_time = 5 + pre_christmas_bubble_time

    #Initialise households visited per timestep for each person
    n_households = length(info_per_household)
    hh_mixed_with_by_day = Array{Array{Int64,1},1}(undef,n_households)
    for hh_itr = 1:n_households
        hh_mixed_with_by_day[hh_itr] = []
    end

    # Construct Christmas bubbles and update visit schedule
    construct_christmas_bubbles_and_update_visit_schedule!(rng,
                                                            n_households,
                                                            info_per_household,
                                                            hh_mixed_with_by_day)

    # Populate each timestep with relevant bubble household ID
    update_hh_visit_plan!(hh_visit_plan,
                            hh_mixed_with_by_day,
                            endtime,
                            n_persons,
                            info_per_person,
                            info_per_household,
                            mixing_start_time,
                            mixing_end_time)

    return nothing
end

function three_household_fixed_bubble_visit_plan!(hh_visit_plan::Array{Array{Array{Int64,1},1},2},
                                        endtime::Int64,
                                        pre_christmas_bubble_time::Int64,
                                        n_persons::Int64,
                                        info_per_person::Array{person_info,1},
                                        info_per_household::Array{household_info,1},
                                        rng::MersenneTwister)

    # Assign timesteps Christmas bubbles occur to variables
    mixing_start_time = 1 + pre_christmas_bubble_time
    mixing_end_time = 5 + pre_christmas_bubble_time

    #Initialise households visited per timestep for each person
    n_households = length(info_per_household)
    hh_mixed_with_by_day = Array{Array{Int64,1},1}(undef,n_households)
    for hh_itr = 1:n_households
        hh_mixed_with_by_day[hh_itr] = []
    end

    # Get list of all households
    n_hh = length(info_per_household)
    hh_ID_vec = collect(1:n_hh)

    #=============================
    Bubble formation structure

    # Iterate over each household
       - Check if triplet is already assigned.
       - If not:
            - If have one spot spare, sample a household from hh_bubble_assignment_IDs
            that also has a spare slot
            - If has two spots spare, sample a household from hh_bubble_assignment_IDs
            that also has a spare slot. And then repeat.
    =============================#

    # Have a spare slot counter. Begins at 2, will be incremented down
    # in bubble assignment loop
    household_bubble_spaces_available = 2*ones(Int64,n_hh)

    # Iterate over each household
    for hh_ID = 1:n_hh
        # Only enter loop if index household has spaces to be assigned
        if household_bubble_spaces_available[hh_ID] > 0

            # Get households with spare slot
            eligible_households = hh_ID_vec[household_bubble_spaces_available .> 0]

            # Remove hh_ID from eligible households
            filter!(x->x≠hh_ID,eligible_households)

            # For index household, check if in support bubble
            # Remove ID of support bubble household from vector used to sample
            # christmas bubble household IDs
            if info_per_household[hh_ID].in_support_bubble == true
                support_bubble_hh_ID = info_per_household[hh_ID].support_bubble_household_ID
                filter!(x->x≠support_bubble_hh_ID,eligible_households)
            end

            # If elgibile household > 0, sample from those eligible households
            if length(eligible_households) > 0

                sampled_hh_ID_vec = rand(rng,eligible_households,household_bubble_spaces_available[hh_ID])

                # Check. If two households sampled, check they are not in the same support bubble.
                # If true, resample second household (if enough eligible households)
                if household_bubble_spaces_available[hh_ID] == 2
                    first_proposed_bubble_hh_ID = sampled_hh_ID_vec[1]

                    if info_per_household[first_proposed_bubble_hh_ID].in_support_bubble == true
                        first_proposed_bubble_hh_support_bubble_hh_ID = info_per_household[first_proposed_bubble_hh_ID].support_bubble_household_ID
                        if length(eligible_households) > 2
                            while first_proposed_bubble_hh_support_bubble_hh_ID == sampled_hh_ID_vec[2]
                                sampled_hh_ID_vec[2] = rand(rng,eligible_households)
                            end
                        else
                            # No other households may be sampled. Restrict list to first sampled household
                            sampled_hh_ID_vec = first_proposed_bubble_hh_ID
                        end
                    end
                end

                for sampled_hh_itr = 1:length(sampled_hh_ID_vec)

                    # Initialise christmas_bubble_hh_IDs
                    christmas_bubble_hh_IDs = Int64[]

                    # Add household selected to be part of bubble to christmas_bubble_hh_IDs
                    append!(christmas_bubble_hh_IDs,hh_ID)

                    # For current loop, get ID of sampled household to join bubble
                    sampled_hh_ID = sampled_hh_ID_vec[sampled_hh_itr]

                    # Decrease spare household slot counter for index household
                    # and sampled household
                    household_bubble_spaces_available[hh_ID] -= 1
                    household_bubble_spaces_available[sampled_hh_ID] -= 1

                    # Add sampled household selected to be part of bubble to christmas_bubble_hh_IDs
                    append!(christmas_bubble_hh_IDs,sampled_hh_ID)

                    # For index household, check if in support bubble
                    if (info_per_household[hh_ID].in_support_bubble == true)
                        # If true, add support bubble households to be part of bubble to christmas_bubble_hh_IDs
                        append!(christmas_bubble_hh_IDs,support_bubble_hh_ID)

                        # Decrease spare household slot counter for support bubble household
                        household_bubble_spaces_available[support_bubble_hh_ID] -= 1
                    end

                    # Check if sampled household in support bubble
                    if info_per_household[sampled_hh_ID].in_support_bubble == true
                        sampled_hh_support_bubble_household_ID = info_per_household[sampled_hh_ID].support_bubble_household_ID

                        # Decrease spare household slot counter for support bubble household
                        household_bubble_spaces_available[sampled_hh_support_bubble_household_ID] -= 1

                        # Add household selected to be part of bubble to christmas_bubble_hh_IDs
                        append!(christmas_bubble_hh_IDs,sampled_hh_support_bubble_household_ID)
                    end

                    # Do contact assignment for household triplet
                    # Update hh_mixed_with_by_day
                    update_visitation_schedule_scenario_D!(christmas_bubble_hh_IDs,
                                                            hh_mixed_with_by_day)
                end
            end
        end
    end

    # Populate each timestep with relevant bubble household ID
    update_hh_visit_plan!(hh_visit_plan,
                            hh_mixed_with_by_day,
                            endtime,
                            n_persons,
                            info_per_person,
                            info_per_household,
                            mixing_start_time,
                            mixing_end_time)
    return nothing
end

function three_household_unfaithful_bubble_visit_plan!(hh_visit_plan::Array{Array{Array{Int64,1},1},2},
                                        endtime::Int64,
                                        pre_christmas_bubble_time::Int64,
                                        n_persons::Int64,
                                        info_per_person::Array{person_info,1},
                                        info_per_household::Array{household_info,1},
                                        rng::MersenneTwister)

    # Assign timesteps Christmas bubbles occur to variables
    mixing_start_time = 1 + pre_christmas_bubble_time
    mixing_end_time = 5 + pre_christmas_bubble_time

    # Get number of households
    n_households = length(info_per_household)

    # Within Christmas bubble period, create new household triplets each day
    # and assign to visit schedule.
    # If outside Christmas bubble period, gatherings will be with support bubble households only
    for timestep_itr = 1:endtime
        if mixing_start_time <= timestep_itr <= mixing_end_time
            # A day within the Christmas bubble period

            # Initialise household visit schedule vector of vectors
            hh_mixed_with_by_day = Array{Array{Int64,1},1}(undef,n_households)
            for hh_itr = 1:n_households
                hh_mixed_with_by_day[hh_itr] = []
            end

            # Construct Christmas bubbles and update visit schedule
            construct_christmas_bubbles_and_update_visit_schedule!(rng,
                                                                    n_households,
                                                                    info_per_household,
                                                                    hh_mixed_with_by_day)

            # Update household visit schedule for each person
            for person_itr = 1:n_persons
                hh_ID::Int64 = info_per_person[person_itr].household_ID

                # If household not meeting with any others, intiailise entry in
                # hh_mixed_with_by_day
                if isassigned(hh_mixed_with_by_day,hh_ID) == false
                    println("Warning: hh_mixed_with_by_day[$hh_ID] not assigned on timestep $timestep_itr.")
                    println("No more than one household should have unassigned visit schedule.")
                    hh_mixed_with_by_day[hh_ID] = []
                end

                # Update visit plan
                hh_visit_plan[timestep_itr,person_itr] = [hh_mixed_with_by_day[hh_ID]]
            end
        else
            # A day where only meeting support bubble is allowed
            # Iterate over each person and check if in support bubble
            for person_itr = 1:n_persons
                hh_ID::Int64 = info_per_person[person_itr].household_ID
                if info_per_household[hh_ID].in_support_bubble == true
                    # ID of other household in the support bubble (if index household is in one)
                    support_bubble_hh_ID = info_per_household[hh_ID].support_bubble_household_ID

                    # Set up as a vector of vectors
                    # (multiple vector for multiple different visits on the same timestep)
                    hh_visit_plan[timestep_itr,person_itr] = [[support_bubble_hh_ID]]
                else
                    # If not in a support bubble, then no other households met
                    hh_visit_plan[timestep_itr,person_itr] = [[]]
                end
            end
        end
    end

    return nothing
end
