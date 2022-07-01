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
function get_hh_support_bubble_status_vecs(info_per_household::Array{household_info,1})

    # Get households IDs
    hh_IDs = getproperty.(info_per_household, :household_ID)

    # Get support bubble status for each household
    hh_support_bubble_status = getproperty.(info_per_household, :in_support_bubble)

    # Produce lists
    hh_in_support_bubble = hh_IDs[hh_support_bubble_status.== true]
    hh_not_in_support_bubble = hh_IDs[hh_support_bubble_status.== false]

    return hh_in_support_bubble::Array{Int64,1},
            hh_not_in_support_bubble::Array{Int64,1}
end

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

        # Error check. Should not be more than two other households met.
        if length(hh_mixed_with_by_day[index_hh_ID]) > 2
            error("More than two households visited by hh $index_hh_ID: $(hh_mixed_with_by_day[index_hh_ID])")
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

    # Support bubble paris need a non-support bubble household.
    # Get lists of households in support bubble and not in support bubble
    hh_in_support_bubble,
    hh_not_in_support_bubble =  get_hh_support_bubble_status_vecs(info_per_household)

    # Shuffle order of hh_not_in_support_bubble
    shuffle!(rng,hh_not_in_support_bubble)

    # Initialise index counter for assigning non-support bubble households to a triplet
    # with two households that are in a support bubble
    hh_not_in_support_bubble_idx = 1

    # Support bubble pairs can only form christmas bubble with a non-support bubble household.
    support_bubble_allocated_to_christmas_bubble = zeros(Bool,length(hh_in_support_bubble))
    for support_bubble_hh_itr = 1:length(hh_in_support_bubble)
        # Get household ID of index househould in this loop iteration
        support_bubble_hh_ID = hh_in_support_bubble[support_bubble_hh_itr]

        # Error check
        if info_per_household[support_bubble_hh_ID].household_ID != support_bubble_hh_ID
            error("support_bubble_hh_ID, $support_bubble_hh_ID, and info_per_household[support_bubble_hh_ID].household_ID $(info_per_household[support_bubble_hh_ID].household_ID) do not match.")
        end

        # Check if already assigned to christmas bubble
        if support_bubble_allocated_to_christmas_bubble[support_bubble_hh_itr] == false

            # If yet to be assigned, take next household in list of hh_not_in_support_bubble
            third_hh_in_bubble_ID = hh_not_in_support_bubble[hh_not_in_support_bubble_idx]
            hh_not_in_support_bubble_idx += 1 # Increment index access counter

            # Update assignment status of other household in support bubble
            other_hh_in_support_bubble_ID = info_per_household[support_bubble_hh_ID].support_bubble_household_ID # Get ID of other household in support bubble
            support_bubble_allocated_to_christmas_bubble[hh_in_support_bubble .== other_hh_in_support_bubble_ID] .= true
                # Find relevant index in list of household IDs in support bubbles that matches
                # the ID of household in this support bubble

            # Update assignment status of index household
            support_bubble_allocated_to_christmas_bubble[support_bubble_hh_itr] = true

            # Do contact assignment for household triplet
            # Update hh_mixed_with_by_day
            christmas_bubble_hh_IDs = [support_bubble_hh_ID;
                                        other_hh_in_support_bubble_ID;
                                        third_hh_in_bubble_ID]
            update_visitation_schedule!(christmas_bubble_hh_IDs,
                                        hh_mixed_with_by_day)

        end
    end

    # Once all support bubbles allocated a triple, remaining households can be allocated as triples.
    remaining_hh_to_assign_to_bubble = length(hh_not_in_support_bubble) - (hh_not_in_support_bubble_idx - 1)
    n_triples = ceil(Int64,remaining_hh_to_assign_to_bubble/3)
    for hh_grp_itr = 1:n_triples
        start_idx = hh_not_in_support_bubble_idx + ((hh_grp_itr-1)*3) # Don't need to +1 as hh_not_in_support_bubble_idx one larger than actual number of non-support bubble hh allocated to christmas bubbles
        end_idx = min((hh_not_in_support_bubble_idx - 1) + (hh_grp_itr*3),length(hh_not_in_support_bubble))
        household_ID_slice = hh_not_in_support_bubble[start_idx:end_idx]
        update_visitation_schedule!(household_ID_slice,
                                    hh_mixed_with_by_day)
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

    # Support bubble paris need a non-support bubble household.
    # Get lists of households in support bubble and not in support bubble
    hh_in_support_bubble,
    hh_not_in_support_bubble =  get_hh_support_bubble_status_vecs(info_per_household)

    # Shuffle order of hh_not_in_support_bubble
    shuffle!(rng,hh_not_in_support_bubble)

    # Initialise index counter for assigning non-support bubble households to a triplet
    # with two households that are in a support bubble
    hh_not_in_support_bubble_idx = 1

    # Support bubble pairs can only form christmas bubble with a non-support bubble household.
    support_bubble_allocated_to_christmas_bubble = zeros(Bool,length(hh_in_support_bubble))
    for support_bubble_hh_itr = 1:length(hh_in_support_bubble)
        # Get household ID of index househould in this loop iteration
        support_bubble_hh_ID = hh_in_support_bubble[support_bubble_hh_itr]

        # Error check
        if info_per_household[support_bubble_hh_ID].household_ID != support_bubble_hh_ID
            error("support_bubble_hh_ID, $support_bubble_hh_ID, and info_per_household[support_bubble_hh_ID].household_ID $(info_per_household[support_bubble_hh_ID].household_ID) do not match.")
        end

        # Check if already assigned to christmas bubble
        if support_bubble_allocated_to_christmas_bubble[support_bubble_hh_itr] == false

            # If yet to be assigned, take next household in list of hh_not_in_support_bubble
            third_hh_in_bubble_ID = hh_not_in_support_bubble[hh_not_in_support_bubble_idx]
            hh_not_in_support_bubble_idx += 1 # Increment index access counter

            # Update assignment status of other household in support bubble
            other_hh_in_support_bubble_ID = info_per_household[support_bubble_hh_ID].support_bubble_household_ID # Get ID of other household in support bubble
            support_bubble_allocated_to_christmas_bubble[hh_in_support_bubble .== other_hh_in_support_bubble_ID] .= true
                # Find relevant index in list of household IDs in support bubbles that matches
                # the ID of household in this support bubble

            # Update assignment status of index household
            support_bubble_allocated_to_christmas_bubble[support_bubble_hh_itr] = true

            # Do contact assignment for household triplet
            # Update hh_mixed_with_by_day
            christmas_bubble_hh_IDs = [support_bubble_hh_ID;
                                        other_hh_in_support_bubble_ID;
                                        third_hh_in_bubble_ID]
            update_visitation_schedule!(christmas_bubble_hh_IDs,
                                        hh_mixed_with_by_day)

        end
    end

    # Once all support bubbles allocated a triple, remaining households can be allocated as triples.
    remaining_hh_to_assign_to_bubble = length(hh_not_in_support_bubble) - (hh_not_in_support_bubble_idx - 1)
    n_triples = ceil(Int64,remaining_hh_to_assign_to_bubble/3)
    for hh_grp_itr = 1:n_triples
        start_idx = hh_not_in_support_bubble_idx + ((hh_grp_itr-1)*3) # Don't need to +1 as hh_not_in_support_bubble_idx one larger than actual number of non-support bubble hh allocated to christmas bubbles
        end_idx = min((hh_not_in_support_bubble_idx - 1) + (hh_grp_itr*3),length(hh_not_in_support_bubble))
        household_ID_slice = hh_not_in_support_bubble[start_idx:end_idx]
        update_visitation_schedule!(household_ID_slice,
                                    hh_mixed_with_by_day)
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

    # Get lists of households in support bubble and not in support bubble
    hh_in_support_bubble,
    hh_not_in_support_bubble =  get_hh_support_bubble_status_vecs(info_per_household)

    # Get list of all households
    n_hh = length(info_per_household)
    hh_ID_vec = collect(1:n_hh)

    # Get number of households in & not in support bubbles
    n_hh_in_support_bubble = length(hh_in_support_bubble)
    n_hh_not_in_support_bubble = length(hh_not_in_support_bubble)

    # Get visit schedule entries for just households not in support bubbles
    hh_not_in_support_bubble_mixed_with_by_day = hh_mixed_with_by_day[hh_not_in_support_bubble]

    #=============================
    Bubble formation structure

    # Iterate over each household

    # If in support bubble,
    #   - Check if triplet is already assigned.
    #   - If not:
    #       - sample one household that is not in a support bubble
    #       - That completes the triplet

    # If not in a support bubble
      - Check if triplet is already assigned.
      - If not:
        - If have one spot spare, sample a household from hh_not_in_support_bubble
           that also has a spare slot
        - If has two spots spare, sample a household from hh_not_in_support_bubble
           that also has a spare slot. And then repeat.
    =============================#

    # Intialise vectors for allocation
    n_hh_assigned_to_xmas_bubbles = zeros(Int64,n_hh)
    n_hh_assigned_to_xmas_bubbles_hh_not_in_support_bubble = zeros(Int64,n_hh_not_in_support_bubble)

    # Initialise BitVectors for allocation
    bitvector_hh_with_xmas_bubble_space = BitVector(undef,n_hh)
    bitvector_hh_not_in_support_bubble_with_two_xmas_bubble_spaces = BitVector(undef,n_hh_not_in_support_bubble)
    bitvector_hh_not_in_support_bubble_with_one_xmas_bubble_space = BitVector(undef,n_hh_not_in_support_bubble)

    # Iterate over each household
    for first_hh_in_bubble_ID = 1:n_hh

        # Populate temporary arrays to be used to get eligibile households for bubble selection
        n_hh_assigned_to_xmas_bubbles_hh_not_in_support_bubble .= length.(hh_not_in_support_bubble_mixed_with_by_day)
        bitvector_hh_not_in_support_bubble_with_two_xmas_bubble_spaces .= n_hh_assigned_to_xmas_bubbles_hh_not_in_support_bubble .== 0
        bitvector_hh_not_in_support_bubble_with_one_xmas_bubble_space .= n_hh_assigned_to_xmas_bubbles_hh_not_in_support_bubble .== 1

        n_hh_assigned_to_xmas_bubbles .= length.(hh_mixed_with_by_day)
        bitvector_hh_with_xmas_bubble_space .= n_hh_assigned_to_xmas_bubbles .< 2

        # Update eligible households to be selected.
        hh_not_in_support_bubble_with_current_size_zero_xmas_bubble = hh_not_in_support_bubble[bitvector_hh_not_in_support_bubble_with_two_xmas_bubble_spaces]
        hh_not_in_support_bubble_with_current_size_one_xmas_bubble = hh_not_in_support_bubble[bitvector_hh_not_in_support_bubble_with_one_xmas_bubble_space]
        hh_with_space_in_xmas_bubble = hh_ID_vec[bitvector_hh_with_xmas_bubble_space]

        # If in support bubble, check if triplet is already assigned
        if info_per_household[first_hh_in_bubble_ID].in_support_bubble == true
            # If not, sample one household that is not in a support bubble
            if length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) == 0
                # Get ID of other household in support bubble
                second_hh_in_bubble_ID = info_per_household[first_hh_in_bubble_ID].support_bubble_household_ID

                # Sample a household from hh_not_in_support_bubble
                #   - If no spare slots, no more assignments are made.
                if length(hh_not_in_support_bubble_with_current_size_zero_xmas_bubble) > 0
                    third_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_zero_xmas_bubble)

                    # Error if household does not have two spare slots
                    if length(hh_mixed_with_by_day[third_hh_in_bubble_ID]) > 0
                        error("Invalid entry in hh_not_in_support_bubble_with_current_size_zero_xmas_bubble.")
                        #third_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_zero_xmas_bubble)
                    end

                    # Do contact assignment for household triplet
                    # Update hh_mixed_with_by_day
                    christmas_bubble_hh_IDs = [first_hh_in_bubble_ID;
                                                second_hh_in_bubble_ID;
                                                third_hh_in_bubble_ID]
                    update_visitation_schedule!(christmas_bubble_hh_IDs,
                                                hh_mixed_with_by_day)
                end
            elseif length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) != 2
                # Error check. Not possible for support bubble household to only be allocated
                # one other household to meet
                error("Support bubble household (ID $first_hh_in_bubble_ID) not previously assigned a full triplet: $(hh_mixed_with_by_day[first_hh_in_bubble_ID])")
            end
                # If length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) == 2 then no assignments needed
        else
            # Error check
            if length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) > 2
                error("More than two households visited by hh $first_hh_in_bubble_ID: $(hh_mixed_with_by_day[first_hh_in_bubble_ID])")
            end

            # Check if triplet is already assigned.
            # - If not, it will enter loop
            if length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) == 1
                #   - If have one spot spare, sample a household from hh_not_in_support_bubble
                #      that also has a spare slot
                #   - If no spare slots, no more assignments are made.
                if (length(hh_not_in_support_bubble_with_current_size_one_xmas_bubble) > 1) ||
                        ( (length(hh_not_in_support_bubble_with_current_size_one_xmas_bubble) == 1) &&
                           (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[1] != first_hh_in_bubble_ID) )
                    second_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_one_xmas_bubble)

                    # Resample if drawn value is same as first_hh_in_bubble_ID
                    while (second_hh_in_bubble_ID == first_hh_in_bubble_ID)
                        second_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_one_xmas_bubble)
                    end

                    # Error if household does not have spare slots
                    if length(hh_mixed_with_by_day[second_hh_in_bubble_ID]) == 2
                        error("Invalid entry in hh_not_in_support_bubble_with_current_size_one_xmas_bubble.")
                        #second_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_one_xmas_bubble)
                    end

                    # Update visitation schedules for both first_hh_in_bubble_ID & second_hh_in_bubble_ID
                    append!(hh_mixed_with_by_day[first_hh_in_bubble_ID],second_hh_in_bubble_ID)
                    append!(hh_mixed_with_by_day[second_hh_in_bubble_ID],first_hh_in_bubble_ID)
                end
            elseif length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) == 0

                # First sample a household
                #   - If no spare slots, no more assignments are made.
                if (length(hh_with_space_in_xmas_bubble) > 1) ||
                    ( (length(hh_with_space_in_xmas_bubble) == 1) &&
                       (hh_with_space_in_xmas_bubble[1] != first_hh_in_bubble_ID) )

                    # Sample a households
                    second_hh_in_bubble_ID = rand(rng,hh_with_space_in_xmas_bubble)

                    # Resample if drawn value is same as first_hh_in_bubble_ID
                    while (second_hh_in_bubble_ID == first_hh_in_bubble_ID)
                        second_hh_in_bubble_ID = rand(rng,hh_with_space_in_xmas_bubble)
                    end

                    # Error if household does not have spare slots
                    if length(hh_mixed_with_by_day[second_hh_in_bubble_ID]) == 2
                        error("Invalid entry in hh_with_space_in_xmas_bubble.")
                        #second_hh_in_bubble_ID = rand(rng,hh_with_space_in_xmas_bubble)
                    end

                    #   - If sampled household in a support bubble, that will complete the triplet
                    if info_per_household[second_hh_in_bubble_ID].in_support_bubble == true

                        # Error check. As household is in support bubble, it should either have fully assigned
                        # triplet or not have any household meets yet assigned
                        if length(hh_mixed_with_by_day[second_hh_in_bubble_ID]) == 1
                            error("Support bubble household (ID $second_hh_in_bubble_ID) not previously assigned a full triplet: $(hh_mixed_with_by_day[second_hh_in_bubble_ID])")
                        end

                        # Get ID of other household in support bubble
                        third_hh_in_bubble_ID = info_per_household[second_hh_in_bubble_ID].support_bubble_household_ID

                        # Do contact assignment for household triplet
                        # Update hh_mixed_with_by_day
                        christmas_bubble_hh_IDs = [first_hh_in_bubble_ID;
                                                    second_hh_in_bubble_ID;
                                                    third_hh_in_bubble_ID]
                        update_visitation_schedule!(christmas_bubble_hh_IDs,
                                                    hh_mixed_with_by_day)
                    else

                        #   - If sampled household not from a suppot bubble, assign pair of households as vitising each other
                        append!(hh_mixed_with_by_day[first_hh_in_bubble_ID],second_hh_in_bubble_ID)
                        append!(hh_mixed_with_by_day[second_hh_in_bubble_ID],first_hh_in_bubble_ID)

                        # # Update hh_not_in_support_bubble_with_current_size_one_xmas_bubble
                        # hh_not_in_support_bubble_with_current_size_one_xmas_bubble = hh_not_in_support_bubble[length.(hh_mixed_with_by_day[hh_not_in_support_bubble]) .== 1]

                        # Populate temporary arrays to be used to get eligibile households for bubble selection
                        n_hh_assigned_to_xmas_bubbles_hh_not_in_support_bubble .= length.(hh_not_in_support_bubble_mixed_with_by_day)
                        bitvector_hh_not_in_support_bubble_with_one_xmas_bubble_space .= n_hh_assigned_to_xmas_bubbles_hh_not_in_support_bubble .== 1
                        hh_not_in_support_bubble_with_current_size_one_xmas_bubble = hh_not_in_support_bubble[bitvector_hh_not_in_support_bubble_with_one_xmas_bubble_space]

                        #   - Then sample second household from hh_not_in_support_bubble that also has a spare slot.
                        #   - If no spare slots, no more assignments are made.
                        if (length(hh_not_in_support_bubble_with_current_size_one_xmas_bubble) > 2) ||  # If more than two entries, is an eligible household to be selected
                               ( (length(hh_not_in_support_bubble_with_current_size_one_xmas_bubble) == 1) && # If only one eligible household in list, check it has not already been selected
                                   (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[1] != first_hh_in_bubble_ID) &&
                                   (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[1] != second_hh_in_bubble_ID)) ||
                                ( (length(hh_not_in_support_bubble_with_current_size_one_xmas_bubble) == 2) && # If two eligible households in list, check both have not been selected
                                    ( ( (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[1] != first_hh_in_bubble_ID) && (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[2] != first_hh_in_bubble_ID) ) ||
                                      ( (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[1] != second_hh_in_bubble_ID) && (hh_not_in_support_bubble_with_current_size_one_xmas_bubble[2] != second_hh_in_bubble_ID) ) ))

                            # Sample a household to be third household on non-exclusive bubble
                            # for index household
                            third_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_one_xmas_bubble)

                            # Resample if drawn value is same as either first_hh_in_bubble_ID
                            # or second__hh_in_bubble_ID
                            while (third_hh_in_bubble_ID == first_hh_in_bubble_ID) ||
                                    (third_hh_in_bubble_ID == second_hh_in_bubble_ID)
                                third_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_one_xmas_bubble)
                            end


                            # Error if household does not have spare slots
                            if length(hh_mixed_with_by_day[third_hh_in_bubble_ID]) == 2
                                error("Invalid entry in hh_not_in_support_bubble_with_current_size_one_xmas_bubble.")
                                #third_hh_in_bubble_ID = rand(rng,hh_not_in_support_bubble_with_current_size_one_xmas_bubble)
                            end

                            # Once sampled, assign pair of households as visiting each other
                            append!(hh_mixed_with_by_day[first_hh_in_bubble_ID],third_hh_in_bubble_ID)
                            append!(hh_mixed_with_by_day[third_hh_in_bubble_ID],first_hh_in_bubble_ID)
                        end
                    end
                end
            end
        end

        # Error check at end of assignment
        if length(hh_mixed_with_by_day[first_hh_in_bubble_ID]) > 2
            error("Too many households assigned to bubble for hh $first_hh_in_bubble_ID: $(hh_mixed_with_by_day[first_hh_in_bubble_ID])")
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

    # Support bubble paris need a non-support bubble household.
    # Get lists of households in support bubble and not in support bubble
    hh_in_support_bubble,
    hh_not_in_support_bubble =  get_hh_support_bubble_status_vecs(info_per_household)

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

            # Shuffle order of hh_not_in_support_bubble
            shuffle!(rng,hh_not_in_support_bubble)

            # Initialise index counter for assigning non-support bubble households to a triplet
            # with two households that are in a support bubble
            hh_not_in_support_bubble_idx = 1

            # Support bubble pairs can only form christmas bubble with a non-support bubble household.
            support_bubble_allocated_to_christmas_bubble = zeros(Bool,length(hh_in_support_bubble))
            for support_bubble_hh_itr = 1:length(hh_in_support_bubble)
                # Get household ID of index househould in this loop iteration
                support_bubble_hh_ID = hh_in_support_bubble[support_bubble_hh_itr]

                # Error check
                if info_per_household[support_bubble_hh_ID].household_ID != support_bubble_hh_ID
                    error("support_bubble_hh_ID, $support_bubble_hh_ID, and info_per_household[support_bubble_hh_ID].household_ID $(info_per_household[support_bubble_hh_ID].household_ID) do not match.")
                end

                # Check if already assigned to christmas bubble
                if support_bubble_allocated_to_christmas_bubble[support_bubble_hh_itr] == false

                    # If yet to be assigned, take next household in list of hh_not_in_support_bubble
                    third_hh_in_bubble_ID = hh_not_in_support_bubble[hh_not_in_support_bubble_idx]
                    hh_not_in_support_bubble_idx += 1 # Increment index access counter

                    # Update assignment status of other household in support bubble
                    other_hh_in_support_bubble_ID = info_per_household[support_bubble_hh_ID].support_bubble_household_ID # Get ID of other household in support bubble
                    support_bubble_allocated_to_christmas_bubble[hh_in_support_bubble .== other_hh_in_support_bubble_ID] .= true
                        # Find relevant index in list of household IDs in support bubbles that matches
                        # the ID of household in this support bubble

                    # Update assignment status of index household
                    support_bubble_allocated_to_christmas_bubble[support_bubble_hh_itr] = true

                    # Do contact assignment for household triplet
                    # Update hh_mixed_with_by_day
                    christmas_bubble_hh_IDs = [support_bubble_hh_ID;
                                                other_hh_in_support_bubble_ID;
                                                third_hh_in_bubble_ID]
                    update_visitation_schedule!(christmas_bubble_hh_IDs,
                                                hh_mixed_with_by_day)

                end
            end

            # Once all support bubbles allocated a triple, remaining households can be allocated as triples.
            remaining_hh_to_assign_to_bubble = length(hh_not_in_support_bubble) - (hh_not_in_support_bubble_idx - 1)
            n_triples = ceil(Int64,remaining_hh_to_assign_to_bubble/3)
            for hh_grp_itr = 1:n_triples
                start_idx = hh_not_in_support_bubble_idx + ((hh_grp_itr-1)*3) # Don't need to +1 as hh_not_in_support_bubble_idx one larger than actual number of non-support bubble hh allocated to christmas bubbles
                end_idx = min((hh_not_in_support_bubble_idx - 1) + (hh_grp_itr*3),length(hh_not_in_support_bubble))
                household_ID_slice = hh_not_in_support_bubble[start_idx:end_idx]
                update_visitation_schedule!(household_ID_slice,
                                            hh_mixed_with_by_day)
            end

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
