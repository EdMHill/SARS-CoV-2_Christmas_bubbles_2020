#=
Purpose:
Stash functions that are used with the stochastic household IBM
for performing contact tracing
=#


"""
Perform forward CT from an identified infector
"""
# Comment break
function forwardCT_from_infector!(infector_ID::Int64,
                                CT_vars::contact_tracing_vars,
                                CT_parameters::CT_params,
                                engage_with_CT::Array{Bool,1},
                                time::Int64,
                                count::Int64,
                                num_CT::Array{Int64,2},
                                hh_isolation::Array{Int64,1},
                                timeisol_CTcause::Array{Int64,1},
                                other_hh_contact_by_timestep::Array{Array{Int64,1},2},
                                rng::MersenneTwister,
                                infector_trace_count::Array{Int64,1})

@unpack dynamic_contacts_recalled_propn, social_contacts_recalled_propn, infector_engage_with_CT_prob = CT_parameters
@unpack Inds_to_be_contacted, Test_result_false_negative = CT_vars

    # Already contact traced? If not, do it now
    if (isassigned(Inds_to_be_contacted,infector_ID)) # Checks if infector_ID reported infection.
                                                      # If not, no forward contact tracing will be done from infector
        if !(length(Inds_to_be_contacted[infector_ID])>0) # Inds_to_be_contacted[infector_ID] being empty signifies
                                                         # infector reported symptoms, but returned false neg and/or
                                                         # did not engage in contact tracing
            # check if the infector will engage with CT
            if (((engage_with_CT[infector_ID] == false) && (rand(rng)<infector_engage_with_CT_prob))
                || (engage_with_CT[infector_ID] == true)) &&
                (Test_result_false_negative[infector_ID] == false)
                        # didn't engage before but does now or already willing to engage
                        # and then checks infector has not previously tested negative

                infector_trace_count[1] += 1

                trace_node!(infector_ID,time,CT_vars,CT_parameters,other_hh_contact_by_timestep,rng)
            end
        end
    end
end

function trace_node!(person_itr::Int64,
                        time::Int64,
                        CT_vars::contact_tracing_vars,
                        CT_parameters::CT_params,
                        other_hh_contact_by_timestep::Array{Array{Int64,1},2},
                        rng)

    # Overview:
    # - Load history of contacts
    # - Look back over the relevant time window
    # - Get contacts made on those days
    # - Append to vector tracking traceable contacts

    # Get contacts that will be contacted
    for time_itr = 1:CT_vars.relevant_prev_days_for_CT[person_itr]
        if time-time_itr > 0 # Can't look back before the simulation started

            # Get previous time being checked
            time_to_check = time-time_itr

            # Get non-immedaite household contacts made on time_to_check timestep
            other_hh_contacts = other_hh_contact_by_timestep[time_to_check,person_itr]

            # If there are contacts,
            # add to vector tracking traceable contacts
            if !isempty(other_hh_contacts)
                append!(CT_vars.Inds_to_be_contacted[person_itr],other_hh_contacts)
            end
        end
    end
end
