#=
Purpose:
Store functions used to seed individuals in non-susceptible states at beginning of
simulation replicates

Fns to select individuals to begin in given non-susceptible state
- choose_from_all_popn  (select from the entire population)

Fns to get number of individuals to be seeded in each non-susceptible state
- set_five_symp_five_asymp_initial
- seed_states_with_uncertainty
- seed_states_using_age
=#


#=
Fns to select individuals to begin in given non-susceptible state
=#
# select from the entire population. Sample time already elapsed in disease state
function choose_from_all_popn(rng::MersenneTwister,
                                    relevant_person_IDs::Array{Int64,1},
                                    info_per_person::Array{person_info,1},
                                    states::indiv_states,
                                    probasymp_adults::Float64,
                                    probasymp_children::Float64,
                                    infected_by::Array{Int64,1},
                                    n_initial_latent::Int64,
                                    n_initial_asymp::Int64,
                                    n_initial_symp::Int64,
                                    n_initial_rec::Int64,
                                    initialise_start_disease_state_flag::Bool,
                                    count::Int64,
                                    output::sim_outputs)
# Inputs:
# rng::MersenneTwister - The random number generator
# relevant_person_IDs::Array{Int64,1} - Vector of IDs for those individuals that are being sampled from
# states::indiv_states - Record of info per individual
# probasymp::Float64 - Probability of infected case being asymptomatic
# infected_by::Array{Int64,1} - Record of who each individual was infected by
# n_initial_latent, n_initial_asymp, n_initial_symp, n_initial_rec
#   -  Numbers to be seeded in the named disease state
# initialise_start_disease_state_flag::Bool - If true, those assigned to non-susceptible state are initialised
#                   as having just entered that disease state (no elapsed time prior occurred)
# count::Int64 - Replicate ID
# output::sim_outputs - record of outputs from the simulations

# Ouputs:
# Directly mutates inputs states and infected_by

   # Initialise initial asymp infectious individuals
   for seed_inf_itr = 1:n_initial_asymp
      chosen_indiv = rand(rng,relevant_person_IDs)

      valid_asymp_indiv = false
      while valid_asymp_indiv == false

          if (states.timeinf[chosen_indiv] > 0) ||
              (states.timesymp[chosen_indiv] > 0) # Chosen individual already selected as an asymptomatic
             # Resample
             chosen_indiv = rand(rng,relevant_person_IDs)
          else
             # Chosen individual has asymptomatic status.
             valid_asymp_indiv = true
          end
      end

      # Set individual as an asymptomatic
      states.asymp[chosen_indiv] = 1

      # Update states for chosen individual
      states.timelat[chosen_indiv] = -1
      infected_by[chosen_indiv] = -1

      # Update time elapsed in state
      if initialise_start_disease_state_flag == true
          # Just entered infectious state
          states.timeinf[chosen_indiv] = 1
          states.acquired_infection[chosen_indiv] = -states.lattime[chosen_indiv]
      else
          # Sample time already elapsed as infected
          infected_period_length = states.inftime + states.symptime
          infected_time_elapsed::Int64 = rand(rng,1:infected_period_length)

          if infected_time_elapsed > states.inftime  # Elapsed infection time beyond pre-symptomatic stage
              states.timesymp[chosen_indiv] = infected_time_elapsed - states.inftime
              states.timeinf[chosen_indiv] = -1
          else
              states.timeinf[chosen_indiv] = infected_time_elapsed
          end
          states.acquired_infection[chosen_indiv] = -states.lattime[chosen_indiv] - (infected_time_elapsed - 1)
      end
   end

   # Initialise initial symptomatic infectious individuals
   for seed_inf_itr = 1:n_initial_symp
        chosen_indiv = rand(rng,relevant_person_IDs)

           valid_symp_individual = false
           while valid_symp_individual == false
               if (states.asymp[chosen_indiv] == 1) ||  # Chosen individual does not have symptomatic status.
                       (states.timeinf[chosen_indiv] > 0) ||
                       (states.timesymp[chosen_indiv] > 0) # Chosen individual already selected as infectious
                  # Resample
                  chosen_indiv = rand(rng,relevant_person_IDs)
               else
                  # Chosen individual has symptomatic status.
                  valid_symp_individual = true
               end
           end

           # Update states for chosen individual
           states.timelat[chosen_indiv] = -1
           infected_by[chosen_indiv] = -1

           # Update time elapsed in state
           if initialise_start_disease_state_flag == true
               # Just entered infectious state
               states.timeinf[chosen_indiv] = 1
               states.acquired_infection[chosen_indiv] = -states.lattime[chosen_indiv]
           else
               # Sample time already elapsed as infected (from pre-symptomatic period)
               presymp_time_elapsed::Int64 = rand(rng,1:states.inftime)

               states.timeinf[chosen_indiv] = presymp_time_elapsed
               states.acquired_infection[chosen_indiv] = -states.lattime[chosen_indiv] - (presymp_time_elapsed - 1)
            end
   end

   # Initialise initial latent individuals
   for seed_latent_itr = 1:n_initial_latent
      chosen_latent_indiv = rand(rng,relevant_person_IDs)

      # Check if individual already set to be an initial infected or initial latent
      valid_latent_individual = false
      while valid_latent_individual == false
             if (states.timeinf[chosen_latent_indiv] > 0) ||  # Initial infected condition
                  (states.timesymp[chosen_latent_indiv] > 0) ||
                  (states.timelat[chosen_latent_indiv] > 0)  # Initial latent condition
                  # Redraw sample if already set to be an initial infected or initial recovered
                  chosen_latent_indiv = rand(rng,relevant_person_IDs)
             else
                  # Valid individual now drawn. Update flag
                  valid_latent_individual = true
             end
      end

      # Update states variables
      infected_by[chosen_latent_indiv] = -2               # Unknown infector (someone outside student popn). Set to -2 (rather than -1) to differentiate from initial infectious individuals.

      # Update time elapsed in state
      if initialise_start_disease_state_flag == true
          # Just entered infectious state
          states.timelat[chosen_latent_indiv] = 1
          states.acquired_infection[chosen_latent_indiv] = 0
      else
          # Sample time spent latent from length of lattime
          latent_period_length = states.lattime[chosen_latent_indiv]
          latent_time_elapsed::Int64 = rand(rng,1:latent_period_length)

          # Update states variables
          states.timelat[chosen_latent_indiv] = latent_time_elapsed
          states.acquired_infection[chosen_latent_indiv] = -(latent_time_elapsed - 1)
      end

      # Check if infection will be asymptomatic
      # Determine asymp probability based on age group
      chosen_person_age_grp = info_per_person[chosen_latent_indiv].age_group
      if chosen_person_age_grp == 1
          probasymp = probasymp_children
      else
          probasymp = probasymp_adults
      end
      if rand(rng) < probasymp # Now check if case will be asymptomatic
          states.asymp[chosen_latent_indiv] = 1
      end
   end

   # Initialise individuals that are already recovered
   for seed_rec_itr = 1:n_initial_rec
       chosen_rec_indiv = rand(rng,relevant_person_IDs)

       # Check if individual already set to be an initial infected, initial latent
       # or initial recovered
       valid_individual = false
       while valid_individual == false
           if (states.timeinf[chosen_rec_indiv] > 0) ||  # Initial infected condition
                 (states.timesymp[chosen_rec_indiv] > 0) ||  # Initial infected condition
                 (states.timelat[chosen_rec_indiv] > 0) || # Initial latent condition
               ( (states.timelat[chosen_rec_indiv] == -1) && # Initial recovered conditions
                   (states.timeinf[chosen_rec_indiv] == -1) &&
                   (states.timesymp[chosen_rec_indiv] == -1) )

               # Redraw sample if already set to be an initial infected or initial recovered
               chosen_rec_indiv = rand(rng,relevant_person_IDs)
           else
               # Valid individual now drawn. Update flag
               valid_individual = true
           end
       end

       # Set infection state variables to -1 (as previously progressed through them)
       states.timelat[chosen_rec_indiv] = -1
       states.timeinf[chosen_rec_indiv] = -1
       states.timesymp[chosen_rec_indiv] = -1
       infected_by[chosen_rec_indiv] = -3      # Unknown infector (someone outside student popn).
   end

   return nothing
end


function choose_from_age_groups(rng::MersenneTwister,
                                    n_persons::Int64,
                                    n_ppl_by_age_grp::Array{Int64,1},
                                    info_per_person::Array{person_info,1},
                                    states::indiv_states,
                                    probasymp_adults::Float64,
                                    probasymp_children::Float64,
                                    infected_by::Array{Int64,1},
                                    n_initial_latent_by_age_group::Array{Int64,1},
                                    n_initial_asymp_by_age_group::Array{Int64,1},
                                    n_initial_symp_by_age_group::Array{Int64,1},
                                    n_initial_rec_by_age_group::Array{Int64,1},
                                    initialise_start_disease_state_flag::Bool,
                                    count::Int64,
                                    output::sim_outputs)

    # If assigning initial status by age, get vector of age group by person
    age_grp_by_person = zeros(Int64,n_persons)
    for person_itr = 1:n_persons
        age_grp_by_person[person_itr] = info_per_person[person_itr].age_group
    end

    # Get vector of IDs per age group
    person_IDs_by_age_group = Array{Array{Int64,1},1}(undef,3)
    for age_grp_itr = 1:3
        person_IDs_by_age_group[age_grp_itr] = findall(age_grp_by_person.==age_grp_itr)
    end

    # Iterate over each age group and assign initial states
    for age_grp_itr = 1:3
        choose_from_all_popn(rng,
                                person_IDs_by_age_group[age_grp_itr],
                                info_per_person,
                                states,
                                probasymp_adults,
                                probasymp_children,
                                infected_by,
                                n_initial_latent_by_age_group[age_grp_itr],
                                n_initial_asymp_by_age_group[age_grp_itr],
                                n_initial_symp_by_age_group[age_grp_itr],
                                n_initial_rec_by_age_group[age_grp_itr],
                                initialise_start_disease_state_flag,
                                count,
                                output)
    end

    return
end

#=
Fns to get number of individuals to be seeded in each non-susceptible state
=#

# Standard function setup
# Inputs
# rng::MersenneTwister - The random number generator
# n_persons::Int64 - Numer of individuals in the network
# n_ppl_by_age_grp::Array{Int64,1} - Number of people of age: [0-19yrs,20-64yrs,65+yrs]
# info_per_person::Array{person_info,1} - Structure will data specific to each person
# states::indiv_states - Record of info per individual
# probasymp_adults::Float64, probasymp_children::Float64 _ Asymp probability by age group
# infected_by::Array{Int64,1} - Record of who each individual was infected by
# recov_propn::Float64 - Proportion of population to begin in recovered state

# Outputs
# n_initial_latent, n_initial_asymp, n_initial_symp, n_initial_rec
#   -  Numbers to be seeded in the named disease state

# Set up with 5 asymps & 5 symp
function set_five_symp_five_asymp_initial(rng::MersenneTwister,
                                    n_persons::Int64,
                                    n_ppl_by_age_grp::Array{Int64,1},
                                    count::Int64,
                                    info_per_person::Array{person_info,1},
                                    states::indiv_states,
                                    probasymp_adults::Float64,
                                    probasymp_children::Float64,
                                    infected_by::Array{Int64,1},
                                    output::sim_outputs)

    # Set initial infected counts
    n_initial_latent = 0
    n_initial_asymp = 5
    n_initial_symp = 5

    # Set initial recovereds based on recov_propn
    recov_propn = 0.
    n_initial_rec = ceil(Int64,recov_propn*n_persons)

    # Select individuals from the population
    initialise_start_disease_state_flag = true
    relevant_person_IDs = collect(1:n_persons)
    choose_from_all_popn(rng,
                            relevant_person_IDs,
                            info_per_person,
                            states,
                            probasymp_adults,
                            probasymp_children,
                            infected_by,
                            n_initial_latent,
                            n_initial_asymp,
                            n_initial_symp,
                            n_initial_rec,
                            initialise_start_disease_state_flag,
                            count,
                            output)

    return n_initial_latent::Int64,
        n_initial_asymp::Int64,
        n_initial_symp::Int64,
        n_initial_rec::Int64
end

# Example to generate counts from distributions
function seed_states_with_uncertainty(rng::MersenneTwister,
                                    n_persons::Int64,
                                    n_ppl_by_age_grp::Array{Int64,1},
                                    count::Int64,
                                    info_per_person::Array{person_info,1},
                                    states::indiv_states,
                                    probasymp_adults::Float64,
                                    probasymp_children::Float64,
                                    infected_by::Array{Int64,1},
                                    output::sim_outputs)

    # Set distributions to draw counts from
    d_asymp = Uniform(0,10)
    d_symp = Uniform(0,3)

    # Set initial infected counts
    n_initial_latent = 0
    n_initial_asymp = round(Int64,rand(rng,d_asymp))
    n_initial_symp = round(Int64,rand(rng,d_symp))

    # Set initial recovereds based on recov_propn
    recov_propn = 0.
    n_initial_rec = ceil(Int64,recov_propn*n_persons)

    # Select individuals from the population
    initialise_start_disease_state_flag = true
    relevant_person_IDs = collect(1:n_persons)
    choose_from_all_popn(rng,
                            relevant_person_IDs,
                            info_per_person,
                            states,
                            probasymp_adults,
                            probasymp_children,
                            infected_by,
                            n_initial_latent,
                            n_initial_asymp,
                            n_initial_symp,
                            n_initial_rec,
                            initialise_start_disease_state_flag,
                            count,
                            output)

    return n_initial_latent,
        n_initial_asymp,
        n_initial_symp,
        n_initial_rec
end

# Generate counts using info on age group prevalence
function seed_states_using_age(rng::MersenneTwister,
                                    n_persons::Int64,
                                    n_ppl_by_age_grp::Array{Int64,1},
                                    count::Int64,
                                    info_per_person::Array{person_info,1},
                                    states::indiv_states,
                                    probasymp_adults::Float64,
                                    probasymp_children::Float64,
                                    infected_by::Array{Int64,1},
                                    output::sim_outputs)

    # Set proportion in each disease state by age group
    propn_initial_latent = [0.01,0.005,0.0025] # Latent: 1%, 0.5%, 0.25%
    propn_initial_asymp = [0.003,0.001,0.0005] # Overall infectious: 1%, 0.5%, 0.25%
    propn_initial_symp = [0.007,0.004,0.002]
    recov_propns = [0.25,0.25,0.15]

    # Get initial count in each disease state per age group
    n_initial_latent_by_age_group = ceil.(Int64,propn_initial_latent.*n_ppl_by_age_grp)
    n_initial_asymp_by_age_group = ceil.(Int64,propn_initial_asymp.*n_ppl_by_age_grp)
    n_initial_symp_by_age_group = ceil.(Int64,propn_initial_symp.*n_ppl_by_age_grp)
    n_initial_rec_by_age_group = ceil.(Int64,recov_propns.*n_ppl_by_age_grp)
    println("n_initial_latent_by_age_group: $n_initial_latent_by_age_group")
    println("n_initial_rec_by_age_group: $n_initial_rec_by_age_group")

    # Select individuals from the population based on age group
    initialise_start_disease_state_flag = false
    choose_from_age_groups(rng,
                            n_persons,
                            n_ppl_by_age_grp,
                            info_per_person,
                            states,
                            probasymp_adults,
                            probasymp_children,
                            infected_by,
                            n_initial_latent_by_age_group,
                            n_initial_asymp_by_age_group,
                            n_initial_symp_by_age_group,
                            n_initial_rec_by_age_group,
                            initialise_start_disease_state_flag,
                            count,
                            output)

    # # Get counts summed over age groups
    # n_initial_latent = sum(n_initial_latent_by_age_group)
    # n_initial_asymp = sum(n_initial_asymp_by_age_group)
    # n_initial_symp = sum(n_initial_symp_by_age_group)
    # n_initial_rec = sum(n_initial_rec_by_age_group)

    return n_initial_latent_by_age_group::Array{Int64,1},
            n_initial_asymp_by_age_group::Array{Int64,1},
            n_initial_symp_by_age_group::Array{Int64,1},
            n_initial_rec_by_age_group::Array{Int64,1}
end
