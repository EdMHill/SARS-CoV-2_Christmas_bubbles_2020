# SARS-CoV-2_christmas_bubbles_2020

This repository contains files for performing computational simulations of a stochastic individual based model of a population of households with three-age classes (0-19yrs, 20-64yrs, 65+yrs) to explore SARS-CoV-2 transmission for different household bubble strategies.

[![DOI](https://zenodo.org/badge/509207602.svg)](https://zenodo.org/badge/latestdoi/509207602)

The code was developed for the analysis presented in the scientific paper "Modelling the epidemiological implications for SARS-CoV-2 of Christmas household bubbles in England" by Edward M. Hill.

Model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project

Please find below an explainer of the directory structure within this repository.

## data
Directory containing household composition data from the 2011 census for England and Wales.

## results
Directory containing simulation outputs and the plot script (household_bubble_impact_plots.m). The sub-directory 'matlab_packages' has files that are loaded by the plot script file

## src

**run_stochastic_householdIBM.jl**  
The main model run file. Pass requested configuration to run_stochastic_household_IBM_fn (in include_files_stochastic_householdIBM/main_function.jl) and saves outputs into an MAT file.

**cluster_submit_job_files**  
Submission scripts for a computing cluster using the Slurm Workload Manager.

**include_files_stochastic_householdIBM**  
Houses function files to be used when running the uni model.


- **main_function.jl**   
    Outline of the code structure:.
    * Unpack required variables
    * Set the random number generator
    * Assign households (supporting functions in stochastic_householdIBM_support_fns.jl)
    * Generate household specific transmission (supporting functions in stochastic_householdIBM_support_fns.jl)
    * Initialise variables
    * Iterate over replicates
        - Reinitialisation phase
        - Set course of infection times (supporting functions in stochastic_householdIBM_support_fns.jl)
        - Update output time series with initial conditions
        - Reset contact tracing variables (supporting functions in stochastic_householdIBM_support_fns.jl)
        - Iterate over time
            - Reinitialise variables at start of timestep
            - Assign outputs
            - Increment counters (supporting functions in stochastic_householdIBM_support_fns.jl)
            - Increment infection process (if come to the end of latent time, move to infectious time etc.). Includes household infection loop. (supporting functions in stochastic_householdIBM_support_fns.jl)
            - Record whether individuals are in isolation.
            - Transmit infections (supporting functions in stochastic_householdIBM_support_fns.jl)
            - Perform contact tracing (supporting functions in contact_tracing_fns.jl)
            - Assign prevalence & isolation outputs
    * Return outputs

- **stochastic_householdIBM_support_fns.jl**   
    Supporting functions for running the stochastic IBM household bubble model.
    * Functions to generate households & person info:
        - sample_households (Sample the requested number of households)
        - generate_persons (Assign required number of people to each household)

    * Functions to increment state tracking variables
        - increment_counters! (Used each timestep to increase the value of time in state variables)

    * Transmission functions contained within this file include:
        - transmit_over!  (check infection to possible contacts)

    * Functions to reinitialise states at start of each run
        - reinitialise_node_states! (multiply time series vectors by 0)
        - reinitialise_CT_vars! (reset contact tracing associated variables)

    * Miscellaneous functions contained within this file include:
        - draw_sample_from_pmf!
        - set_infection_related_times! (set times to infection etc.: returns inftime, symptime, lattime, hh_isolation and delay_adherence)

- **seed_initial_states_fn.jl**   
    Function specifying how the initial disease states will be assigned in each simulation replicate.

- **contact_tracing_fns.jl**  
    Functions that are used for performing forward contact tracing from an identified infector.

- **stochastic_householdIBM_parametertypes.jl**  
    Defines the parameter types. Fields accessible with dot notation. Example using type household_info, with a variable named households. households.household_ID accesses the value in the household_ID field.    

- **populate_hh_visit_plan_fns.jl**  
    Functions to generate household bubbles and visitation schedules. Includes:
    * support_bubble_visit_plan! (Scenario A)
    * three_household_faithful_bubble_shorter_visit_plan! (Scenario B)
    * three_household_faithful_bubble_visit_plan! (Scenario C)
    * three_household_fixed_bubble_visit_plan! (Scenario D)
    * three_household_unfaithful_bubble_visit_plan! (Scenario E)
