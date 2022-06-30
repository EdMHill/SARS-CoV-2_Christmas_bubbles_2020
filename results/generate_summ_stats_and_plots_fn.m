%Purpose:
% Produce outputs for epidemiological outcomes under different 
% household bubbling assumptions that occur during the permitted five day period
%--------------------------------------------------------------------------

function [incidence_percentile_percentages_output,...
          age_split_incidence_percentile_percentages_output,...
          cumulative_infections_percentile_percentages_output,...
          age_split_cumulative_infections_percentile_percentages_output] =...
            generate_summ_stats_and_plots_fn(main_analysis_flag,...
                                             main_adherence_flag,...
                                             plot_percentages_flag,...
                                             generate_figures_flag,...
                                             save_file_flag)


    %% Load the relevant data files
    if main_adherence_flag == true % 70% adherence probability runs
        if main_analysis_flag == true
            support_bubble_only_data = load('MAT_files/support_bubble_visit_plan_main_analysis.mat');
            short_faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_shorter_visit_plan_main_analysis.mat');
            faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_visit_plan_main_analysis.mat');
            fixed_bubble_data = load('MAT_files/three_household_fixed_bubble_visit_plan_main_analysis.mat');
            unfaithful_bubble_data = load('MAT_files/three_household_unfaithful_bubble_visit_plan_main_analysis.mat');
        
            % Specify time index (row of arrays) to access
            row_idx = 6; % e.g. Entry 6 corresponds to day 5
        
        else
            support_bubble_only_data = load('MAT_files/support_bubble_visit_plan_alternative_analysis.mat');
            short_faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_shorter_visit_plan_alternative_analysis.mat');
            faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_visit_plan_alternative_analysis.mat');
            fixed_bubble_data = load('MAT_files/three_household_fixed_bubble_visit_plan_alternative_analysis.mat');
            unfaithful_bubble_data = load('MAT_files/three_household_unfaithful_bubble_visit_plan_alternative_analysis.mat');
        
            % Specify time index (row of arrays) to access
            row_idx = 16; % e.g. Entry 6 corresponds to day 5
        end
    else % 30% adherence probability runs
        if main_analysis_flag == true
            support_bubble_only_data = load('MAT_files/support_bubble_visit_plan_main_analysis_lower_adherence.mat');
            short_faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_shorter_visit_plan_main_analysis_lower_adherence.mat');
            faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_visit_plan_main_analysis_lower_adherence.mat');
            fixed_bubble_data = load('MAT_files/three_household_fixed_bubble_visit_plan_main_analysis_lower_adherence.mat');
            unfaithful_bubble_data = load('MAT_files/three_household_unfaithful_bubble_visit_plan_main_analysis_lower_adherence.mat');
        
            % Specify time index (row of arrays) to access
            row_idx = 6; % e.g. Entry 6 corresponds to day 5
        
        else
            support_bubble_only_data = load('MAT_files/support_bubble_visit_plan_alternative_analysis_lower_adherence.mat');
            short_faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_shorter_visit_plan_alternative_analysis_lower_adherence.mat');
            faithful_bubble_data = load('MAT_files/three_household_faithful_bubble_visit_plan_alternative_analysis_lower_adherence.mat');
            fixed_bubble_data = load('MAT_files/three_household_fixed_bubble_visit_plan_alternative_analysis_lower_adherence.mat');
            unfaithful_bubble_data = load('MAT_files/three_household_unfaithful_bubble_visit_plan_alternative_analysis_lower_adherence.mat');
        
            % Specify time index (row of arrays) to access
            row_idx = 16; % e.g. Entry 6 corresponds to day 5
        end
    end
    
    %% Set up variables used throughout script
    
    % Suffix for save_filename associated with calendar date
    save_filename_calender_date_suffixes = {'_23Dec2020';...
                                             '_24Dec2020';...
                                             '_25Dec2020';...
                                             '_26Dec2020';...
                                             '_27Dec2020'};
    save_filename_prevalence_suffixes = {'_23Dec2020_27Dec2020';...
                                             '_whole_time_horizon'};
    
    %% Compute summary statistics
    %% Age distribution of new latent infections on 23-27 December 2020
    
    % Get relevant values
    new_latent_23Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    
    new_latent_24Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    
    new_latent_25Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    
    new_latent_26Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    
    new_latent_27Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    
    % Normalise the values
    new_latent_23Dec2020_support_bubble_only_normalised = new_latent_23Dec2020_support_bubble_only./sum(new_latent_23Dec2020_support_bubble_only,2);
    new_latent_23Dec2020_short_faithful_bubble_normalised = new_latent_23Dec2020_short_faithful_bubble./sum(new_latent_23Dec2020_short_faithful_bubble,2);
    new_latent_23Dec2020_faithful_bubble_normalised = new_latent_23Dec2020_faithful_bubble./sum(new_latent_23Dec2020_faithful_bubble,2);
    new_latent_23Dec2020_fixed_bubble_normalised = new_latent_23Dec2020_fixed_bubble./sum(new_latent_23Dec2020_fixed_bubble,2);
    new_latent_23Dec2020_unfaithful_bubble_normalised = new_latent_23Dec2020_unfaithful_bubble./sum(new_latent_23Dec2020_unfaithful_bubble,2);
    
    new_latent_24Dec2020_support_bubble_only_normalised = new_latent_24Dec2020_support_bubble_only./sum(new_latent_24Dec2020_support_bubble_only,2);
    new_latent_24Dec2020_short_faithful_bubble_normalised = new_latent_24Dec2020_short_faithful_bubble./sum(new_latent_24Dec2020_short_faithful_bubble,2);
    new_latent_24Dec2020_faithful_bubble_normalised = new_latent_24Dec2020_faithful_bubble./sum(new_latent_24Dec2020_faithful_bubble,2);
    new_latent_24Dec2020_fixed_bubble_normalised = new_latent_24Dec2020_fixed_bubble./sum(new_latent_24Dec2020_fixed_bubble,2);
    new_latent_24Dec2020_unfaithful_bubble_normalised = new_latent_24Dec2020_unfaithful_bubble./sum(new_latent_24Dec2020_unfaithful_bubble,2);
    
    new_latent_25Dec2020_support_bubble_only_normalised = new_latent_25Dec2020_support_bubble_only./sum(new_latent_25Dec2020_support_bubble_only,2);
    new_latent_25Dec2020_short_faithful_bubble_normalised = new_latent_25Dec2020_short_faithful_bubble./sum(new_latent_25Dec2020_short_faithful_bubble,2);
    new_latent_25Dec2020_faithful_bubble_normalised = new_latent_25Dec2020_faithful_bubble./sum(new_latent_25Dec2020_faithful_bubble,2);
    new_latent_25Dec2020_fixed_bubble_normalised = new_latent_25Dec2020_fixed_bubble./sum(new_latent_25Dec2020_fixed_bubble,2);
    new_latent_25Dec2020_unfaithful_bubble_normalised = new_latent_25Dec2020_unfaithful_bubble./sum(new_latent_25Dec2020_unfaithful_bubble,2);
    
    new_latent_26Dec2020_support_bubble_only_normalised = new_latent_26Dec2020_support_bubble_only./sum(new_latent_26Dec2020_support_bubble_only,2);
    new_latent_26Dec2020_short_faithful_bubble_normalised = new_latent_26Dec2020_short_faithful_bubble./sum(new_latent_26Dec2020_short_faithful_bubble,2);
    new_latent_26Dec2020_faithful_bubble_normalised = new_latent_26Dec2020_faithful_bubble./sum(new_latent_26Dec2020_faithful_bubble,2);
    new_latent_26Dec2020_fixed_bubble_normalised = new_latent_26Dec2020_fixed_bubble./sum(new_latent_26Dec2020_fixed_bubble,2);
    new_latent_26Dec2020_unfaithful_bubble_normalised = new_latent_26Dec2020_unfaithful_bubble./sum(new_latent_26Dec2020_unfaithful_bubble,2);
    
    new_latent_27Dec2020_support_bubble_only_normalised = new_latent_27Dec2020_support_bubble_only./sum(new_latent_27Dec2020_support_bubble_only,2);
    new_latent_27Dec2020_short_faithful_bubble_normalised = new_latent_27Dec2020_short_faithful_bubble./sum(new_latent_27Dec2020_short_faithful_bubble,2);
    new_latent_27Dec2020_faithful_bubble_normalised = new_latent_27Dec2020_faithful_bubble./sum(new_latent_27Dec2020_faithful_bubble,2);
    new_latent_27Dec2020_fixed_bubble_normalised = new_latent_27Dec2020_fixed_bubble./sum(new_latent_27Dec2020_fixed_bubble,2);
    new_latent_27Dec2020_unfaithful_bubble_normalised = new_latent_27Dec2020_unfaithful_bubble./sum(new_latent_27Dec2020_unfaithful_bubble,2);
    
    
    %% Compute summary statistics
    %% New latent infections as a proportion of specified groups
    
    % Convert population count arrays to a double
    support_bubble_scen_n_ppl_by_age_grp = double(support_bubble_only_data.n_ppl_by_age_grp);
    short_faithful_bubble_scen_n_ppl_by_age_grp = double(short_faithful_bubble_data.n_ppl_by_age_grp);
    faithful_bubble_scen_n_ppl_by_age_grp = double(faithful_bubble_data.n_ppl_by_age_grp);
    fixed_bubble_scen_n_ppl_by_age_grp = double(fixed_bubble_data.n_ppl_by_age_grp);
    unfaithful_bubble_scen_n_ppl_by_age_grp = double(unfaithful_bubble_data.n_ppl_by_age_grp);
    
    % Get overall population per simulation
    support_bubble_scen_overall_popn = sum(support_bubble_scen_n_ppl_by_age_grp,2);
    short_faithful_bubble_scen_overall_popn = sum(short_faithful_bubble_scen_n_ppl_by_age_grp,2);
    faithful_bubble_scen_overall_popn = sum(faithful_bubble_scen_n_ppl_by_age_grp,2);
    fixed_bubble_scen_overall_popn = sum(fixed_bubble_scen_n_ppl_by_age_grp,2);
    unfaithful_bubble_scen_overall_popn = sum(unfaithful_bubble_scen_n_ppl_by_age_grp,2);
    
    % Get relevant latent values
    new_latent_23Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    new_latent_23Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-4,:,:)));
    
    all_new_latent_23Dec2020_support_bubble_only = double(support_bubble_only_data.newlat(row_idx-4,:)');
    all_new_latent_23Dec2020_short_faithful_bubble = double(short_faithful_bubble_data.newlat(row_idx-4,:)');
    all_new_latent_23Dec2020_faithful_bubble = double(faithful_bubble_data.newlat(row_idx-4,:)');
    all_new_latent_23Dec2020_fixed_bubble = double(faithful_bubble_data.newlat(row_idx-4,:)');
    all_new_latent_23Dec2020_unfaithful_bubble = double(unfaithful_bubble_data.newlat(row_idx-4,:)');
    
    new_latent_24Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    new_latent_24Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-3,:,:)));
    
    all_new_latent_24Dec2020_support_bubble_only = double(support_bubble_only_data.newlat(row_idx-3,:)');
    all_new_latent_24Dec2020_short_faithful_bubble = double(short_faithful_bubble_data.newlat(row_idx-3,:)');
    all_new_latent_24Dec2020_faithful_bubble = double(faithful_bubble_data.newlat(row_idx-3,:)');
    all_new_latent_24Dec2020_fixed_bubble = double(faithful_bubble_data.newlat(row_idx-3,:)');
    all_new_latent_24Dec2020_unfaithful_bubble = double(unfaithful_bubble_data.newlat(row_idx-3,:)');
    
    new_latent_25Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    new_latent_25Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-2,:,:)));
    
    all_new_latent_25Dec2020_support_bubble_only = double(support_bubble_only_data.newlat(row_idx-2,:)');
    all_new_latent_25Dec2020_short_faithful_bubble = double(short_faithful_bubble_data.newlat(row_idx-2,:)');
    all_new_latent_25Dec2020_faithful_bubble = double(faithful_bubble_data.newlat(row_idx-2,:)');
    all_new_latent_25Dec2020_fixed_bubble = double(faithful_bubble_data.newlat(row_idx-2,:)');
    all_new_latent_25Dec2020_unfaithful_bubble = double(unfaithful_bubble_data.newlat(row_idx-2,:)');
    
    new_latent_26Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    new_latent_26Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx-1,:,:)));
    
    all_new_latent_26Dec2020_support_bubble_only = double(support_bubble_only_data.newlat(row_idx-1,:)');
    all_new_latent_26Dec2020_short_faithful_bubble = double(short_faithful_bubble_data.newlat(row_idx-1,:)');
    all_new_latent_26Dec2020_faithful_bubble = double(faithful_bubble_data.newlat(row_idx-1,:)');
    all_new_latent_26Dec2020_fixed_bubble = double(faithful_bubble_data.newlat(row_idx-1,:)');
    all_new_latent_26Dec2020_unfaithful_bubble = double(unfaithful_bubble_data.newlat(row_idx-1,:)');
    
    new_latent_27Dec2020_support_bubble_only = double(squeeze(support_bubble_only_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_short_faithful_bubble = double(squeeze(short_faithful_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_faithful_bubble = double(squeeze(faithful_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_fixed_bubble = double(squeeze(fixed_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    new_latent_27Dec2020_unfaithful_bubble = double(squeeze(unfaithful_bubble_data.newlat_by_age_grp(row_idx,:,:)));
    
    all_new_latent_27Dec2020_support_bubble_only = double(support_bubble_only_data.newlat(row_idx,:)');
    all_new_latent_27Dec2020_short_faithful_bubble = double(short_faithful_bubble_data.newlat(row_idx,:)');
    all_new_latent_27Dec2020_faithful_bubble = double(faithful_bubble_data.newlat(row_idx,:)');
    all_new_latent_27Dec2020_fixed_bubble = double(faithful_bubble_data.newlat(row_idx,:)');
    all_new_latent_27Dec2020_unfaithful_bubble = double(unfaithful_bubble_data.newlat(row_idx,:)');
    
    % Group the input data
    new_latent_23Dec2020_propn_agegrp_support_bubble = new_latent_23Dec2020_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    new_latent_23Dec2020_propn_agegrp_short_faithful_bubble = new_latent_23Dec2020_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_23Dec2020_propn_agegrp_faithful_bubble = new_latent_23Dec2020_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_23Dec2020_propn_agegrp_fixed_bubble = new_latent_23Dec2020_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    new_latent_23Dec2020_propn_agegrp_unfaithful_bubble = new_latent_23Dec2020_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    new_latent_24Dec2020_propn_agegrp_support_bubble = new_latent_24Dec2020_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    new_latent_24Dec2020_propn_agegrp_short_faithful_bubble = new_latent_24Dec2020_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_24Dec2020_propn_agegrp_faithful_bubble = new_latent_24Dec2020_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_24Dec2020_propn_agegrp_fixed_bubble = new_latent_24Dec2020_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    new_latent_24Dec2020_propn_agegrp_unfaithful_bubble = new_latent_24Dec2020_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    new_latent_25Dec2020_propn_agegrp_support_bubble = new_latent_25Dec2020_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    new_latent_25Dec2020_propn_agegrp_short_faithful_bubble = new_latent_25Dec2020_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_25Dec2020_propn_agegrp_faithful_bubble = new_latent_25Dec2020_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_25Dec2020_propn_agegrp_fixed_bubble = new_latent_25Dec2020_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    new_latent_25Dec2020_propn_agegrp_unfaithful_bubble = new_latent_25Dec2020_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    new_latent_26Dec2020_propn_agegrp_support_bubble = new_latent_26Dec2020_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    new_latent_26Dec2020_propn_agegrp_short_faithful_bubble = new_latent_26Dec2020_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_26Dec2020_propn_agegrp_faithful_bubble = new_latent_26Dec2020_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_26Dec2020_propn_agegrp_fixed_bubble = new_latent_26Dec2020_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    new_latent_26Dec2020_propn_agegrp_unfaithful_bubble = new_latent_26Dec2020_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    new_latent_27Dec2020_propn_agegrp_support_bubble = new_latent_27Dec2020_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    new_latent_27Dec2020_propn_agegrp_short_faithful_bubble = new_latent_27Dec2020_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_27Dec2020_propn_agegrp_faithful_bubble = new_latent_27Dec2020_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    new_latent_27Dec2020_propn_agegrp_fixed_bubble = new_latent_27Dec2020_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    new_latent_27Dec2020_propn_agegrp_unfaithful_bubble = new_latent_27Dec2020_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    % Overall
    new_latent_23Dec2020_propn_overall_popn_support_bubble = all_new_latent_23Dec2020_support_bubble_only./support_bubble_scen_overall_popn;
    new_latent_23Dec2020_propn_overall_popn_short_faithful_bubble = all_new_latent_23Dec2020_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    new_latent_23Dec2020_propn_overall_popn_faithful_bubble = all_new_latent_23Dec2020_faithful_bubble./faithful_bubble_scen_overall_popn;
    new_latent_23Dec2020_propn_overall_popn_fixed_bubble = all_new_latent_23Dec2020_fixed_bubble./fixed_bubble_scen_overall_popn;
    new_latent_23Dec2020_propn_overall_popn_unfaithful_bubble = all_new_latent_23Dec2020_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
    
    new_latent_24Dec2020_propn_overall_popn_support_bubble = all_new_latent_24Dec2020_support_bubble_only./support_bubble_scen_overall_popn;
    new_latent_24Dec2020_propn_overall_popn_short_faithful_bubble = all_new_latent_24Dec2020_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    new_latent_24Dec2020_propn_overall_popn_faithful_bubble = all_new_latent_24Dec2020_faithful_bubble./faithful_bubble_scen_overall_popn;
    new_latent_24Dec2020_propn_overall_popn_fixed_bubble = all_new_latent_24Dec2020_fixed_bubble./fixed_bubble_scen_overall_popn;
    new_latent_24Dec2020_propn_overall_popn_unfaithful_bubble = all_new_latent_24Dec2020_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
    
    new_latent_25Dec2020_propn_overall_popn_support_bubble = all_new_latent_25Dec2020_support_bubble_only./support_bubble_scen_overall_popn;
    new_latent_25Dec2020_propn_overall_popn_short_faithful_bubble = all_new_latent_25Dec2020_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    new_latent_25Dec2020_propn_overall_popn_faithful_bubble = all_new_latent_25Dec2020_faithful_bubble./faithful_bubble_scen_overall_popn;
    new_latent_25Dec2020_propn_overall_popn_fixed_bubble = all_new_latent_25Dec2020_fixed_bubble./fixed_bubble_scen_overall_popn;
    new_latent_25Dec2020_propn_overall_popn_unfaithful_bubble = all_new_latent_25Dec2020_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
    
    new_latent_26Dec2020_propn_overall_popn_support_bubble = all_new_latent_26Dec2020_support_bubble_only./support_bubble_scen_overall_popn;
    new_latent_26Dec2020_propn_overall_popn_short_faithful_bubble = all_new_latent_26Dec2020_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    new_latent_26Dec2020_propn_overall_popn_faithful_bubble = all_new_latent_26Dec2020_faithful_bubble./faithful_bubble_scen_overall_popn;
    new_latent_26Dec2020_propn_overall_popn_fixed_bubble = all_new_latent_26Dec2020_fixed_bubble./fixed_bubble_scen_overall_popn;
    new_latent_26Dec2020_propn_overall_popn_unfaithful_bubble = all_new_latent_26Dec2020_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
    
    new_latent_27Dec2020_propn_overall_popn_support_bubble = all_new_latent_27Dec2020_support_bubble_only./support_bubble_scen_overall_popn;
    new_latent_27Dec2020_propn_overall_popn_short_faithful_bubble = all_new_latent_27Dec2020_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    new_latent_27Dec2020_propn_overall_popn_faithful_bubble = all_new_latent_27Dec2020_faithful_bubble./faithful_bubble_scen_overall_popn;
    new_latent_27Dec2020_propn_overall_popn_fixed_bubble = all_new_latent_27Dec2020_fixed_bubble./fixed_bubble_scen_overall_popn;
    new_latent_27Dec2020_propn_overall_popn_unfaithful_bubble = all_new_latent_27Dec2020_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
    
    %% Compute summary statistics
    %% Total new infections over the time horizon as a proportion of specified groups
    
    % Convert population count arrays to a double
    support_bubble_scen_n_ppl_by_age_grp = double(support_bubble_only_data.n_ppl_by_age_grp);
    short_faithful_bubble_scen_n_ppl_by_age_grp = double(short_faithful_bubble_data.n_ppl_by_age_grp);
    faithful_bubble_scen_n_ppl_by_age_grp = double(faithful_bubble_data.n_ppl_by_age_grp);
    fixed_bubble_scen_n_ppl_by_age_grp = double(fixed_bubble_data.n_ppl_by_age_grp);
    unfaithful_bubble_scen_n_ppl_by_age_grp = double(unfaithful_bubble_data.n_ppl_by_age_grp);
    
    % Get overall population per simulation
    support_bubble_scen_overall_popn = sum(support_bubble_scen_n_ppl_by_age_grp,2);
    short_faithful_bubble_scen_overall_popn = sum(short_faithful_bubble_scen_n_ppl_by_age_grp,2);
    faithful_bubble_scen_overall_popn = sum(faithful_bubble_scen_n_ppl_by_age_grp,2);
    fixed_bubble_scen_overall_popn = sum(fixed_bubble_scen_n_ppl_by_age_grp,2);
    unfaithful_bubble_scen_overall_popn = sum(unfaithful_bubble_scen_n_ppl_by_age_grp,2);
    
    % Get relevant new infected values
    newinf_support_bubble_only = squeeze(sum(support_bubble_only_data.newlat_by_age_grp));
    newinf_short_faithful_bubble = squeeze(sum(short_faithful_bubble_data.newlat_by_age_grp));
    newinf_faithful_bubble = squeeze(sum(faithful_bubble_data.newlat_by_age_grp));
    newinf_fixed_bubble = squeeze(sum(fixed_bubble_data.newlat_by_age_grp));
    newinf_unfaithful_bubble = squeeze(sum(unfaithful_bubble_data.newlat_by_age_grp));
    
    newinf_christmas_bubble_window_support_bubble_only = squeeze(sum(support_bubble_only_data.newlat_by_age_grp(row_idx-4:row_idx,:,:)));
    newinf_christmas_bubble_window_short_faithful_bubble = squeeze(sum(short_faithful_bubble_data.newlat_by_age_grp(row_idx-4:row_idx,:,:)));
    newinf_christmas_bubble_window_faithful_bubble = squeeze(sum(faithful_bubble_data.newlat_by_age_grp(row_idx-4:row_idx,:,:)));
    newinf_christmas_bubble_window_fixed_bubble = squeeze(sum(fixed_bubble_data.newlat_by_age_grp(row_idx-4:row_idx,:,:)));
    newinf_christmas_bubble_window_unfaithful_bubble = squeeze(sum(unfaithful_bubble_data.newlat_by_age_grp(row_idx-4:row_idx,:,:)));
    
    all_newinf_support_bubble_only = sum(support_bubble_only_data.newlat)';
    all_newinf_short_faithful_bubble = sum(short_faithful_bubble_data.newlat)';
    all_newinf_faithful_bubble = sum(faithful_bubble_data.newlat)';
    all_newinf_fixed_bubble = sum(fixed_bubble_data.newlat)';
    all_newinf_unfaithful_bubble = sum(unfaithful_bubble_data.newlat)';
    
    all_newinf_christmas_bubble_window_support_bubble_only = sum(support_bubble_only_data.newlat(row_idx-4:row_idx,:))';
    all_newinf_christmas_bubble_window_short_faithful_bubble = sum(short_faithful_bubble_data.newlat(row_idx-4:row_idx,:))';
    all_newinf_christmas_bubble_window_faithful_bubble = sum(faithful_bubble_data.newlat(row_idx-4:row_idx,:))';
    all_newinf_christmas_bubble_window_fixed_bubble = sum(fixed_bubble_data.newlat(row_idx-4:row_idx,:))';
    all_newinf_christmas_bubble_window_unfaithful_bubble = sum(unfaithful_bubble_data.newlat(row_idx-4:row_idx,:))';
    
    % By age
    newinf_propn_agegrp_support_bubble = newinf_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    newinf_propn_agegrp_short_faithful_bubble = newinf_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    newinf_propn_agegrp_faithful_bubble = newinf_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    newinf_propn_agegrp_fixed_bubble = newinf_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    newinf_propn_agegrp_unfaithful_bubble = newinf_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    newinf_xmas_bubble_window_propn_agegrp_support_bubble = newinf_christmas_bubble_window_support_bubble_only./support_bubble_scen_n_ppl_by_age_grp;
    newinf_xmas_bubble_window_propn_agegrp_short_bubble = newinf_christmas_bubble_window_short_faithful_bubble./short_faithful_bubble_scen_n_ppl_by_age_grp;
    newinf_xmas_bubble_window_propn_agegrp_faithful_bubble = newinf_christmas_bubble_window_faithful_bubble./faithful_bubble_scen_n_ppl_by_age_grp;
    newinf_xmas_bubble_window_propn_agegrp_fixed_bubble = newinf_christmas_bubble_window_fixed_bubble./fixed_bubble_scen_n_ppl_by_age_grp;
    newinf_xmas_bubble_window_propn_agegrp_unfaithful_bubble = newinf_christmas_bubble_window_unfaithful_bubble./unfaithful_bubble_scen_n_ppl_by_age_grp;
    
    % Overall
    newinf_propn_overall_popn_support_bubble = all_newinf_support_bubble_only./support_bubble_scen_overall_popn;
    newinf_propn_overall_popn_short_faithful_bubble = all_newinf_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    newinf_propn_overall_popn_faithful_bubble = all_newinf_faithful_bubble./faithful_bubble_scen_overall_popn;
    newinf_propn_overall_popn_fixed_bubble = all_newinf_fixed_bubble./fixed_bubble_scen_overall_popn;
    newinf_propn_overall_popn_unfaithful_bubble = all_newinf_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
            
    newinf_xmas_bubble_window_propn_overall_popn_support_bubble = all_newinf_christmas_bubble_window_support_bubble_only./support_bubble_scen_overall_popn;
    newinf_xmas_bubble_window_propn_overall_popn_short_bubble = all_newinf_christmas_bubble_window_short_faithful_bubble./short_faithful_bubble_scen_overall_popn;
    newinf_xmas_bubble_window_propn_overall_popn_faithful_bubble = all_newinf_christmas_bubble_window_faithful_bubble./faithful_bubble_scen_overall_popn;
    newinf_xmas_bubble_window_propn_overall_popn_fixed_bubble = all_newinf_christmas_bubble_window_fixed_bubble./fixed_bubble_scen_overall_popn;
    newinf_xmas_bubble_window_propn_overall_popn_unfaithful_bubble = all_newinf_christmas_bubble_window_unfaithful_bubble./unfaithful_bubble_scen_overall_popn;
    
    
    %% Calculate summary statistics for tables
    
    % Set percentile values to compute
    percentile_vals = [2.5 50 97.5];
    
    % Incidence on 27Dec2020
    scen_A_incidence_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_support_bubble new_latent_27Dec2020_propn_overall_popn_support_bubble];
    scen_B_incidence_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_short_faithful_bubble new_latent_27Dec2020_propn_overall_popn_short_faithful_bubble];
    scen_C_incidence_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_faithful_bubble new_latent_27Dec2020_propn_overall_popn_faithful_bubble];
    scen_D_incidence_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_fixed_bubble new_latent_27Dec2020_propn_overall_popn_fixed_bubble];
    scen_E_incidence_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_unfaithful_bubble new_latent_27Dec2020_propn_overall_popn_unfaithful_bubble];
    
    scen_A_incidence_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_support_bubble new_latent_26Dec2020_propn_overall_popn_support_bubble];
    scen_B_incidence_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_short_faithful_bubble new_latent_26Dec2020_propn_overall_popn_short_faithful_bubble];
    scen_C_incidence_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_faithful_bubble new_latent_26Dec2020_propn_overall_popn_faithful_bubble];
    scen_D_incidence_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_fixed_bubble new_latent_26Dec2020_propn_overall_popn_fixed_bubble];
    scen_E_incidence_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_unfaithful_bubble new_latent_26Dec2020_propn_overall_popn_unfaithful_bubble];
    
    scen_A_incidence_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_support_bubble new_latent_25Dec2020_propn_overall_popn_support_bubble];
    scen_B_incidence_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_short_faithful_bubble new_latent_25Dec2020_propn_overall_popn_short_faithful_bubble];
    scen_C_incidence_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_faithful_bubble new_latent_25Dec2020_propn_overall_popn_faithful_bubble];
    scen_D_incidence_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_fixed_bubble new_latent_25Dec2020_propn_overall_popn_fixed_bubble];
    scen_E_incidence_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_unfaithful_bubble new_latent_25Dec2020_propn_overall_popn_unfaithful_bubble];
    
    scen_A_incidence_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_support_bubble new_latent_24Dec2020_propn_overall_popn_support_bubble];
    scen_B_incidence_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_short_faithful_bubble new_latent_24Dec2020_propn_overall_popn_short_faithful_bubble];
    scen_C_incidence_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_faithful_bubble new_latent_24Dec2020_propn_overall_popn_faithful_bubble];
    scen_D_incidence_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_fixed_bubble new_latent_24Dec2020_propn_overall_popn_fixed_bubble];
    scen_E_incidence_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_unfaithful_bubble new_latent_24Dec2020_propn_overall_popn_unfaithful_bubble];
    
    scen_A_incidence_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_support_bubble new_latent_23Dec2020_propn_overall_popn_support_bubble];
    scen_B_incidence_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_short_faithful_bubble new_latent_23Dec2020_propn_overall_popn_short_faithful_bubble];
    scen_C_incidence_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_faithful_bubble new_latent_23Dec2020_propn_overall_popn_faithful_bubble];
    scen_D_incidence_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_fixed_bubble new_latent_23Dec2020_propn_overall_popn_fixed_bubble];
    scen_E_incidence_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_unfaithful_bubble new_latent_23Dec2020_propn_overall_popn_unfaithful_bubble];
    
    % Age group split of new cases on 27Dec2020
    scen_A_age_split_incidence_27Dec2020 = new_latent_27Dec2020_support_bubble_only_normalised;
    scen_B_age_split_incidence_27Dec2020 = new_latent_27Dec2020_short_faithful_bubble_normalised;
    scen_C_age_split_incidence_27Dec2020 = new_latent_27Dec2020_faithful_bubble_normalised;
    scen_D_age_split_incidence_27Dec2020 = new_latent_27Dec2020_fixed_bubble_normalised;
    scen_E_age_split_incidence_27Dec2020 = new_latent_27Dec2020_unfaithful_bubble_normalised;
     
    scen_A_age_split_incidence_26Dec2020 = new_latent_26Dec2020_support_bubble_only_normalised;
    scen_B_age_split_incidence_26Dec2020 = new_latent_26Dec2020_short_faithful_bubble_normalised;
    scen_C_age_split_incidence_26Dec2020 = new_latent_26Dec2020_faithful_bubble_normalised;
    scen_D_age_split_incidence_26Dec2020 = new_latent_26Dec2020_fixed_bubble_normalised;
    scen_E_age_split_incidence_26Dec2020 = new_latent_26Dec2020_unfaithful_bubble_normalised;
        
    scen_A_age_split_incidence_25Dec2020 = new_latent_25Dec2020_support_bubble_only_normalised;
    scen_B_age_split_incidence_25Dec2020 = new_latent_25Dec2020_short_faithful_bubble_normalised;
    scen_C_age_split_incidence_25Dec2020 = new_latent_25Dec2020_faithful_bubble_normalised;
    scen_D_age_split_incidence_25Dec2020 = new_latent_25Dec2020_fixed_bubble_normalised;
    scen_E_age_split_incidence_25Dec2020 = new_latent_25Dec2020_unfaithful_bubble_normalised;
                 
    scen_A_age_split_incidence_24Dec2020 = new_latent_24Dec2020_support_bubble_only_normalised;
    scen_B_age_split_incidence_24Dec2020 = new_latent_24Dec2020_short_faithful_bubble_normalised;
    scen_C_age_split_incidence_24Dec2020 = new_latent_24Dec2020_faithful_bubble_normalised;
    scen_D_age_split_incidence_24Dec2020 = new_latent_24Dec2020_fixed_bubble_normalised;
    scen_E_age_split_incidence_24Dec2020 = new_latent_24Dec2020_unfaithful_bubble_normalised;
         
    scen_A_age_split_incidence_23Dec2020 = new_latent_23Dec2020_support_bubble_only_normalised;
    scen_B_age_split_incidence_23Dec2020 = new_latent_23Dec2020_short_faithful_bubble_normalised;
    scen_C_age_split_incidence_23Dec2020 = new_latent_23Dec2020_faithful_bubble_normalised;
    scen_D_age_split_incidence_23Dec2020 = new_latent_23Dec2020_fixed_bubble_normalised;
    scen_E_age_split_incidence_23Dec2020 = new_latent_23Dec2020_unfaithful_bubble_normalised;
                      
    % Cumulative infections
    scen_A_cumulative_infections = [newinf_propn_agegrp_support_bubble newinf_propn_overall_popn_support_bubble];
    scen_B_cumulative_infections = [newinf_propn_agegrp_short_faithful_bubble newinf_propn_overall_popn_short_faithful_bubble];
    scen_C_cumulative_infections = [newinf_propn_agegrp_faithful_bubble newinf_propn_overall_popn_faithful_bubble];
    scen_D_cumulative_infections = [newinf_propn_agegrp_fixed_bubble newinf_propn_overall_popn_fixed_bubble];
    scen_E_cumulative_infections = [newinf_propn_agegrp_unfaithful_bubble newinf_propn_overall_popn_unfaithful_bubble];
    
    scen_A_cumulative_infections_xmas_window = [newinf_xmas_bubble_window_propn_agegrp_support_bubble newinf_xmas_bubble_window_propn_overall_popn_support_bubble];
    scen_B_cumulative_infections_xmas_window = [newinf_xmas_bubble_window_propn_agegrp_short_bubble newinf_xmas_bubble_window_propn_overall_popn_short_bubble];
    scen_C_cumulative_infections_xmas_window = [newinf_xmas_bubble_window_propn_agegrp_faithful_bubble newinf_xmas_bubble_window_propn_overall_popn_faithful_bubble];
    scen_D_cumulative_infections_xmas_window = [newinf_xmas_bubble_window_propn_agegrp_fixed_bubble newinf_xmas_bubble_window_propn_overall_popn_fixed_bubble];
    scen_E_cumulative_infections_xmas_window = [newinf_xmas_bubble_window_propn_agegrp_unfaithful_bubble newinf_xmas_bubble_window_propn_overall_popn_unfaithful_bubble];
    
    % Age split for cumulative infections
    scen_A_age_split_cumulative_infections = newinf_propn_agegrp_support_bubble./sum(newinf_propn_agegrp_support_bubble,2);
    scen_B_age_split_cumulative_infections = newinf_propn_agegrp_short_faithful_bubble./sum(newinf_propn_agegrp_short_faithful_bubble,2);
    scen_C_age_split_cumulative_infections = newinf_propn_agegrp_faithful_bubble./sum(newinf_propn_agegrp_faithful_bubble,2);
    scen_D_age_split_cumulative_infections = newinf_propn_agegrp_fixed_bubble./sum(newinf_propn_agegrp_fixed_bubble,2);
    scen_E_age_split_cumulative_infections = newinf_propn_agegrp_unfaithful_bubble./sum(newinf_propn_agegrp_unfaithful_bubble,2);
    
    scen_A_age_split_cumulative_infections_xmas_window = newinf_xmas_bubble_window_propn_agegrp_support_bubble./sum(newinf_xmas_bubble_window_propn_agegrp_support_bubble,2);
    scen_B_age_split_cumulative_infections_xmas_window = newinf_xmas_bubble_window_propn_agegrp_short_bubble./sum(newinf_xmas_bubble_window_propn_agegrp_short_bubble,2);
    scen_C_age_split_cumulative_infections_xmas_window = newinf_xmas_bubble_window_propn_agegrp_faithful_bubble./sum(newinf_xmas_bubble_window_propn_agegrp_faithful_bubble,2);
    scen_D_age_split_cumulative_infections_xmas_window = newinf_xmas_bubble_window_propn_agegrp_fixed_bubble./sum(newinf_xmas_bubble_window_propn_agegrp_fixed_bubble,2);
    scen_E_age_split_cumulative_infections_xmas_window = newinf_xmas_bubble_window_propn_agegrp_unfaithful_bubble./sum(newinf_xmas_bubble_window_propn_agegrp_unfaithful_bubble,2);
    
    % Stack 2D arrays into 3D arrays: Incidence
    all_scen_incidence_27Dec2020 = cat(3,scen_A_incidence_27Dec2020,...
                                            scen_B_incidence_27Dec2020,...
                                            scen_C_incidence_27Dec2020,...
                                            scen_D_incidence_27Dec2020,...
                                            scen_E_incidence_27Dec2020); 
    
    all_scen_incidence_26Dec2020 = cat(3,scen_A_incidence_26Dec2020,...
                                            scen_B_incidence_26Dec2020,...
                                            scen_C_incidence_26Dec2020,...
                                            scen_D_incidence_26Dec2020,...
                                            scen_E_incidence_26Dec2020); 
    
    all_scen_incidence_25Dec2020 = cat(3,scen_A_incidence_25Dec2020,...
                                            scen_B_incidence_25Dec2020,...
                                            scen_C_incidence_25Dec2020,...
                                            scen_D_incidence_25Dec2020,...
                                            scen_E_incidence_25Dec2020); 
    
    all_scen_incidence_24Dec2020 = cat(3,scen_A_incidence_24Dec2020,...
                                            scen_B_incidence_24Dec2020,...
                                            scen_C_incidence_24Dec2020,...
                                            scen_D_incidence_24Dec2020,...
                                            scen_E_incidence_24Dec2020); 
    
    all_scen_incidence_23Dec2020 = cat(3,scen_A_incidence_23Dec2020,...
                                            scen_B_incidence_23Dec2020,...
                                            scen_C_incidence_23Dec2020,...
                                            scen_D_incidence_23Dec2020,...
                                            scen_E_incidence_23Dec2020); 
    
    all_scen_age_split_incidence_27Dec2020 = cat(3,scen_A_age_split_incidence_27Dec2020,...
                                                    scen_B_age_split_incidence_27Dec2020,...
                                                    scen_C_age_split_incidence_27Dec2020,...
                                                    scen_D_age_split_incidence_27Dec2020,...
                                                    scen_E_age_split_incidence_27Dec2020); 
    
    all_scen_age_split_incidence_26Dec2020 = cat(3,scen_A_age_split_incidence_26Dec2020,...
                                                    scen_B_age_split_incidence_26Dec2020,...
                                                    scen_C_age_split_incidence_26Dec2020,...
                                                    scen_D_age_split_incidence_26Dec2020,...
                                                    scen_E_age_split_incidence_26Dec2020);
    
    all_scen_age_split_incidence_25Dec2020 = cat(3,scen_A_age_split_incidence_25Dec2020,...
                                                    scen_B_age_split_incidence_25Dec2020,...
                                                    scen_C_age_split_incidence_25Dec2020,...
                                                    scen_D_age_split_incidence_25Dec2020,...
                                                    scen_E_age_split_incidence_25Dec2020);
    
    
    all_scen_age_split_incidence_24Dec2020 = cat(3,scen_A_age_split_incidence_24Dec2020,...
                                                    scen_B_age_split_incidence_24Dec2020,...
                                                    scen_C_age_split_incidence_24Dec2020,...
                                                    scen_D_age_split_incidence_24Dec2020,...
                                                    scen_E_age_split_incidence_24Dec2020);
    
    all_scen_age_split_incidence_23Dec2020 = cat(3,scen_A_age_split_incidence_23Dec2020,...
                                                    scen_B_age_split_incidence_23Dec2020,...
                                                    scen_C_age_split_incidence_23Dec2020,...
                                                    scen_D_age_split_incidence_23Dec2020,...
                                                    scen_E_age_split_incidence_23Dec2020);
    
    % Stack 2D arrays into 3D arrays: Prevalence
    all_scen_cumulative_infections = cat(3,scen_A_cumulative_infections,...
                                            scen_B_cumulative_infections,...
                                            scen_C_cumulative_infections,...
                                            scen_D_cumulative_infections,...
                                            scen_E_cumulative_infections);  
    
    all_scen_cumulative_infections_xmas_window = cat(3,scen_A_cumulative_infections_xmas_window,...
                                                        scen_B_cumulative_infections_xmas_window,...
                                                        scen_C_cumulative_infections_xmas_window,...
                                                        scen_D_cumulative_infections_xmas_window,...
                                                        scen_E_cumulative_infections_xmas_window);
    
    all_scen_age_split_cumulative_infections = cat(3,scen_A_age_split_cumulative_infections,...
                                            scen_B_age_split_cumulative_infections,...
                                            scen_C_age_split_cumulative_infections,...
                                            scen_D_age_split_cumulative_infections,...
                                            scen_E_age_split_cumulative_infections);  
    
    all_scen_age_split_cumulative_infections_xmas_window = cat(3,scen_A_age_split_cumulative_infections_xmas_window,...
                                                        scen_B_age_split_cumulative_infections_xmas_window,...
                                                        scen_C_age_split_cumulative_infections_xmas_window,...
                                                        scen_D_age_split_cumulative_infections_xmas_window,...
                                                        scen_E_age_split_cumulative_infections_xmas_window);
    
    % Get percentiles outputs for proportion summary statistics
    incidence_27Dec2020_percentiles = prctile(all_scen_incidence_27Dec2020,percentile_vals);
    ages_split_27Dec2020_incidence_percentiles = prctile(all_scen_age_split_incidence_27Dec2020,percentile_vals);
    
    incidence_26Dec2020_percentiles = prctile(all_scen_incidence_26Dec2020,percentile_vals);
    ages_split_26Dec2020_incidence_percentiles = prctile(all_scen_age_split_incidence_26Dec2020,percentile_vals);
    
    incidence_25Dec2020_percentiles = prctile(all_scen_incidence_25Dec2020,percentile_vals);
    ages_split_25Dec2020_incidence_percentiles = prctile(all_scen_age_split_incidence_25Dec2020,percentile_vals);
    
    incidence_24Dec2020_percentiles = prctile(all_scen_incidence_24Dec2020,percentile_vals);
    ages_split_24Dec2020_incidence_percentiles = prctile(all_scen_age_split_incidence_24Dec2020,percentile_vals);
    
    incidence_23Dec2020_percentiles = prctile(all_scen_incidence_23Dec2020,percentile_vals);
    ages_split_23Dec2020_incidence_percentiles = prctile(all_scen_age_split_incidence_23Dec2020,percentile_vals);
    
    cumulative_infections_percentiles = prctile(all_scen_cumulative_infections,percentile_vals);
    ages_split_cumulative_infections_percentiles = prctile(all_scen_age_split_cumulative_infections,percentile_vals);
    
    cumulative_infections_xmas_window_percentiles = prctile(all_scen_cumulative_infections_xmas_window,percentile_vals);
    ages_split_cumulative_infections_xmas_window_percentiles = prctile(all_scen_age_split_cumulative_infections_xmas_window,percentile_vals);
    
    % To get percentage output, multiply percentile values by 100
    incidence_27Dec2020_percentile_percentages = incidence_27Dec2020_percentiles*100;
    ages_split_27Dec2020_incidence_percentile_percentages = ages_split_27Dec2020_incidence_percentiles*100;
    
    incidence_26Dec2020_percentile_percentages = incidence_26Dec2020_percentiles*100;
    ages_split_26Dec2020_incidence_percentile_percentages = ages_split_26Dec2020_incidence_percentiles*100;
    
    incidence_25Dec2020_percentile_percentages = incidence_25Dec2020_percentiles*100;
    ages_split_25Dec2020_incidence_percentile_percentages = ages_split_25Dec2020_incidence_percentiles*100;
    
    incidence_24Dec2020_percentile_percentages = incidence_24Dec2020_percentiles*100;
    ages_split_24Dec2020_incidence_percentile_percentages = ages_split_24Dec2020_incidence_percentiles*100;
    
    incidence_23Dec2020_percentile_percentages = incidence_23Dec2020_percentiles*100;
    ages_split_23Dec2020_incidence_percentile_percentages = ages_split_23Dec2020_incidence_percentiles*100;
    
    cumulative_infections_percentile_percentages = cumulative_infections_percentiles*100;
    ages_split_cumulative_infections_percentile_percentages = ages_split_cumulative_infections_percentiles*100;
    
    cumulative_infections_xmas_window_percentile_percentages = cumulative_infections_xmas_window_percentiles*100;
    ages_split_cumul_infections_xmas_window_percentile_percentages = ages_split_cumulative_infections_xmas_window_percentiles*100;
    
    %% Produce summary table output
    
    % Specify number of scenarios per age group
    n_scens = 5;
    
    % Initialise cells to write percentage results into
    incidence_23Dec2020_percentile_percentages_cell = cell(size(incidence_23Dec2020_percentile_percentages,2),5);
    incidence_24Dec2020_percentile_percentages_cell = cell(size(incidence_24Dec2020_percentile_percentages,2),5);
    incidence_25Dec2020_percentile_percentages_cell = cell(size(incidence_25Dec2020_percentile_percentages,2),5);
    incidence_26Dec2020_percentile_percentages_cell = cell(size(incidence_26Dec2020_percentile_percentages,2),5);
    incidence_27Dec2020_percentile_percentages_cell = cell(size(incidence_27Dec2020_percentile_percentages,2),5);
    
    ages_split_23Dec2020_incidence_percentages_cell = cell(size(ages_split_23Dec2020_incidence_percentile_percentages,2),5);
    ages_split_24Dec2020_incidence_percentages_cell = cell(size(ages_split_24Dec2020_incidence_percentile_percentages,2),5);
    ages_split_25Dec2020_incidence_percentages_cell = cell(size(ages_split_25Dec2020_incidence_percentile_percentages,2),5);
    ages_split_26Dec2020_incidence_percentages_cell = cell(size(ages_split_26Dec2020_incidence_percentile_percentages,2),5);
    ages_split_27Dec2020_incidence_percentages_cell = cell(size(ages_split_27Dec2020_incidence_percentile_percentages,2),5);
    
    cumulative_infections_percentile_percentages_cell = cell(size(cumulative_infections_percentile_percentages,2),5);
    ages_split_cumulative_infections_percentile_percentages_cell = cell(size(ages_split_cumulative_infections_percentile_percentages,2),5);
    
    cumulative_infections_xmas_window_percentile_percentages_cell = cell(size(cumulative_infections_xmas_window_percentile_percentages,2),5);
    ages_split_cumul_infections_xmas_window_percentile_percent_cell = cell(size(ages_split_cumul_infections_xmas_window_percentile_percentages,2),5);
    
    % Initialise cells to write proportion results into
    incidence_27Dec2020_percentile_propns_cell = cell(size(incidence_27Dec2020_percentiles,2),5);
    incidence_26Dec2020_percentile_propns_cell = cell(size(incidence_26Dec2020_percentiles,2),5);
    incidence_25Dec2020_percentile_propns_cell = cell(size(incidence_25Dec2020_percentiles,2),5);
    incidence_24Dec2020_percentile_propns_cell = cell(size(incidence_24Dec2020_percentiles,2),5);
    incidence_23Dec2020_percentile_propns_cell = cell(size(incidence_23Dec2020_percentiles,2),5);
    
    ages_split_27Dec2020_incidence_propns_cell = cell(size(ages_split_27Dec2020_incidence_percentiles,2),5);
    ages_split_26Dec2020_incidence_propns_cell = cell(size(ages_split_26Dec2020_incidence_percentiles,2),5);
    ages_split_25Dec2020_incidence_propns_cell = cell(size(ages_split_25Dec2020_incidence_percentiles,2),5);
    ages_split_24Dec2020_incidence_propns_cell = cell(size(ages_split_24Dec2020_incidence_percentiles,2),5);
    ages_split_23Dec2020_incidence_propns_cell = cell(size(ages_split_23Dec2020_incidence_percentiles,2),5);
    
    cumulative_infections_percentile_propns_cell = cell(size(cumulative_infections_percentiles,2),5);
    ages_split_cumulative_infections_percentile_propns_cell = cell(size(ages_split_cumulative_infections_percentiles,2),5);
    
    cumulative_infections_xmas_window_percentile_propns_cell = cell(size(cumulative_infections_xmas_window_percentiles,2),5);
    ages_split_cumul_infections_xmas_window_percentile_propns_cell = cell(size(ages_split_cumulative_infections_xmas_window_percentiles,2),5);
    
    % Populate cell arrays
    % For each statistic, populates percentages cell array then proportions
    % cell array
    for slice_itr = 1:n_scens
    
        % Incidence on 27Dec2020
        for age_grp_itr = 1:size(incidence_27Dec2020_percentile_percentages,2)
            incidence_27Dec2020_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.3f%%(%0.3f%%,%0.3f%%)',...
                            incidence_27Dec2020_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            incidence_27Dec2020_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_27Dec2020_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             incidence_27Dec2020_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.5f(%0.5f,%0.5f)',...
                            incidence_27Dec2020_percentiles(2,age_grp_itr,slice_itr),... % Median
                            incidence_27Dec2020_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_27Dec2020_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Incidence on 26Dec2020
        for age_grp_itr = 1:size(incidence_26Dec2020_percentile_percentages,2)
            incidence_26Dec2020_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.3f%%(%0.3f%%,%0.3f%%)',...
                            incidence_26Dec2020_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            incidence_26Dec2020_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_26Dec2020_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             incidence_26Dec2020_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.5f(%0.5f,%0.5f)',...
                            incidence_26Dec2020_percentiles(2,age_grp_itr,slice_itr),... % Median
                            incidence_26Dec2020_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_26Dec2020_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Incidence on 25Dec2020
        for age_grp_itr = 1:size(incidence_25Dec2020_percentile_percentages,2)
            incidence_25Dec2020_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.3f%%(%0.3f%%,%0.3f%%)',...
                            incidence_25Dec2020_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            incidence_25Dec2020_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_25Dec2020_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             incidence_25Dec2020_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.5f(%0.5f,%0.5f)',...
                            incidence_25Dec2020_percentiles(2,age_grp_itr,slice_itr),... % Median
                            incidence_25Dec2020_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_25Dec2020_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Incidence on 24Dec2020
        for age_grp_itr = 1:size(incidence_24Dec2020_percentile_percentages,2)
            incidence_24Dec2020_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.3f%%(%0.3f%%,%0.3f%%)',...
                            incidence_24Dec2020_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            incidence_24Dec2020_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_24Dec2020_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             incidence_24Dec2020_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.5f(%0.5f,%0.5f)',...
                            incidence_24Dec2020_percentiles(2,age_grp_itr,slice_itr),... % Median
                            incidence_24Dec2020_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_24Dec2020_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Incidence on 23Dec2020
        for age_grp_itr = 1:size(incidence_23Dec2020_percentile_percentages,2)
            incidence_23Dec2020_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.3f%%(%0.3f%%,%0.3f%%)',...
                            incidence_23Dec2020_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            incidence_23Dec2020_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_23Dec2020_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             incidence_23Dec2020_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.5f(%0.5f,%0.5f)',...
                            incidence_23Dec2020_percentiles(2,age_grp_itr,slice_itr),... % Median
                            incidence_23Dec2020_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                            incidence_23Dec2020_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age group split of new cases on 27Dec2020
        for age_grp_itr = 1:size(ages_split_27Dec2020_incidence_percentile_percentages,2)
             ages_split_27Dec2020_incidence_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                            ages_split_27Dec2020_incidence_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            ages_split_27Dec2020_incidence_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_27Dec2020_incidence_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             ages_split_27Dec2020_incidence_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                            ages_split_27Dec2020_incidence_percentiles(2,age_grp_itr,slice_itr),... % Median
                            ages_split_27Dec2020_incidence_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_27Dec2020_incidence_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age group split of new cases on 26Dec2020
        for age_grp_itr = 1:size(ages_split_26Dec2020_incidence_percentile_percentages,2)
             ages_split_26Dec2020_incidence_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                            ages_split_26Dec2020_incidence_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            ages_split_26Dec2020_incidence_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_26Dec2020_incidence_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             ages_split_26Dec2020_incidence_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                            ages_split_26Dec2020_incidence_percentiles(2,age_grp_itr,slice_itr),... % Median
                            ages_split_26Dec2020_incidence_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_26Dec2020_incidence_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age group split of new cases on 25Dec2020
        for age_grp_itr = 1:size(ages_split_25Dec2020_incidence_percentile_percentages,2)
             ages_split_25Dec2020_incidence_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                            ages_split_25Dec2020_incidence_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            ages_split_25Dec2020_incidence_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_25Dec2020_incidence_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             ages_split_25Dec2020_incidence_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                            ages_split_25Dec2020_incidence_percentiles(2,age_grp_itr,slice_itr),... % Median
                            ages_split_25Dec2020_incidence_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_25Dec2020_incidence_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age group split of new cases on 24Dec2020
        for age_grp_itr = 1:size(ages_split_24Dec2020_incidence_percentile_percentages,2)
             ages_split_24Dec2020_incidence_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                            ages_split_24Dec2020_incidence_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            ages_split_24Dec2020_incidence_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_24Dec2020_incidence_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             ages_split_24Dec2020_incidence_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                            ages_split_24Dec2020_incidence_percentiles(2,age_grp_itr,slice_itr),... % Median
                            ages_split_24Dec2020_incidence_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_24Dec2020_incidence_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age group split of new cases on 23Dec2020
        for age_grp_itr = 1:size(ages_split_23Dec2020_incidence_percentile_percentages,2)
             ages_split_23Dec2020_incidence_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                            ages_split_23Dec2020_incidence_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                            ages_split_23Dec2020_incidence_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_23Dec2020_incidence_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
             ages_split_23Dec2020_incidence_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                            ages_split_23Dec2020_incidence_percentiles(2,age_grp_itr,slice_itr),... % Median
                            ages_split_23Dec2020_incidence_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                            ages_split_23Dec2020_incidence_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Cumulative infections over time horizon
        for age_grp_itr = 1:size(cumulative_infections_percentile_percentages,2)
            cumulative_infections_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f%%(%0.2f%%,%0.2f%%)',...
                cumulative_infections_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                cumulative_infections_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                cumulative_infections_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
            cumulative_infections_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.4f(%0.4f,%0.4f)',...
                cumulative_infections_percentiles(2,age_grp_itr,slice_itr),... % Median
                cumulative_infections_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                cumulative_infections_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Cumulative infections over Christmas bubble window
        for age_grp_itr = 1:size(cumulative_infections_xmas_window_percentile_percentages,2)
            cumulative_infections_xmas_window_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f%%(%0.2f%%,%0.2f%%)',...
                cumulative_infections_xmas_window_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                cumulative_infections_xmas_window_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                cumulative_infections_xmas_window_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
            cumulative_infections_xmas_window_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.4f(%0.4f,%0.4f)',...
                cumulative_infections_xmas_window_percentiles(2,age_grp_itr,slice_itr),... % Median
                cumulative_infections_xmas_window_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                cumulative_infections_xmas_window_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age split of cumulative infections over time horizon
        for age_grp_itr = 1:size(ages_split_cumulative_infections_percentile_percentages,2)
            ages_split_cumulative_infections_percentile_percentages_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                ages_split_cumulative_infections_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                ages_split_cumulative_infections_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                ages_split_cumulative_infections_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
            ages_split_cumulative_infections_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                ages_split_cumulative_infections_percentiles(2,age_grp_itr,slice_itr),... % Median
                ages_split_cumulative_infections_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                ages_split_cumulative_infections_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    
        % Age split of cumulative infections over Christmas bubble window
        for age_grp_itr = 1:size(ages_split_cumul_infections_xmas_window_percentile_percentages,2)
            ages_split_cumul_infections_xmas_window_percentile_percent_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.0f%%(%0.0f%%,%0.0f%%)',...
                ages_split_cumul_infections_xmas_window_percentile_percentages(2,age_grp_itr,slice_itr),... % Median
                ages_split_cumul_infections_xmas_window_percentile_percentages(1,age_grp_itr,slice_itr),... % Lower uncertainty bound 
                ages_split_cumul_infections_xmas_window_percentile_percentages(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
    
            ages_split_cumul_infections_xmas_window_percentile_propns_cell{age_grp_itr,slice_itr} =...
                sprintf('%0.2f(%0.2f,%0.2f)',...
                ages_split_cumulative_infections_xmas_window_percentiles(2,age_grp_itr,slice_itr),... % Median
                ages_split_cumulative_infections_xmas_window_percentiles(1,age_grp_itr,slice_itr),... % Lower uncertainty bound
                ages_split_cumulative_infections_xmas_window_percentiles(3,age_grp_itr,slice_itr));   % Upper uncertainty bound
        end
    end

    % Consolidate summary stats into output cells.
    incidence_percentile_percentages_output =...
        {incidence_23Dec2020_percentile_percentages_cell,...
        incidence_24Dec2020_percentile_percentages_cell,...
        incidence_25Dec2020_percentile_percentages_cell,...
        incidence_26Dec2020_percentile_percentages_cell,...
        incidence_27Dec2020_percentile_percentages_cell};

    age_split_incidence_percentile_percentages_output =...
        {ages_split_23Dec2020_incidence_percentages_cell,...
        ages_split_24Dec2020_incidence_percentages_cell,...
        ages_split_25Dec2020_incidence_percentages_cell,...
        ages_split_26Dec2020_incidence_percentages_cell,...
        ages_split_27Dec2020_incidence_percentages_cell};

    cumulative_infections_percentile_percentages_output =...
        {cumulative_infections_percentile_percentages_cell,...
        cumulative_infections_xmas_window_percentile_percentages_cell};

    age_split_cumulative_infections_percentile_percentages_output =...
        {ages_split_cumulative_infections_percentile_percentages_cell,...
        ages_split_cumul_infections_xmas_window_percentile_percent_cell};

    %% Generate figures if requested
    if generate_figures_flag == true
    
        %% Produce violin plot with grouping by age group and shading by scenario
        %% Age distribution of new latent infecteds on given day
        
        % Set x-positions data will be plotted at
        x_plot_pos = [1:5;7:11;13:17;19:23];
        
        % Set colour values for each plotting group
        colour_vals = [0, 0.4470, 0.7410;
                        0.8500, 0.3250, 0.0980;
                        0.9290, 0.6940, 0.1250;
                        0.4940, 0.1840, 0.5560;
                        0.4660, 0.6740, 0.1880;
                        0.3010, 0.7450, 0.9330;
                        0.5, 0.5, 0.5];
                    
        % Set shading intensity of each violin
        alpha_vec = [1. 0.8 0.6 0.4 0.2];
                    
        % Set up legend labels
        legend_label = {'Scenario A', 'Scenario B','Scenario C','Scenario D','Scenario E'};
        
        % Set x-axis labels
        xaxis_label = 'Age group';
        xticks_vals = [3 9 15 21];
        xticks_labels = {'0-19yrs','20-64yrs','65+yrs','All'};
        
        % Set y-axis label
        if plot_percentages_flag == true
            yaxis_label_vec = {'Percentage newly infected on 23 December 2020 (%)';...
                                'Percentage newly infected on 24 December 2020 (%)';...
                                'Percentage newly infected on 25 December 2020 (%)';...
                                'Percentage newly infected on 26 December 2020 (%)';...
                                'Percentage newly infected on 27 December 2020 (%)'};
        else
            yaxis_label_vec = {'Proportion newly infected on 23 December 2020';...
                                'Proportion newly infected on 24 December 2020';...
                                'Proportion newly infected on 25 December 2020';...
                                'Proportion newly infected on 26 December 2020';...
                                'Proportion newly infected on 27 December 2020'};
        end
        
        %Set axes limits
        x_limits = [0 24];
        if main_analysis_flag == true
            if plot_percentages_flag == true
                y_limits = [0.0 0.21];
            else
                y_limits = [0.0 0.0021];
            end
        else
            if plot_percentages_flag == true
                y_limits = [0.0 0.1];
            else
                y_limits = [0.0 0.001];
            end
        end
        
        % Set plot fontsize
        plot_fontsize = 22;   
        
        %%
        
        % Set up plot data
        input_data_incidence = cell(5,1);
        
        % Bring together age-group data and all ages data
        scen_A_array_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_support_bubble new_latent_27Dec2020_propn_overall_popn_support_bubble];
        scen_B_array_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_short_faithful_bubble new_latent_27Dec2020_propn_overall_popn_short_faithful_bubble];
        scen_C_array_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_faithful_bubble new_latent_27Dec2020_propn_overall_popn_faithful_bubble];
        scen_D_array_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_fixed_bubble new_latent_27Dec2020_propn_overall_popn_fixed_bubble];
        scen_E_array_27Dec2020 = [new_latent_27Dec2020_propn_agegrp_unfaithful_bubble new_latent_27Dec2020_propn_overall_popn_unfaithful_bubble];
        
        scen_A_array_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_support_bubble new_latent_26Dec2020_propn_overall_popn_support_bubble];
        scen_B_array_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_short_faithful_bubble new_latent_26Dec2020_propn_overall_popn_short_faithful_bubble];
        scen_C_array_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_faithful_bubble new_latent_26Dec2020_propn_overall_popn_faithful_bubble];
        scen_D_array_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_fixed_bubble new_latent_26Dec2020_propn_overall_popn_fixed_bubble];
        scen_E_array_26Dec2020 = [new_latent_26Dec2020_propn_agegrp_unfaithful_bubble new_latent_26Dec2020_propn_overall_popn_unfaithful_bubble];
        
        scen_A_array_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_support_bubble new_latent_25Dec2020_propn_overall_popn_support_bubble];
        scen_B_array_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_short_faithful_bubble new_latent_25Dec2020_propn_overall_popn_short_faithful_bubble];
        scen_C_array_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_faithful_bubble new_latent_25Dec2020_propn_overall_popn_faithful_bubble];
        scen_D_array_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_fixed_bubble new_latent_25Dec2020_propn_overall_popn_fixed_bubble];
        scen_E_array_25Dec2020 = [new_latent_25Dec2020_propn_agegrp_unfaithful_bubble new_latent_25Dec2020_propn_overall_popn_unfaithful_bubble];
        
        scen_A_array_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_support_bubble new_latent_24Dec2020_propn_overall_popn_support_bubble];
        scen_B_array_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_short_faithful_bubble new_latent_24Dec2020_propn_overall_popn_short_faithful_bubble];
        scen_C_array_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_faithful_bubble new_latent_24Dec2020_propn_overall_popn_faithful_bubble];
        scen_D_array_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_fixed_bubble new_latent_24Dec2020_propn_overall_popn_fixed_bubble];
        scen_E_array_24Dec2020 = [new_latent_24Dec2020_propn_agegrp_unfaithful_bubble new_latent_24Dec2020_propn_overall_popn_unfaithful_bubble];
        
        scen_A_array_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_support_bubble new_latent_23Dec2020_propn_overall_popn_support_bubble];
        scen_B_array_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_short_faithful_bubble new_latent_23Dec2020_propn_overall_popn_short_faithful_bubble];
        scen_C_array_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_faithful_bubble new_latent_23Dec2020_propn_overall_popn_faithful_bubble];
        scen_D_array_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_fixed_bubble new_latent_23Dec2020_propn_overall_popn_fixed_bubble];
        scen_E_array_23Dec2020 = [new_latent_23Dec2020_propn_agegrp_unfaithful_bubble new_latent_23Dec2020_propn_overall_popn_unfaithful_bubble];
        
        % Add to input data cells for plotting
        input_data_incidence{1} = cat(3,scen_A_array_23Dec2020,...
                                    scen_B_array_23Dec2020,...
                                    scen_C_array_23Dec2020,...
                                    scen_D_array_23Dec2020,...
                                    scen_E_array_23Dec2020);
        
        input_data_incidence{2} = cat(3,scen_A_array_24Dec2020,...
                                    scen_B_array_24Dec2020,...
                                    scen_C_array_24Dec2020,...
                                    scen_D_array_24Dec2020,...
                                    scen_E_array_24Dec2020);
        
        input_data_incidence{3} = cat(3,scen_A_array_25Dec2020,...
                                    scen_B_array_25Dec2020,...
                                    scen_C_array_25Dec2020,...
                                    scen_D_array_25Dec2020,...
                                    scen_E_array_25Dec2020);
        
        input_data_incidence{4} = cat(3,scen_A_array_26Dec2020,...
                                    scen_B_array_26Dec2020,...
                                    scen_C_array_26Dec2020,...
                                    scen_D_array_26Dec2020,...
                                    scen_E_array_26Dec2020);
        
        input_data_incidence{5} = cat(3,scen_A_array_27Dec2020,...
                                    scen_B_array_27Dec2020,...
                                    scen_C_array_27Dec2020,...
                                    scen_D_array_27Dec2020,...
                                    scen_E_array_27Dec2020);
        
        
        % Set save filename
        if main_analysis_flag == true
            if plot_percentages_flag == true
                save_filename = 'figure_files/new_latent_age_propn_of_age_grps_violins_percentages_main_analysis';
            else
                save_filename = 'figure_files/new_latent_age_propn_of_age_grps_violins_main_analysis';
            end
        else
            if plot_percentages_flag == true
                save_filename = 'figure_files/new_latent_age_propn_of_age_grps_violins_percentages_alternative_analysis';
            else
                save_filename = 'figure_files/new_latent_age_propn_of_age_grps_violins_alternative_analysis';
            end
        end
        
        % Amend save_filename if using lower adherence data
        if main_adherence_flag == false
            save_filename = strcat(save_filename,'_lower_adherence');
        end
        
        % Generate the plots
        scenario_label_x_offset = 0.25; % Variable used for placing scenario text label above violin plots  
        for plot_itr = 1:5
            % Add calendar date suffix to save_filename
            save_filename_final = strcat(save_filename,save_filename_calender_date_suffixes{plot_itr});
        
            % Generate figure
            if plot_percentages_flag == true
                % Multiply input data by 100, to plot as percentages
                plot_violin_panels_group_by_age(input_data_incidence{plot_itr}*100,...
                                                        x_plot_pos,...
                                                        colour_vals,...
                                                        alpha_vec,...
                                                        legend_label,...
                                                        xaxis_label,...
                                                        yaxis_label_vec{plot_itr},...
                                                        xticks_vals,...
                                                        xticks_labels,...
                                                        x_limits,...
                                                        y_limits,...
                                                        scenario_label_x_offset,...
                                                        plot_fontsize,...
                                                        save_filename_final,...
                                                        save_file_flag)     
            else
                plot_violin_panels_group_by_age(input_data_incidence{plot_itr},...
                                                        x_plot_pos,...
                                                        colour_vals,...
                                                        alpha_vec,...
                                                        legend_label,...
                                                        xaxis_label,...
                                                        yaxis_label_vec{plot_itr},...
                                                        xticks_vals,...
                                                        xticks_labels,...
                                                        x_limits,...
                                                        y_limits,...
                                                        scenario_label_x_offset,...
                                                        plot_fontsize,...
                                                        save_filename_final,...
                                                        save_file_flag)
            end
        end
        
        %% Produce violin plot with grouping by age group and shading by scenario
        %% Age distribution of new infections
        
        % Set x-positions data will be plotted at
        x_plot_pos = [1:5;7:11;13:17];
        
        % Set colour values for each plotting group
        colour_vals = [0, 0.4470, 0.7410;
                        0.8500, 0.3250, 0.0980;
                        0.9290, 0.6940, 0.1250;
                        0.4940, 0.1840, 0.5560;
                        0.4660, 0.6740, 0.1880;
                        0.3010, 0.7450, 0.9330;
                        0.5, 0.5, 0.5];
                    
        % Set shading intensity of each violin
        alpha_vec = [1. 0.8 0.6 0.4 0.2];
                    
        % Set up legend labels
        legend_label = {'Scenario A','Scenario B','Scenario C','Scenario D','Scenario E'};
        
        % Set x-axis labels
        xaxis_label = 'Age group';
        xticks_vals = [3 9 15];
        xticks_labels = {'0-19yrs','20-64yrs','65+yrs'};
        
        % Set y-axis label
        if plot_percentages_flag == true
            yaxis_label_vec = {'Share of new infections between 23-27 December 2020 (%)';...
                            'Share of new infections over time horizon (%)'};
        else
            yaxis_label_vec = {'Share of new infections between 23-27 December 2020';...
                            'Share of new infections over time horizon'};
        end
        
        %Set axes limits
        x_limits = [0 18];
        
        if plot_percentages_flag == true
            y_limits = [0 60];
        else
            y_limits = [0 0.6];
        end
        
        % Set plot fontsize
        plot_fontsize = 22;   
        
        % Set up plot data
        input_data_share_of_new_infections = cell(2,1);
        
        input_data_share_of_new_infections{1} = cat(3,scen_A_age_split_cumulative_infections_xmas_window,...
                                                scen_B_age_split_cumulative_infections_xmas_window,...
                                                scen_C_age_split_cumulative_infections_xmas_window,...
                                                scen_D_age_split_cumulative_infections_xmas_window,...
                                                scen_E_age_split_cumulative_infections_xmas_window);
        
        input_data_share_of_new_infections{2} = cat(3,scen_A_age_split_cumulative_infections,...
                                                scen_B_age_split_cumulative_infections,...
                                                scen_C_age_split_cumulative_infections,...
                                                scen_D_age_split_cumulative_infections,...
                                                scen_E_age_split_cumulative_infections);
        
        % Set save filename
        if main_analysis_flag == true
            if plot_percentages_flag == true
                save_filename = 'figure_files/new_infections_age_breakdown_violins_percentages_main_analysis';
            else
                save_filename = 'figure_files/new_infections_age_breakdown_violins_main_analysis';
            end
        else
            if plot_percentages_flag == true
                save_filename = 'figure_files/new_infections_age_breakdown_violins_percentages_alternative_analysis';
            else
                save_filename = 'figure_files/new_infections_age_breakdown_violins_alternative_analysis';
            end
        end
        
        % Amend save_filename if using lower adherence data
        if main_adherence_flag == false
            save_filename = strcat(save_filename,'_lower_adherence');
        end
        
        % Generate the plots
        scenario_label_x_offset = 0.15; % Variable used for placing scenario text label above violin plots  
        for plot_itr = 1:length(input_data_share_of_new_infections)
            % Add calendar date suffix to save_filename
            save_filename_final = strcat(save_filename,save_filename_prevalence_suffixes{plot_itr});
        
            % Generate figure
            if plot_percentages_flag == true
                % Multiply input data by 100, to plot as percentages
                plot_violin_panels_group_by_age(input_data_share_of_new_infections{plot_itr}*100,...
                                                    x_plot_pos,...
                                                    colour_vals,...
                                                    alpha_vec,...
                                                    legend_label,...
                                                    xaxis_label,...
                                                    yaxis_label_vec{plot_itr},...
                                                    xticks_vals,...
                                                    xticks_labels,...
                                                    x_limits,...
                                                    y_limits,...
                                                    scenario_label_x_offset,...
                                                    plot_fontsize,...
                                                    save_filename_final,...
                                                    save_file_flag)
            else
                plot_violin_panels_group_by_age(input_data_share_of_new_infections{plot_itr},...
                                                    x_plot_pos,...
                                                    colour_vals,...
                                                    alpha_vec,...
                                                    legend_label,...
                                                    xaxis_label,...
                                                    yaxis_label_vec{plot_itr},...
                                                    xticks_vals,...
                                                    xticks_labels,...
                                                    x_limits,...
                                                    y_limits,...
                                                    scenario_label_x_offset,...
                                                    plot_fontsize,...
                                                    save_filename_final,...
                                                    save_file_flag)
            end
        end
        
        
        %% Produce violin plot with grouping by age group and shading by scenario
        %% Total new infections over the time horizon as a proportion of specified groups
        
        % Set x-positions data will be plotted at
        x_plot_pos = [1:5;7:11;13:17;19:23];
        
        % Set colour values for each plotting group
        colour_vals = [0, 0.4470, 0.7410;
                        0.8500, 0.3250, 0.0980;
                        0.9290, 0.6940, 0.1250;
                        0.4940, 0.1840, 0.5560;
                        0.4660, 0.6740, 0.1880;
                        0.3010, 0.7450, 0.9330;
                        0.5, 0.5, 0.5];
                    
        % Set shading intensity of each violin
        alpha_vec = [1. 0.8 0.6 0.4 0.2];
                    
        % Set up legend labels
        legend_label = {'Scenario A', 'Scenario B','Scenario C','Scenario D','Scenario E'};
        
        % Set x-axis labels
        xaxis_label = 'Age group';
        xticks_vals = [3 9 15 21];
        xticks_labels = {'0-19yrs','20-64yrs','65+yrs','All'};
        
        % Set y-axis label
        if plot_percentages_flag == true
            yaxis_label_vec = {'Cumulative infection between 23-27 December 2020 (%)';...
                            'Cumulative infection over time horizon (%)'};
        else
            yaxis_label_vec = {'Cumulative fraction newly infected between 23-27 December 2020';...
                            'Cumulative fraction newly infected over time horizon'};
        end
        
        %Set axes limits
        x_limits = [0 24];
        
        if plot_percentages_flag == true
            y_limits = [0.0 1.4];
        else
            y_limits = [0.0 0.014];
        end
        
        % Set plot fontsize
        plot_fontsize = 22;   
        
        % Collect together relevant data
        input_data_cumulative_infections = cell(2,1);
        
        scen_A_array_cumulative_infections_23Dec2020_27Dec2020 = [newinf_xmas_bubble_window_propn_agegrp_support_bubble newinf_xmas_bubble_window_propn_overall_popn_support_bubble];
        scen_B_array_cumulative_infections_23Dec2020_27Dec2020 = [newinf_xmas_bubble_window_propn_agegrp_short_bubble newinf_xmas_bubble_window_propn_overall_popn_short_bubble];
        scen_C_array_cumulative_infections_23Dec2020_27Dec2020 = [newinf_xmas_bubble_window_propn_agegrp_faithful_bubble newinf_xmas_bubble_window_propn_overall_popn_faithful_bubble];
        scen_D_array_cumulative_infections_23Dec2020_27Dec2020 = [newinf_xmas_bubble_window_propn_agegrp_fixed_bubble newinf_xmas_bubble_window_propn_overall_popn_fixed_bubble];
        scen_E_array_cumulative_infections_23Dec2020_27Dec2020 = [newinf_xmas_bubble_window_propn_agegrp_unfaithful_bubble newinf_xmas_bubble_window_propn_overall_popn_unfaithful_bubble];
        
        input_data_cumulative_infections{1} = cat(3,scen_A_array_cumulative_infections_23Dec2020_27Dec2020,...
                                        scen_B_array_cumulative_infections_23Dec2020_27Dec2020,...
                                        scen_C_array_cumulative_infections_23Dec2020_27Dec2020,...
                                        scen_D_array_cumulative_infections_23Dec2020_27Dec2020,...
                                        scen_E_array_cumulative_infections_23Dec2020_27Dec2020);
        
        scen_A_array_cumulative_infections_time_horizon = [newinf_propn_agegrp_support_bubble newinf_propn_overall_popn_support_bubble];
        scen_B_array_cumulative_infections_time_horizon = [newinf_propn_agegrp_short_faithful_bubble newinf_propn_overall_popn_short_faithful_bubble];
        scen_C_array_cumulative_infections_time_horizon = [newinf_propn_agegrp_faithful_bubble newinf_propn_overall_popn_faithful_bubble];
        scen_D_array_cumulative_infections_time_horizon = [newinf_propn_agegrp_fixed_bubble newinf_propn_overall_popn_fixed_bubble];
        scen_E_array_cumulative_infections_time_horizon = [newinf_propn_agegrp_unfaithful_bubble newinf_propn_overall_popn_unfaithful_bubble];
        
        input_data_cumulative_infections{2} = cat(3,scen_A_array_cumulative_infections_time_horizon,...
                                        scen_B_array_cumulative_infections_time_horizon,...
                                        scen_C_array_cumulative_infections_time_horizon,...
                                        scen_D_array_cumulative_infections_time_horizon,...
                                        scen_E_array_cumulative_infections_time_horizon);
        
        
        % Set save filename
        if main_analysis_flag == true
            if plot_percentages_flag == true
                save_filename = 'figure_files/newinf_percentage_age_grps_violins_main_analysis';
            else
                save_filename = 'figure_files/newinf_propn_of_age_grps_violins_main_analysis';
            end
        else
            if plot_percentages_flag == true
                save_filename = 'figure_files/newinf_percentage_age_grps_violins_alternative_analysis';
            else
                save_filename = 'figure_files/newinf_propn_of_age_grps_violins_alternative_analysis';
            end
        end
        
        % Amend save_filename if using lower adherence data
        if main_adherence_flag == false
            save_filename = strcat(save_filename,'_lower_adherence');
        end
        
        % Generate the figures
        scenario_label_x_offset = 0.15; % Variable used for placing scenario text label above violin plots  
        for plot_itr = 1:length(input_data_cumulative_infections)
            % Add calendar date suffix to save_filename
            save_filename_final = strcat(save_filename,save_filename_prevalence_suffixes{plot_itr});
        
            if plot_percentages_flag == true
                % Multiply input data by 100, to plot as percentages
                plot_violin_panels_group_by_age(input_data_cumulative_infections{plot_itr}*100,...
                                                    x_plot_pos,...
                                                    colour_vals,...
                                                    alpha_vec,...
                                                    legend_label,...
                                                    xaxis_label,...
                                                    yaxis_label_vec{plot_itr},...
                                                    xticks_vals,...
                                                    xticks_labels,...
                                                    x_limits,...
                                                    y_limits,...
                                                    scenario_label_x_offset,...
                                                    plot_fontsize,...
                                                    save_filename_final,...
                                                    save_file_flag)
            else
                plot_violin_panels_group_by_age(input_data_cumulative_infections{plot_itr},...
                                                    x_plot_pos,...
                                                    colour_vals,...
                                                    alpha_vec,...
                                                    legend_label,...
                                                    xaxis_label,...
                                                    yaxis_label_vec{plot_itr},...
                                                    xticks_vals,...
                                                    xticks_labels,...
                                                    x_limits,...
                                                    y_limits,...
                                                    scenario_label_x_offset,...
                                                    plot_fontsize,...
                                                    save_filename_final,...
                                                    save_file_flag)
            end
        end
    end
end

%% Plotting function
function plot_violin_panels_group_by_age(input_data,...
                                            x_plot_pos,...
                                            colour_vals,...
                                            alpha_vec,...
                                            legend_label,...
                                            xaxis_label,...
                                            yaxis_label,...
                                            xticks_vals,...
                                            xticks_labels,...
                                            x_limits,...
                                            y_limits,...
                                            scenario_label_x_offset,...
                                            plot_fontsize,...
                                            save_filename, ...
                                            save_file_flag)

  % Initialise figures
    position = [10, 10, 1.5*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    %set(0, 'DefaultFigurePosition', position,'Units','points');
    fig = figure();
    clf;
    set(fig,'Color', [1 1 1])
    hold on
    
    % Set text labels for scenarios to be put above each violin plot
    scenario_text_label = {'A','B','C','D','E'};
  
    % Get number of groups of data in use
    n_data_grps = size(input_data,2);
    
    % Get number of scenarios in use. Number of slices of array
    n_scens = size(input_data,3);
    
    % Pick out relevant data for this statistic
    stat_data = input_data;
    for group_bin_itr = 1:n_data_grps
        % Iterate over each non-baseline scenario. Plot data for each
        % scenario
        group_data = stat_data(:,group_bin_itr,:);

        % Get y-position to put scenario labels
        % Go from largest value cross scenario in current age group and 
        % move up by 5% of the span of the y-axis
        scenario_text_label_y_pos = y_limits(2)*0.95;
        %scenario_text_label_y_pos = max(group_data(:)) + (y_limits(2)*0.05);

        for scen_itr = 1:n_scens
            % Get group data for this statistic to be plotted
            config_data = group_data(:,scen_itr);

            % Get colour for the plot
            colour_vec = colour_vals(group_bin_itr,:);
            
            % Set shading intensity
            alpha_val = alpha_vec(scen_itr);
            
            %Create violin plot for each data group
            % Otherwise, create violin plot
            violins = Violin(config_data, x_plot_pos(group_bin_itr,scen_itr),...
                                'ShowData',false);
            
            % Set violin plot properties
            violins.ViolinColor = colour_vec;
            violins.ViolinAlpha = alpha_val; %shading transparency
            
            % Set violin plot region properties
            violins.EdgeColor = [1 1 1];
            
            % Set median marker properties
            violins.MedianColor = [1 1 1];
            violins.MedianPlot.Marker = 's';
            violins.MedianPlot.MarkerEdgeColor = [0 0 0];
            violins.MedianPlot.LineWidth = 1;
            
            % Set whisker line properties
            violins.BoxColor = [0 0 0];

            % Add a label for the scenario above the violin plot
            t = text(x_plot_pos(group_bin_itr,scen_itr)-scenario_label_x_offset,scenario_text_label_y_pos,scenario_text_label{scen_itr});
            t.FontSize = 16;
            t.FontWeight = 'bold';
        end
    end
    
%     % Add line at y=1
%     plot([0 x_limits(end)],[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
%     
    % Set axes limits
    xlim(x_limits)
    ylim(y_limits)
    
    % Add x-axis label, if required
    xlabel(xaxis_label)
    
    % Add xtick labels
    xticks(xticks_vals)
    xticklabels(xticks_labels)
    
    % Add y-axis labels
    ylabel(yaxis_label)

    %Specify general axis properties
    set(gca,'FontSize',plot_fontsize)
    set(gca,'LineWidth',1)
    box on
      
    if save_file_flag == true
        % Save figure to file  
        export_fig(save_filename,'-pdf','-r1200')
    end
    
%     % Add title
%     title(title_vec{stat_itr})

%     % Add legend
%     Initialise patch vector and populate
%     H = gobjects(n_scens,1);
%     for scen_itr = 1:n_scens
%         H(scen_itr) = fill(NaN,NaN,[0 0 0],'facealpha', alpha_vec(scen_itr),'DisplayName',legend_label{scen_itr});
%     end
% 
%     Construct the legend
%     legend(H,...
%         'LineWidth',1.5,...
%         'FontSize',plot_fontsize)   
end

% %% Produce violin plot with grouping by age group and shading by scenario
% %% Age distribution of new latent infects on stated calendar date
% 
% % Set x-positions data will be plotted at
% x_plot_pos = [1:5;7:11;13:17];
% 
% % Set colour values for each plotting group
% colour_vals = [0, 0.4470, 0.7410;
%                 0.8500, 0.3250, 0.0980;
%                 0.9290, 0.6940, 0.1250;
%                 0.4940, 0.1840, 0.5560;
%                 0.4660, 0.6740, 0.1880;
%                 0.3010, 0.7450, 0.9330;
%                 0.5, 0.5, 0.5];
%             
% % Set shading intensity of each violin
% alpha_vec = [1. 0.8 0.6 0.4 0.2];
%             
% % Set up legend labels
% legend_label = {'Scenario A','Scenario B','Scenario C','Scenario D','Scenario E'};
% 
% % Set x-axis labels
% xaxis_label = 'Age group';
% xticks_vals = [3 9 15];
% xticks_labels = {'0-19yrs','20-64yrs','65+yrs'};
% 
% % Set y-axis label
% if plot_percentages_flag == true
%     yaxis_label_vec = {'Percentage of those infected on 23 December 2020 (%)';...
%                     'Percentage of those infected on 24 December 2020 (%)';...
%                     'Percentage of those infected on 25 December 2020 (%)';...
%                     'Percentage of those infected on 26 December 2020 (%)';...
%                     'Percentage of those infected on 27 December 2020 (%)'};
% else
%     yaxis_label_vec = {'Proportion of those infected on 23 December 2020 (%)';...
%                     'Proportion of those infected on 24 December 2020 (%)';...
%                     'Proportion of those infected on 25 December 2020 (%)';...
%                     'Proportion of those infected on 26 December 2020 (%)';...
%                     'Proportion of those infected on 27 December 2020 (%)'};
% end
% 
% %Set axes limits
% x_limits = [0 18];
% 
% if plot_percentages_flag == true
%     y_limits = [0 100];
% else
%     y_limits = [0 1];
% end
% 
% % Set plot fontsize
% plot_fontsize = 22;   
% 
% % Set up plot data
% input_data = cell(5,1);
% 
% input_data{1} = cat(3,new_latent_23Dec2020_support_bubble_only_normalised,...
%                     new_latent_23Dec2020_short_faithful_bubble_normalised,...
%                     new_latent_23Dec2020_faithful_bubble_normalised,...
%                     new_latent_23Dec2020_fixed_bubble_normalised,...
%                     new_latent_23Dec2020_unfaithful_bubble_normalised);
% 
% input_data{2} = cat(3,new_latent_24Dec2020_support_bubble_only_normalised,...
%                     new_latent_24Dec2020_short_faithful_bubble_normalised,...
%                     new_latent_24Dec2020_faithful_bubble_normalised,...
%                     new_latent_24Dec2020_fixed_bubble_normalised,...
%                     new_latent_24Dec2020_unfaithful_bubble_normalised);
% 
% input_data{3} = cat(3,new_latent_25Dec2020_support_bubble_only_normalised,...
%                     new_latent_25Dec2020_short_faithful_bubble_normalised,...
%                     new_latent_25Dec2020_faithful_bubble_normalised,...
%                     new_latent_25Dec2020_fixed_bubble_normalised,...
%                     new_latent_25Dec2020_unfaithful_bubble_normalised);
% 
% input_data{4} = cat(3,new_latent_26Dec2020_support_bubble_only_normalised,...
%                     new_latent_26Dec2020_short_faithful_bubble_normalised,...
%                     new_latent_26Dec2020_faithful_bubble_normalised,...
%                     new_latent_26Dec2020_fixed_bubble_normalised,...
%                     new_latent_26Dec2020_unfaithful_bubble_normalised);
% 
% input_data{5} = cat(3,new_latent_27Dec2020_support_bubble_only_normalised,...
%                     new_latent_27Dec2020_short_faithful_bubble_normalised,...
%                     new_latent_27Dec2020_faithful_bubble_normalised,...
%                     new_latent_27Dec2020_fixed_bubble_normalised,...
%                     new_latent_27Dec2020_unfaithful_bubble_normalised);
% 
% %%
% % Set save filename
% if main_analysis_flag == true
%     if plot_percentages_flag == true
%         save_filename = 'figure_files/new_latent_age_violins_percentages_main_analysis';
%     else
%         save_filename = 'figure_files/new_latent_age_violins_main_analysis';
%     end
% else
%     if plot_percentages_flag == true
%         save_filename = 'figure_files/new_latent_age_violins_percentages_alternative_analysis';
%     else
%         save_filename = 'figure_files/new_latent_age_violins_alternative_analysis';
%     end
% end
% 
% % Amend save_filename if using lower adherence data
% if main_adherence_flag == false
%     save_filename = strcat(save_filename,'_lower_adherence');
% end
% 
% for plot_itr = 1:5
%     % Add calendar date suffix to save_filename
%     save_filename_final = strcat(save_filename,save_filename_calender_date_suffixes{plot_itr});
% 
%     % Generate figure
%     if plot_percentages_flag == true
%         % Multiply input data by 100, to plot as percentages
%         plot_violin_panels_group_by_age(input_data{plot_itr}*100,...
%                                             x_plot_pos,...
%                                             colour_vals,...
%                                             alpha_vec,...
%                                             legend_label,...
%                                             xaxis_label,...
%                                             yaxis_label_vec{plot_itr},...
%                                             xticks_vals,...
%                                             xticks_labels,...
%                                             x_limits,...
%                                             y_limits,...
%                                             plot_fontsize,...
%                                             save_filename_final)
%     else
%         plot_violin_panels_group_by_age(input_data{plot_itr},...
%                                             x_plot_pos,...
%                                             colour_vals,...
%                                             alpha_vec,...
%                                             legend_label,...
%                                             xaxis_label,...
%                                             yaxis_label_vec{plot_itr},...
%                                             xticks_vals,...
%                                             xticks_labels,...
%                                             x_limits,...
%                                             y_limits,...
%                                             plot_fontsize,...
%                                             save_filename_final)
%     end
% end