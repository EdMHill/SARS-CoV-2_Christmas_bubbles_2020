%Purpose:
% Generate summary staistics and plot for the different time horizons and 
% adherence assumptions
%--------------------------------------------------------------------------

clear variables

%% ADD PATH DEPENDENCIES
addpath('matlab_packages/Violinplot-Matlab')
addpath('matlab_packages/export_fig')

%% Flag options for this file

% Specify if figures should be produced
generate_figures_flag = true;

% Specify if generated plots should be saved to file
save_file_flag = false;

% Specify if should plot percentages (flag value = true) or proportions (flag value = false)
plot_percentages_flag = true;

%% Specify the scenarios to be analysed and type of outputs to be produced

%% Main analysis, 70% adherence

% Specify if analysing main scenario runs (flag value = true) or the
% sensitivity analysis with longer time horizon (flag value = false)
main_analysis_flag = true;

% Specify if looking at main adherence of probability of 70% (flag value = true) 
% or lower adherence run with probability of 30% (flag value = false)
main_adherence_flag = true;

% Call function to produce plots and return summary statistics
[incidence_percentile_percentages_output_main_analysis,...
     age_split_incidence_percentile_percentages_output_main_analysis,...
     cumul_infections_percentile_percentages_output_main_analysis,...
     age_split_cumul_infections_percentile_output_main_analysis] = ...
                generate_summ_stats_and_plots_fn(main_analysis_flag,...
                                                main_adherence_flag,...
                                                plot_percentages_flag,...
                                                generate_figures_flag,...
                                                save_file_flag);

% Unpack cell contents
incidence_23Dec2020_main_analysis = incidence_percentile_percentages_output_main_analysis{1};
incidence_24Dec2020_main_analysis = incidence_percentile_percentages_output_main_analysis{2};
incidence_25Dec2020_main_analysis = incidence_percentile_percentages_output_main_analysis{3};
incidence_26Dec2020_main_analysis = incidence_percentile_percentages_output_main_analysis{4};
incidence_27Dec2020_main_analysis = incidence_percentile_percentages_output_main_analysis{5};

age_split_incidence_23Dec2020_main_analysis = age_split_incidence_percentile_percentages_output_main_analysis{1};
age_split_incidence_24Dec2020_main_analysis = age_split_incidence_percentile_percentages_output_main_analysis{2};
age_split_incidence_25Dec2020_main_analysis = age_split_incidence_percentile_percentages_output_main_analysis{3};
age_split_incidence_26Dec2020_main_analysis = age_split_incidence_percentile_percentages_output_main_analysis{4};
age_split_incidence_27Dec2020_main_analysis = age_split_incidence_percentile_percentages_output_main_analysis{5};

cumulative_infections_time_horizon_main_analysis = cumul_infections_percentile_percentages_output_main_analysis{1};
cumulative_infections_23Dec2020_27Dec2020_main_analysis = cumul_infections_percentile_percentages_output_main_analysis{2};

age_split_cumulative_infections_time_horizon_main_analysis = age_split_cumul_infections_percentile_output_main_analysis{1};
age_split_cumul_infections_23Dec2020_27Dec2020_main_analysis = age_split_cumul_infections_percentile_output_main_analysis{2};

%% Alternative analysis, 70% adherence

% Specify if analysing main scenario runs (flag value = true) or the
% sensitivity analysis with longer time horizon (flag value = false)
main_analysis_flag = false;

% Specify if looking at main adherence of probability of 70% (flag value = true) 
% or lower adherence run with probability of 30% (flag value = false)
main_adherence_flag = true;

% Call function to produce plots and return summary statistics
[incidence_percentile_percentages_output_alt_analysis,...
     age_split_incidence_percentile_percentages_output_alt_analysis,...
     cumul_infections_percentile_percentages_output_alt_analysis,...
     age_split_cumul_infections_percentile_output_alt_analysis] = ...
                generate_summ_stats_and_plots_fn(main_analysis_flag,...
                                                main_adherence_flag,...
                                                plot_percentages_flag,...
                                                generate_figures_flag,...
                                                save_file_flag);

% Unpack cell contents
incidence_23Dec2020_alt_analysis = incidence_percentile_percentages_output_alt_analysis{1};
incidence_24Dec2020_alt_analysis = incidence_percentile_percentages_output_alt_analysis{2};
incidence_25Dec2020_alt_analysis = incidence_percentile_percentages_output_alt_analysis{3};
incidence_26Dec2020_alt_analysis = incidence_percentile_percentages_output_alt_analysis{4};
incidence_27Dec2020_alt_analysis = incidence_percentile_percentages_output_alt_analysis{5};

age_split_incidence_23Dec2020_alt_analysis = age_split_incidence_percentile_percentages_output_alt_analysis{1};
age_split_incidence_24Dec2020_alt_analysis = age_split_incidence_percentile_percentages_output_alt_analysis{2};
age_split_incidence_25Dec2020_alt_analysis = age_split_incidence_percentile_percentages_output_alt_analysis{3};
age_split_incidence_26Dec2020_alt_analysis = age_split_incidence_percentile_percentages_output_alt_analysis{4};
age_split_incidence_27Dec2020_alt_analysis = age_split_incidence_percentile_percentages_output_alt_analysis{5};

cumulative_infections_time_horizon_alt_analysis = cumul_infections_percentile_percentages_output_alt_analysis{1};
cumulative_infections_23Dec2020_27Dec2020_alt_analysis = cumul_infections_percentile_percentages_output_alt_analysis{2};

age_split_cumulative_infections_time_horizon_alt_analysis = age_split_cumul_infections_percentile_output_alt_analysis{1};
age_split_cumul_infections_23Dec2020_27Dec2020_alt_analysis = age_split_cumul_infections_percentile_output_alt_analysis{2};

%% Main analysis, 30% adherence

% Specify if analysing main scenario runs (flag value = true) or the
% sensitivity analysis with longer time horizon (flag value = false)
main_analysis_flag = true;

% Specify if looking at main adherence of probability of 70% (flag value = true) 
% or lower adherence run with probability of 30% (flag value = false)
main_adherence_flag = false;

% Call function to produce plots and return summary statistics
[incidence_percentile_percentages_output_main_analysis_low_adherence,...
     age_split_incidence_percentile_percentages_output_main_analysis_low_adherence,...
     cumul_infections_percentile_percentages_output_main_analysis_low_adherence,...
     age_split_cumul_infections_percentile_output_main_analysis_low_adherence] = ...
                generate_summ_stats_and_plots_fn(main_analysis_flag,...
                                                main_adherence_flag,...
                                                plot_percentages_flag,...
                                                generate_figures_flag,...
                                                save_file_flag);

% Unpack cell contents
incidence_23Dec2020_main_analysis_low_adherence = incidence_percentile_percentages_output_main_analysis_low_adherence{1};
incidence_24Dec2020_main_analysis_low_adherence = incidence_percentile_percentages_output_main_analysis_low_adherence{2};
incidence_25Dec2020_main_analysis_low_adherence = incidence_percentile_percentages_output_main_analysis_low_adherence{3};
incidence_26Dec2020_main_analysis_low_adherence = incidence_percentile_percentages_output_main_analysis_low_adherence{4};
incidence_27Dec2020_main_analysis_low_adherence = incidence_percentile_percentages_output_main_analysis_low_adherence{5};

age_split_incidence_23Dec2020_main_analysis_low_adherence = age_split_incidence_percentile_percentages_output_main_analysis_low_adherence{1};
age_split_incidence_24Dec2020_main_analysis_low_adherence = age_split_incidence_percentile_percentages_output_main_analysis_low_adherence{2};
age_split_incidence_25Dec2020_main_analysis_low_adherence = age_split_incidence_percentile_percentages_output_main_analysis_low_adherence{3};
age_split_incidence_26Dec2020_main_analysis_low_adherence = age_split_incidence_percentile_percentages_output_main_analysis_low_adherence{4};
age_split_incidence_27Dec2020_main_analysis_low_adherence = age_split_incidence_percentile_percentages_output_main_analysis_low_adherence{5};

cumulative_infections_time_horizon_main_analysis_low_adherence = cumul_infections_percentile_percentages_output_main_analysis_low_adherence{1};
cumulative_infections_23Dec2020_27Dec2020_main_analysis_low_adherence = cumul_infections_percentile_percentages_output_main_analysis_low_adherence{2};

age_split_cumulative_infections_time_horizon_main_analysis_low_adherence = age_split_cumul_infections_percentile_output_main_analysis_low_adherence{1};
age_split_cumul_infections_23Dec2020_27Dec2020_main_analysis_low_adherence = age_split_cumul_infections_percentile_output_main_analysis_low_adherence{2};

%% Alternative analysis, 30% adherence

% Specify if analysing main scenario runs (flag value = true) or the
% sensitivity analysis with longer time horizon (flag value = false)
main_analysis_flag = false;

% Specify if looking at main adherence of probability of 70% (flag value = true) 
% or lower adherence run with probability of 30% (flag value = false)
main_adherence_flag = false;

% Call function to produce plots and return summary statistics
[incidence_percentile_percentages_output_alt_analysis_low_adherence,...
     age_split_incidence_percentile_percentages_output_alt_analysis_low_adherence,...
     cumul_infections_percentile_percentages_output_alt_analysis_low_adherence,...
     age_split_cumul_infections_percentile_output_alt_analysis_low_adherence] = ...
                generate_summ_stats_and_plots_fn(main_analysis_flag,...
                                                main_adherence_flag,...
                                                plot_percentages_flag,...
                                                generate_figures_flag,...
                                                save_file_flag);

% Unpack cell contents
incidence_23Dec2020_alt_analysis_low_adherence = incidence_percentile_percentages_output_alt_analysis_low_adherence{1};
incidence_24Dec2020_alt_analysis_low_adherence = incidence_percentile_percentages_output_alt_analysis_low_adherence{2};
incidence_25Dec2020_alt_analysis_low_adherence = incidence_percentile_percentages_output_alt_analysis_low_adherence{3};
incidence_26Dec2020_alt_analysis_low_adherence = incidence_percentile_percentages_output_alt_analysis_low_adherence{4};
incidence_27Dec2020_alt_analysis_low_adherence = incidence_percentile_percentages_output_alt_analysis_low_adherence{5};

age_split_incidence_23Dec2020_alt_analysis_low_adherence = age_split_incidence_percentile_percentages_output_alt_analysis_low_adherence{1};
age_split_incidence_24Dec2020_alt_analysis_low_adherence = age_split_incidence_percentile_percentages_output_alt_analysis_low_adherence{2};
age_split_incidence_25Dec2020_alt_analysis_low_adherence = age_split_incidence_percentile_percentages_output_alt_analysis_low_adherence{3};
age_split_incidence_26Dec2020_alt_analysis_low_adherence = age_split_incidence_percentile_percentages_output_alt_analysis_low_adherence{4};
age_split_incidence_27Dec2020_alt_analysis_low_adherence = age_split_incidence_percentile_percentages_output_alt_analysis_low_adherence{5};

cumulative_infections_time_horizon_alt_analysis_low_adherence = cumul_infections_percentile_percentages_output_alt_analysis_low_adherence{1};
cumulative_infections_23Dec2020_27Dec2020_alt_analysis_low_adherence = cumul_infections_percentile_percentages_output_alt_analysis_low_adherence{2};

age_split_cumulative_infections_time_horizon_alt_analysis_low_adherence = age_split_cumul_infections_percentile_output_alt_analysis_low_adherence{1};
age_split_cumul_infections_23Dec2020_27Dec2020_alt_analysis_low_adherence = age_split_cumul_infections_percentile_output_alt_analysis_low_adherence{2};



