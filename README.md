# TF_PAC
Computes PAC analysis over varying time windows, having been mainly used on resting-state EEG data
NOTE: EEGLAB must be installed and added to the folder you are working with, the version you are using will change things for the script, this one was built around eeglab2019_0.
NOTE: The MP plugin must also be installed. All of the encessary code for the plugin is found on this GitHub page: https://github.com/karolaug/mp-eeglab-plugin

Key Functions/Scripts: <br/>
MP_dynamicPAC.m <br/>
MP_dynamic_PAC_Animated.m <br/> 
extract_and_analyze_eeg_time_resolved.m <br/>

MP_dynamicPAC.m: This script uses matching pursuit to create clusters for the input dataset, which are then used to compute phase-amplitude coupling. The inputs for this script are all at the beginning, with necessary specifications being one channel of data over a certain time in a vector format (.mat file), and the sampling frequency of the data. This will display three commodulograms, one for low frequencies over time, one for high frequencies over time, and then one comparing which high and low frequencies coupled.

MP_dynamic_PAC_Animated.m: This script does the exact same thing, except the final result is an animated gif that is downloaded and can be reviewed, which shows how the high and low frequencies couple over the time frame by animating the heat map over the time window.

These scripts are limited though in time, and especially for resting state data were found to be much less effective for datasets encompassing larger than 3 minutes of data.

For longer windows, 
extract_and_analyze_eeg_time_resolved.m: This script specifies immediately what the required inputs are and it takes a single channel's data and shows how the high frequencies couple with specified low frequencies over the window length and how that compares to the actual EEG data waveform itself. The mp_syndata.mat file was included to as an example that can be tested with any of these scripts.
