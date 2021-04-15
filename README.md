# motion-calcium-trace-data-extraction-analysis

image_projection_sudi is used to generate max_min field of view from the selected .tif files

SemiSeg_improve is used to select ROIs from the max_min figure

t_extract_trace_hdf5 is used to extract calcium traces from all the .tif files with the ROI_list

Clear_Velocity is used to extract speed trace from the raw .txt file

t_filter is a function within Clear_Velocity

var_low_high_speed is used to determine low and high speed threshold

motion_event is used to detect other locomotion features with Sigmoid-based method

New_motion_event is used to detect motion event with search box method

