clear;
load('dataset.mat')
%%% phase_series_MLE_MPPL         estimated phase series
%%% pol_detR                      pol_detR metric
%%% vv_cell                       VV SLC images stack of sentinel-1 satellite
%%% vh_cell                       VH SLC images stack of sentinel-1 satellite
%%% homo_num                      numbers of homogeneous pixels
%%% homo_num_thresh               threshhold for DS pixels identification
%%% homo_index                    homogeneous pixel index inside the
%%%                               homogeneous selection window 
%%% homo_wnd_rg                   windowsize for homogeneous pixels selection in
%%%                               range direction (must be odd)
%%% homo_wnd_az                   windowsize for homogeneous pixels selection in
%%%                               azimuth direction (must be odd)
%%% n_images                      number of images in the SLC stack
[phase_series_MLE_MPPL,pol_detR] = MLE_MPPL( vv_cell,vh_cell,homo_num, homo_num_thresh, homo_index,...
    homo_wnd_rg, homo_wnd_az, n_images );