function [match_grid,nonmatch_grid,approach_grid,dig_end_grid,leaveM_grid,leaveNM_grid] = local_EOIs_map(match_times, nonmatch_times, approach_times, ...
                                                                                                         dig_end_times, leaveM_times, leaveNM_times, ...
                                                                                                         n_datapoints, timestamps, order)

% EOI grids
match_grid =    zeros(n_datapoints, 1);
nonmatch_grid = zeros(n_datapoints, 1);
approach_grid = zeros(n_datapoints, 1);
dig_end_grid =  zeros(n_datapoints, 1);
leaveM_grid =   zeros(n_datapoints, 1);
leaveNM_grid =  zeros(n_datapoints, 1);

% Indices given provided EOI times
match_idx =    dsearchn(timestamps, match_times);
nonmatch_idx = dsearchn(timestamps, nonmatch_times);
approach_idx = dsearchn(timestamps, approach_times);
%dig_end_idx =  dsearchn(timestamps, dig_end_times);
%leaveM_idx =   dsearchn(timestamps, leaveM_times);
%leaveNM_idx =  dsearchn(timestamps, leaveNM_times);

% Populate the grids
for i = 0:(order-1)
    match_grid(match_idx+i,1) = 1;
    nonmatch_grid(nonmatch_idx+i,1) = 1;
    approach_grid(approach_idx+i,1) = 1;


end

return








