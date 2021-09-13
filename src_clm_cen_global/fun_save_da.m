function fun_save_da(data_path, sd_controling_factor, iworker, soc_mod_trace, litter_mod_trace, cwd_mod_trace, hr_mod_trace, parameters_keep2)

% da_summary.soc_mod_trace = soc_mod_trace;
% da_summary.litter_mod_trace = litter_mod_trace;
% da_summary.cwd_mod_trace = cwd_mod_trace;
% da_summary.hr_mod_trace = hr_mod_trace;
% da_summary.parameters_keep2 = parameters_keep2;

% save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(grid_num), '_reconstruct_worker_', num2str(iworker), '_soc_mod_trace.mat'], 'soc_mod_trace', '-v7.3');
% 
% save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(grid_num), '_reconstruct_worker_', num2str(iworker), '_litter_mod_trace.mat'], 'litter_mod_trace', '-v7.3');
% 
% save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(grid_num), '_reconstruct_worker_', num2str(iworker), '_cwd_mod_trace.mat'], 'cwd_mod_trace', '-v7.3');
% 
% save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(grid_num), '_reconstruct_worker_', num2str(iworker), '_hr_mod_trace.mat'], 'hr_mod_trace', '-v7.3');

save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_global_reconstruct_test3_s3_matrix_', num2str(sd_controling_factor), '_worker_', num2str(iworker), '_parameters_keep2.mat'], 'parameters_keep2', '-v7.3');



end
