function [ output_args ] = Select_Highest_SR_ZScore( input_args )

    mat = mat(all(~isnan(mat),2),:); % for nan - rows
    sorted_mat = sortrows(mat, 7);
    l = size(sorted_mat,1);
    sorted_mat = sorted_mat((l - 50):l,:);
    sorted_mat = sorted_mat(all(~isnan(sorted_mat),2),:); % for nan - rows
    sorted_mat = sorted_mat(:,all(~isnan(sorted_mat)));   % for nan - columns
    mean(sorted_mat, 1)
    
end

