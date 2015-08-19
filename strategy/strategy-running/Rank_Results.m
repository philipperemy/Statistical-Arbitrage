function [ sorted_mat ] = Rank_Results( mat, idx_col_sort, elt_count)

    mat = mat(all(~isnan(mat),2),:); % for nan - rows
    sorted_mat = sortrows(mat, idx_col_sort);
    l = size(sorted_mat,1);
    sorted_mat = sorted_mat((l - elt_count):l,:);
    sorted_mat = sorted_mat(all(~isnan(sorted_mat),2),:); % for nan - rows
    sorted_mat = sorted_mat(:,all(~isnan(sorted_mat)));   % for nan - columns

end

