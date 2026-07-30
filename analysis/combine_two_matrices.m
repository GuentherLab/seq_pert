function target_matrix = combine_two_matrices(add_matrix, target_matrix, target_i,ndims, just_add_NaNs)
    % collapse add_matrix if needed
    if ndims == 2
        sz = size(add_matrix);
        if sz(1) > 1
            add_matrix = mean(add_matrix,1,"omitnan");
        end

        if size(add_matrix,2) < size(target_matrix,2) % if the new row is smaller than the matrix
            temp_nans = NaN(1,[size(target_matrix,2)-size(add_matrix,2)]);
            add_matrix = cat(2, add_matrix,temp_nans);

            target_matrix(target_i,:) = add_matrix;
        elseif size(add_matrix,2) > size(target_matrix,2) % if the new row is larger than the matrix
            % add NaNs to the end of each of the existing rows
            start_col = size(target_matrix,2)+1;
            end_col = size(add_matrix,2);
            target_matrix(:,start_col:end_col) = NaN;

            target_matrix(target_i,:) = add_matrix;
        else % if the new row is the same size as the matrix
            target_matrix(target_i,:) = add_matrix;
        end
    elseif ndims == 3
        if size(add_matrix,2) < size(target_matrix,2) % if the new row is smaller than the matrix
            temp_nans = NaN(size(target_matrix,1),[size(target_matrix,2)-size(add_matrix,2)]);
            add_matrix = cat(2, add_matrix,temp_nans);

            if just_add_NaNs
                error('Function used incorrectly. Target matrix should be smaller than the add matrix for the "just add NaNs" functionality.');
            end

            target_matrix(:,:,target_i) = add_matrix;
        elseif size(add_matrix,2) > size(target_matrix,2) % if the new row is larger than the matrix
            % add NaNs to the end of each of the existing rows
            start_col = size(target_matrix,2)+1;
            end_col = size(add_matrix,2);
            target_matrix(:,start_col:end_col,:) = NaN;

            if just_add_NaNs
                return
            end
    
            target_matrix(:,:,target_i) = add_matrix;
            %novel_first_f1comp(:,:,i) = cat(2, novel_first_f1comp,temp);
        else % if the new row is the same size as the matrix
            target_matrix(:,:,target_i) = add_matrix;
        end
    else 
        error('unrecognized number of dims')
    end
end