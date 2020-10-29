% this is an implementation of the tree-structured quantizer by a
% recurrence. this implementation is not really computationally efficient,
% but efficient in terms of code lines
function [U_back,ind_back,err_back,tree_searches,CB_searches_out] = Grass_quant_recurrence(U,d_i,bits_vec,percentage_vec,CB2)
%     dim1 = size(U,1);
    dim2 = size(U,2); 
    if d_i <= length(bits_vec)
        CB_size_high = max(ceil(2^bits_vec(d_i)),1);
        CB_size_low = max(floor(2^bits_vec(d_i)),1);
        prob = 2^bits_vec(d_i) - CB_size_low;
        bern = rand(1) <= prob;
        CB_size = CB_size_low*(1-bern) + CB_size_high*bern; % in case the codebook size is not an integer, we randomize between the integer below (floor) and above (ceil)
        CB = CB2{1}(:,:,1:CB_size); 
        quante = zeros(CB_size,1);
        for b_i = 1:CB_size
            temp = U'*CB(:,:,b_i);
            quante(b_i) = real(trace(temp*temp')); % quantization metric
        end
        [sort_inn,sort_ind] = sort(quante,'descend');
        quantization_width = max(round(CB_size*percentage_vec(d_i)),1); % width of the quantization tree at this level
        U_back = [];
        ind_back = [];
        err_back = [];
        tree_searches = 0;
        CB_searches_out = 0;
        for q_i = 1:quantization_width
             Uq = CB(:,:,sort_ind(q_i));
%              [W,L,~] = svd(U'*(Uq*Uq')*U);
%              B = Uq'*U*W*L^(-1/2); % this is a different form of the SQBC combiner, which spans the same subspace 
             B = Uq'*U*(U'*(Uq*Uq')*U)^(-1/2); % this is the SQBC combiner
             [U_in,ind_in,err_in,tree_searches_in,CB_searches_in] = Grass_quant_recurrence(B,d_i+1,bits_vec,percentage_vec,{CB2{2:end}}); % call yourself
             CB_searches_out = CB_searches_out + CB_searches_in;
             for qq_i = 1:size(U_in,3)
                Ut = Uq*U_in(:,:,qq_i);
                U_back = cat(3,U_back,Ut);   
%                 if d_i == 1
%                     err_back = [err_back;1-real(trace(U'*(Ut*Ut')*U))/dim2];
% %                     CB_searches_out = CB_searches_out+1;
%                 end
             end    
             ind_back = [ind_back;cat(2,repmat(sort_ind(q_i),max(size(ind_in,1),1),1),ind_in)];            
             err_back = [err_back;1-sort_inn(q_i)/dim2*(1-err_in)]; % quantization error recursively updated
             tree_searches = tree_searches + tree_searches_in; % this counts the number of branches of the tree visited
        end
        CB_searches_out = CB_searches_out + CB_size; % this counts the total number of codebook searches
    else % this part ends the recurrence
        U_back = 1;
        ind_back = [];    
        err_back = 0;
        tree_searches = 1;
        CB_searches_out = 0;
    end
end