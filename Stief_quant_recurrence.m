function [U_back,ind_back,err_back,tree_searches,CB_searches_out] = Stief_quant_recurrence(U,d_i,bits_vec,percentage_vec,CB2)
%     dim1 = size(U,1);
    dim2 = size(U,2); 
    
    CB_size_high = max(ceil(2^bits_vec(d_i)),1);
    CB_size_low = max(floor(2^bits_vec(d_i)),1);
    prob = 2^bits_vec(d_i) - CB_size_low;
    bern = rand(1) <= prob;
    CB_size = CB_size_low*(1-bern) + CB_size_high*bern; % in case the optimal codebook size is not an integer, we randomize between the integer below and above
    CB = CB2{1}(:,:,1:CB_size);
        
    if d_i < length(bits_vec)        
        quante = zeros(CB_size,1);
        for b_i = 1:CB_size
            temp = U'*CB(:,:,b_i);
            quante(b_i) = sum(svd(temp,'econ'));
        end
        [sort_inn,sort_ind] = sort(quante,'descend');
        quantization_depth = max(round(CB_size*percentage_vec(d_i)),1);
        U_back = [];
        ind_back = [];
        err_back = [];
        tree_searches = 0;
        CB_searches_out = 0;
        for q_i = 1:quantization_depth
             Uq = CB(:,:,sort_ind(q_i));
             [Y,L,V] = svd(Uq'*U,'econ');
             B = Y*V';
             [U_in,ind_in,err_in,tree_searches_in,CB_searches_in] = Stief_quant_recurrence(B,d_i+1,bits_vec,percentage_vec,{CB2{2:end}});
             CB_searches_out = CB_searches_out + CB_searches_in;
             for qq_i = 1:size(U_in,3)
                Ut = Uq*U_in(:,:,qq_i);
                U_back = cat(3,U_back,Ut);   
%                 if d_i == 1
%                     err_back = [err_back;norm(U-Ut,'fro')^2/(2*dim2)];
%                     CB_searches_out = CB_searches_out + 1;
%                 end
             end   
             ind_back = [ind_back;cat(2,repmat(sort_ind(q_i),max(size(ind_in,1),1),1),ind_in)];            
             err_back = [err_back;1-sort_inn(q_i)/(dim2)*(1-err_in)];
             tree_searches = tree_searches + tree_searches_in;
        end
        CB_searches_out = CB_searches_out + CB_size;
    else
        quante = zeros(CB_size,1);
        for b_i = 1:CB_size
            quante(b_i) = norm(U - CB(:,:,b_i),'fro')^2;
        end
        [sort_err,sort_ind] = sort(quante,'ascend');
        U_back = [];
        ind_back = [];        
        quantization_depth = max(round(CB_size*percentage_vec(d_i)),1);
        err_back = sort_err(1:quantization_depth)/(2*dim2);
        for q_i = 1:quantization_depth
            U_back = cat(3,U_back,CB(:,:,sort_ind(q_i))); 
            ind_back = [ind_back;sort_ind(q_i)];
        end
        tree_searches = quantization_depth;
        CB_searches_out = CB_size;
    end
end