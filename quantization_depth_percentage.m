%% Matlab Code for "Tree-Structured Quantization on Grassmann and Stiefel Manifolds", S. Schwarz et al., DCC 2021
% (c) Stefan Schwarz @ Institute of Telecommunications, TU Wien 2020

clc;
clear all;
close all;

n = 16;  % dimension of embedding space
m = 2;   % dimension of source sample
dim_vec = n:-1:m; % intermediate dimension of the quantizer
Dim = length(dim_vec)-1;
NN_in = 25; % inner randomization loop (constant quantization codebooks within this loop)
NN_out = 25; % outer randomization loop (random reinitialization of the codebooks)
bit_scale_vec = 1:9; % number of quantization bits simulated (multiplied by 14 below)
multi_stages = 4; % this determines how many stages/levels of the tree are split

% storage of results
quant_single_out_grass = zeros(NN_in,NN_out,length(bit_scale_vec));
quant_single_out_stief = zeros(NN_in,NN_out,length(bit_scale_vec));
quant_multi_out_grass = zeros(NN_in,NN_out,length(bit_scale_vec));
quant_multi_out_stief = zeros(NN_in,NN_out,length(bit_scale_vec));
grass_num_out = zeros(NN_in,NN_out,length(bit_scale_vec));
grass_search_out = zeros(NN_in,NN_out,length(bit_scale_vec));
stief_num_out = zeros(NN_in,NN_out,length(bit_scale_vec));
stief_search_out = zeros(NN_in,NN_out,length(bit_scale_vec));

for bit_scale_i = 1:length(bit_scale_vec) 
    
    bit_scale_i

bit_base = 14;
Nbs = bit_base*bit_scale_vec(bit_scale_i);  % number of quantization bits of current stage
CB_size_single = 2^(Nbs);
bits_vec = Nbs/Dim*ones(Dim,1);   % number of quantization bits per stage of the multi stage/tree-structured quantizers
% bits_vec(end-1:end) = bits_vec(end-1:end)+1;
percentage_vec = 2^(-bits_vec(1))*ones(Dim,1); % this is used for determining the width of the quantization tree
tree_search = true; % activate/deactivate the tree search
if tree_search    
    percentage_vec(1:multi_stages) = percentage_vec(1:multi_stages)*2; % number of levels of the tree that are split and corresponding width of each level (in terms of percentage of the number of codebook entries of each stage)
end
% percentage_vec(4) = 2^(-3);
% percentage_vec(end-8:end-1) = 2^(-3);
% percentage_vec(end-4:end-1) = percentage_vec(end-4:end-1)*2;
opt_bit = false; % activate this to optimize the bit allocation (only works for Grassmannian and only makes sense if tree_search = false)

% this is required for the calculation of the expected distortion of the individual stages of the Grassmannian multi stage quantizer 
dist_fac = zeros(1,Dim);
for nt = 1:Dim
    pp = m;
    nn = dim_vec(nt);
    qq = dim_vec(nt+1);
    c1 = 1/gamma(pp*(nn-qq)+1);
    for i = 1:pp
        c1 = c1*gamma(nn-i+1)/gamma(qq-i+1);
    end
    pre_factor1 = gamma(1/(pp*(nn-qq)))/(pp*(nn-qq));
    K1 = (c1)^(-1/(pp*(nn-qq)));
    dist_fac(nt) = pre_factor1*K1;
end

if opt_bit % optimized bit allocation according to "Reduced complexity recursive Grassmannian quantization", S. Schwarz et al., IEEE SPL 2020
    dim_diff = abs(diff(dim_vec));
    bits_dim = Nbs/Dim;
    cvx_begin quiet % optimized bit-allocation
    variable bit_alloc(1,Dim)
    %     maximize sum(log(1-dist_fac.*2.^(-bit_alloc./dim_diff)))
    maximize geo_mean(1-dist_fac/m.*2.^(-bit_alloc./(m*dim_diff)))
    subject to
    sum(bit_alloc) <= bits_dim*Dim
    bit_alloc >= 0
    bit_alloc <= 14
    cvx_end
    bits_vec = bit_alloc.';
    percentage_vec = 2.^(-bits_vec);
end

CB_size_vec = ceil(2.^(bits_vec));


% this is the theoretical distortion of single stage quantization
pp = m;
nn = n;
qq = m;
c1 = 1/gamma(pp*(nn-qq)+1);
for i = 1:pp
    c1 = c1*gamma(nn-i+1)/gamma(qq-i+1);
end
pre_factor1 = gamma(1/(pp*(nn-qq)))/(pp*(nn-qq));
K1 = (c1*2^(Nbs))^(-1/(pp*(nn-qq)));
dc_theor_single = pre_factor1*K1/m;

temp_bits = cumsum(bits_vec+log2(percentage_vec));
temp_bits = [0;temp_bits(1:end-1)];
eff_bits_vec = bits_vec + temp_bits;
CB_search_complexity = sum(2.^(eff_bits_vec)); % search complexity of the tree search in terms of the number of searched codebook entries

% Li = 2.^(bits_vec+log2(percentage_vec));
% Dk = 2.^(bits_vec);
% CB_search_complexity2 = 0;
% for kk = 1:Dim
%     L_prod = 1;
%     for ii = 1:kk-1
%         L_prod = L_prod*Li(ii);
%     end
%     CB_search_complexity2 = CB_search_complexity2 + L_prod*Dk(kk);
% end

% this is the theoretical distortion of the recursive multi stage quantizer
eff_bits_vec = bits_vec;
dc_theor = zeros(Dim,1);
for nt = 1:Dim
    pp = m;
    nn = dim_vec(nt);
    qq = dim_vec(nt+1);
    K1 = (2^(eff_bits_vec(nt)))^(-1/(pp*(nn-qq)));
    if eff_bits_vec(nt) <= 1 && pp == 1
        dc_temp1 = 1-qq/nn; % random isotropic projection in case of 0bit codebook
        dc_temp2 =  dist_fac(nt)*K1/m;
        prob = max((2^(eff_bits_vec(nt))-max(floor(2^(eff_bits_vec(nt))),1)),0);
        dc_theor(nt) = (1-prob)*dc_temp1 + prob*dc_temp2;
    else
        dc_theor(nt) = dist_fac(nt)*K1/m;    % theoretic normalized chordal distance (normalized by Nr)
    end
end
dc_theor_full = 1-prod(1-dc_theor);


r_stream = RandStream('mt19937ar','Seed',1);
r_stream_out = RandStream('mt19937ar','Seed',2);

% choose which manifold to simulate
activate_Grassmann = true;
activate_Stiefel = false;
% activate single stage quantization?
activate_single = false; % generally too complex

quant_single_store_grass = zeros(NN_in,NN_out);
quant_single_store_stief = zeros(NN_in,NN_out);
quant_multi_store_grass = zeros(NN_in,NN_out);
quant_multi_store_stief = zeros(NN_in,NN_out);
grass_num_store = zeros(NN_in,NN_out);
grass_search_store = zeros(NN_in,NN_out);
stief_num_store = zeros(NN_in,NN_out);
stief_search_store = zeros(NN_in,NN_out);
parfor nn_o = 1:NN_out     
    if ~mod(nn_o-1,10)
        nn_o
    end
    if activate_single
        CB = RANDOM_MIMO_CB(m,n,CB_size_single,r_stream_out,false,1); % single stage codebook
    end
    CB2 = cell(Dim,1);
    for d_i = 1:Dim
        CB2{d_i} = RANDOM_MIMO_CB(dim_vec(d_i+1),dim_vec(d_i),CB_size_vec(d_i),r_stream_out,false,1); % multi-stage codebooks
    end    
    
    quant_single_grass = zeros(NN_in,1);
    quant_single_stief = zeros(NN_in,1);
    quant_multi_grass = zeros(NN_in,1);
    quant_multi_stief = zeros(NN_in,1);
    grass_num = zeros(NN_in,1);
    grass_search = zeros(NN_in,1);
    stief_num = zeros(NN_in,1);
    stief_search = zeros(NN_in,1);
    for nn_i = 1:NN_in        
        % random point on Grassmann/Stiefel manifold
        [U0,~,~] = svd(randn(r_stream,n,m) + 1i*randn(r_stream,n,m),'econ');
        
        %% Grassmann manifold
        if activate_Grassmann
            % single-stage Grassmann manifold
            if activate_single
                quant_err = zeros(CB_size_single,1);
                for cb_i=1:CB_size_single
                    quant_err(cb_i) = 1 - real(trace(U0'*(CB(:,:,cb_i)*CB(:,:,cb_i)')*U0))/m;
                end
                quant_single_grass(nn_i) = min(quant_err);
            end
            
            % multi-stage Grassmann manifold
            d_ii = 1;
            B = U0;
            [U_back,ind_back,err_back,tree_searches,CB_searches] = Grass_quant_recurrence(B,d_ii,bits_vec,percentage_vec,CB2); % implementation via recurrence is rather slow; for-loop would be more efficient
            grass_num(nn_i) = tree_searches;
            grass_search(nn_i) = CB_searches; % notice, for a full-tree search it is actually inefficient to quantize within each stage
            % CSI reconstruction to test if it works properly
            [min_err,min_ind] = min(err_back);
            Ut = 1;
            for d_i = 1:Dim
                Ut = Ut*CB2{d_i}(:,:,ind_back(min_ind,d_i));
            end
            quant_multi_grass(nn_i) = 1 - real(trace(U0'*(Ut*Ut')*U0))/m;
        end
        %% Stiefel manifold
        if activate_Stiefel
            % single-stage Stiefel manifold
            if activate_single
                quant_err = zeros(CB_size_single,1);
                for cb_i=1:CB_size_single
                    quant_err(cb_i) = norm(U0-CB(:,:,cb_i),'fro')/(2*m);
                end
                quant_single_stief(nn_i) = min(quant_err);
            end
            
            % multi-stage Stiefel manifold
            d_ii = 1;
            B = U0;
            [U_back,ind_back,err_back,tree_searches,CB_searches] = Stief_quant_recurrence(B,d_ii,bits_vec,percentage_vec,CB2); % implementation via recurrence is rather slow; for-loop would be more efficient
            stief_num(nn_i) = tree_searches;
            stief_search(nn_i) = CB_searches;
            % CSI reconstruction to test if it works properly
            [min_err,min_ind] = min(err_back);
            Ut = 1;
            for d_i = 1:Dim
                Ut = Ut*CB2{d_i}(:,:,ind_back(min_ind,d_i));
            end
            quant_multi_stief(nn_i) = norm(U0-Ut,'fro')^2/(2*m);
        end
    end
    quant_single_store_grass(:,nn_o) = quant_single_grass;
    quant_single_store_stief(:,nn_o) = quant_single_stief;
    quant_multi_store_grass(:,nn_o) = quant_multi_grass;
    quant_multi_store_stief(:,nn_o) = quant_multi_stief;
    grass_num_store(:,nn_o) = grass_num;
    grass_search_store(:,nn_o) = grass_search;
    stief_num_store(:,nn_o) = stief_num;
    stief_search_store(:,nn_o) = stief_search;
end
if activate_Grassmann
    mean(quant_single_store_grass(:))
    mean(quant_multi_store_grass(:))
    2^Nbs
    mean(grass_search_store(:))
end
if activate_Stiefel
    mean(quant_single_store_stief(:))
    mean(quant_multi_store_stief(:))
    2^Nbs
    mean(stief_search_store(:))
end
    quant_single_out_grass(:,:,bit_scale_i) = quant_single_store_grass;
    quant_single_out_stief(:,:,bit_scale_i) = quant_single_store_stief;
    quant_multi_out_grass(:,:,bit_scale_i) = quant_multi_store_grass;
    quant_multi_out_stief(:,:,bit_scale_i) = quant_multi_store_stief;
    grass_num_out(:,:,bit_scale_i) = grass_num_store;
    grass_search_out(:,:,bit_scale_i) = grass_search_store;
    stief_num_out(:,:,bit_scale_i) = stief_num_store;
    stief_search_out(:,:,bit_scale_i) = stief_search_store;
    dc_theor_single_out(bit_scale_i) = dc_theor_single;
    percentage_vec_out(:,bit_scale_i) = percentage_vec;
    bits_vec_out(:,bit_scale_i) = bits_vec;
    CB_search_complexity_out(bit_scale_i) = CB_search_complexity;
end
save(['Grass_' num2str(n) 'x' num2str(m) '_' num2str(Nbs) 'bits_' num2str(tree_search) 'tree_' num2str(opt_bit) '_' num2str(multi_stages) 'multi.mat'],'dc_theor_single_out','quant_multi_out_grass','grass_search_out','percentage_vec_out','bits_vec_out','CB_search_complexity_out')
if activate_Grassmann
    figure(1)
    semilogy(bit_base*bit_scale_vec/(m*n),squeeze(mean(mean(quant_multi_out_grass,1),2)),'-o','linewidth',2,'MarkerSize',6);
    hold on
    grid on
    if activate_single
        semilogy(bit_base*bit_scale_vec/(m*n),squeeze(mean(mean(quant_single_out_grass,1),2)),'-o','linewidth',2,'MarkerSize',6);
    else
        semilogy(bit_base*bit_scale_vec/(m*n),dc_theor_single_out,'--o','linewidth',2,'MarkerSize',6);
    end
end
