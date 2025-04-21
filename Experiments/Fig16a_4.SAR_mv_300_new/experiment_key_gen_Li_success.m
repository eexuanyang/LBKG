close all; 
clc;
clear all;

codeword_len = 255; % actual number of bits per BCH codeword (code length)
infor_len = 131;
key_len_final = 128;

data_set_arr = ["230421" "230428" "230602" "230617"];

dev_arr = 1:2;
mat_arr = 1;
mv_len = 200;
thre_num = 4;
ds_ratio = 16;
sample_len = 1;

for data_set = data_set_arr
    if data_set == "230421"
        loc_arr = 1:17;
    elseif data_set == "230428"
        loc_arr = 1:94;
    elseif data_set == "230602"
        loc_arr = 1:15;
    elseif data_set == "230617"
        loc_arr = 1:16;
    end
    for dev = dev_arr
        for mat_index = mat_arr
            key_arr = zeros(1200*length(loc_arr),sample_len);
            sample_idx = 1;
            for loc = loc_arr
                load(['File Path:' num2str(str2num(data_set)) '\Proceed_li\user' num2str(dev) '\data' num2str(loc) '\li_proceed_'  num2str(mat_index)  '.mat']);
                
                %% Power normalization
                DMRS_arr = DMRS_arr_li(1:sample_len,:);
                ampli_tem = sqrt(1./mean((abs(DMRS_arr)).^2));
                DMRS_arr = DMRS_arr * ampli_tem;
                DMRS_arr_smooth = smoothdata(abs(DMRS_arr.'),1,'movmean',mv_len);

                DMRS_arr_proc = DMRS_arr_smooth;
    
                for col_idx = 1:size(DMRS_arr_proc,2)
                    temp_arr_norm = zeros(size(DMRS_arr_proc,1),1);
                    seg_num = 100;
                    seg_len = round(length(temp_arr_norm)/seg_num);
                    for seg_idx = 1:seg_num
                        a = DMRS_arr_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
                        temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
                    end
                    temp_arr_norm_sort = sort(temp_arr_norm);
                    
                    thre_arr(1) = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
                    thre_arr(2) = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
                    thre_arr(3) = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);
                    
                    temp_arr = zeros(length(temp_arr_norm),1);
                    temp_arr(temp_arr_norm<(thre_arr(1))) = 1;
                    temp_arr((temp_arr_norm>(thre_arr(1))) & (temp_arr_norm<(thre_arr(2)))) = 2;
                    temp_arr((temp_arr_norm>(thre_arr(2))) & (temp_arr_norm<(thre_arr(3)))) = 3;
                    temp_arr(temp_arr_norm>(thre_arr(3))) = 4;
            
                    key_arr(1+(sample_idx - 1)*length(temp_arr_norm):sample_idx*length(temp_arr_norm),col_idx) = temp_arr;
                end
                sample_idx = sample_idx + 1;
            end

        end
        if dev == 1
            key_arr_a = key_arr;
        else
            key_arr_b = key_arr;
        end
    end
end

key_arr_a_ds = key_arr_a(1:ds_ratio:end, :);
key_arr_b_ds = key_arr_b(1:ds_ratio:end, :);

key_arr_a_bit = zeros(1, 2*numel(key_arr_a_ds));
real_idx = 1;
for ele_idx = 1:numel(key_arr_a_ds)
    switch key_arr_a_ds(ele_idx)
        case 0
            disp('???')
            key_arr_a_bit(1+(real_idx-1)*2) = randi([0 1],1,1);
            key_arr_a_bit(2+(real_idx-1)*2) = randi([0 1],1,1);
            real_idx = real_idx + 1;
        case 1
            key_arr_a_bit(1+(real_idx-1)*2) = 0;
            key_arr_a_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
        case 2
            key_arr_a_bit(1+(real_idx-1)*2) = 0;
            key_arr_a_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 3
            key_arr_a_bit(1+(real_idx-1)*2) = 1;
            key_arr_a_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 4
            key_arr_a_bit(1+(real_idx-1)*2) = 1;
            key_arr_a_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
    end
end

key_arr_b_bit = zeros(1, 2*numel(key_arr_b_ds));
real_idx = 1;
for ele_idx = 1:numel(key_arr_b_ds)
    switch key_arr_b_ds(ele_idx)
        case 0
            disp('???')
            key_arr_b_bit(1+(real_idx-1)*2) = randi([0 1],1,1);
            key_arr_b_bit(2+(real_idx-1)*2) = randi([0 1],1,1);
            real_idx = real_idx + 1;
        case 1
            key_arr_b_bit(1+(real_idx-1)*2) = 0;
            key_arr_b_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
        case 2
            key_arr_b_bit(1+(real_idx-1)*2) = 0;
            key_arr_b_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 3
            key_arr_b_bit(1+(real_idx-1)*2) = 1;
            key_arr_b_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 4
            key_arr_b_bit(1+(real_idx-1)*2) = 1;
            key_arr_b_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
    end
end

key_arr_a_bit_vec = (key_arr_a_bit(:)).';
key_arr_b_bit_vec = (key_arr_b_bit(:)).';

key_per_len = 1200*2/ds_ratio;
key_num = length(key_arr_a_bit_vec)/key_per_len;

key_arr_a_bit_vec_discard = zeros(1,infor_len*key_num);
key_arr_b_bit_vec_discard = zeros(1,infor_len*key_num);

for key_idx = 1:key_num
    key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:infor_len*key_idx) = key_arr_a_bit_vec((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
    key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:infor_len*key_idx) = key_arr_b_bit_vec((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
end

%     key_arr_b_bit_vec_discard_ori = key_arr_b_bit_vec_discard;

for key_idx = 1:key_num
    key_temp = key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:key_idx*infor_len);
    msgTx = gf(key_temp);
    key_coded = bchenc(msgTx, codeword_len, infor_len);
    key_coded(1:infor_len) = logical(key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:key_idx*infor_len))';
    [key_corr,~] = bchdec(key_coded, codeword_len, infor_len); 
    key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:key_idx*infor_len) = key_corr.x;
end

key_arr_fdd_a_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
key_arr_fdd_b_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);

for key_idx = 1:key_num
    key_arr_fdd_a_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
    key_arr_fdd_b_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
end


key_sar_norm_sort_kdr = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_b_norm_sort_ds_bit_disc_final)))/length(key_arr_fdd_a_norm_sort_ds_bit_disc_final);
key_sar_norm_sort_kgr = 2*length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_b_norm_sort_ds_bit_disc_final)))/length(key_arr_fdd_a_norm_sort_ds_bit_disc_final);

save('key_sar_norm_sort_kdr.mat','key_sar_norm_sort_kdr')
save('key_sar_norm_sort_kgr.mat','key_sar_norm_sort_kgr')