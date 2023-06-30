// typedef struct {
//     int pos;
//     int8_t type = -1;
// } cut_pos;




// // for RRBS
// void parse_cut_rec(char *cut_str, std::vector<cut_rec>& cut_vec)
// {
//     cut_vec.clear();
//     int c, idx = 0;
//     cut_rec tmp_cut;
    
//     for(int i =0; cut_str[i] !='\0'; ++i){
//         if (cut_str[i] == '|'){tmp_cut.idx = idx; continue;}
//         if (cut_str[i] == ','){
//             if (tmp_cut.idx >=0) {
//                 tmp_cut.len = idx;
//                 cut_vec.push_back(tmp_cut);
//                 tmp_cut = {}; idx = 0; continue; // empty the struct, restart
//             } else {
//                 fprintf (stderr, "[%s] ERROR: Invalid enzyme site format (%s). Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);
//             }
//         }
//         ++idx;
//         c = nst_nt4_table[(int)cut_str[i]];
//         tmp_cut.seq.push_back(c);
//     }
//     if(tmp_cut.idx >=0){tmp_cut.len = idx; cut_vec.push_back(tmp_cut);}
// }

// void gen_cut_pos(const kseq_t *ks, std::vector<cut_pos>& cutpos_vec, std::vector<cut_rec>& cut_vec)
// {
//     cutpos_vec.clear(); // clean the container

//     std::vector<mut_t> rseq_ref(ks->seq.l); // TODO: create cut_pos without create rseq_ref

//     for (int i = 0; i != ks->seq.l; ++i) {
//         rseq_ref.push_back((mut_t)nst_nt4_table[(int)ks->seq.s[i]]);
//     }

//     std::vector<mut_t>::iterator ptr_begin = rseq_ref.begin();  // can repalce the type with `auto`
//     std::vector<mut_t>::iterator ptr_end   = rseq_ref.end();
//     std::vector<mut_t>::iterator iter_curr = rseq_ref.begin();  // current
//     std::vector<mut_t>::iterator iter_temp = rseq_ref.begin();  // temp
//     std::vector<mut_t>::iterator iter_save = rseq_ref.begin();  // save
//     //printf("%d %d\n", ptr_begin, ptr_end);

//     int count = 0; 
//     cut_pos tmp_cutpos;
//     while (iter_curr < ptr_end)
//     {
//         // printf("========%i========\n", count);
//         // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
//         //initiate
//         int min_pos = ptr_end - ptr_begin;
//         int site_pos= min_pos;
//         int type_idx = -1;

//         // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
//         for (int j = 0; j < cut_vec.size(); ++j) {
//             iter_temp = std::search(iter_curr, ptr_end, cut_vec[j].seq.begin(), cut_vec[j].seq.end());
//             // printf("%d %d %d\n", iter_curr, iter_temp);
//             site_pos = iter_temp - ptr_begin;
//             if(site_pos < min_pos){
//                 min_pos = site_pos;
//                 type_idx= j;
//                 iter_save= iter_temp;
//             }
//         }
//         // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
//         // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);

//         if(type_idx >=0){
//             tmp_cutpos = {};
//             tmp_cutpos.pos = min_pos + cut_vec[type_idx].idx;
//             tmp_cutpos.type= type_idx;
//             cutpos_vec.push_back(tmp_cutpos);
//             iter_curr = iter_save + cut_vec[type_idx].len; // need to check for other int types
//         }
//         // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
//         // ++count;
//         // if(count >= 20){
//         //     break;
//         // }
//     }
// }


// void collect_len_score_chr(const kseq_t *ks, chr_rec *tmp_len, char *bed_file, std::vector<probe_rec>& probe_vec)
// {
//     uint64_t eff_len;
//     float sum_score;
//     bool bool_bed_set = strcmp(bed_file,"None") && strlen(bed_file);

//     if(bool_bed_set){                       // targeted sequencing or
//         parse_bed_chr(bed_file, ks->name.s, probe_vec);
//         int len, pos_l, pos_r, pos_l_prev, pos_r_prev;
//         float score;
//         for (size_t i = 0; i < probe_vec.size(); ++i){
//             pos_l = probe_vec[i].pos_l;
//             pos_r = probe_vec[i].pos_r;
//             len   = pos_r - pos_r;
//             score = probe_vec[i].score;
//             // deal with overlap regions
//             if ((pos_l < pos_r_prev) && (pos_l > pos_l_prev)){ 
//                 len = (pos_r > pos_r_prev) ? (pos_r - pos_r_prev) : 0;
//             }
//             eff_len += len;
//             sum_score += score;
//             pos_l_prev = pos_l;
//             pos_r_prev = pos_r;
//         }
//     // }else if (tech_mode==1){                // reduced representative sequencing 
//     //     gen_cut_pos(ks, cutpos_vec, cut_vec);
//     //     int len;
//     //     for (int i = -1 ; i <= (int)cutpos_vec.size(); ++i){
//     //         if(i == -1 || i == (int)cutpos_vec.size()){
//     //             len = i==-1 ? cutpos_vec[0].pos : (int) ks->seq.l - cutpos_vec[i-1].pos;
//     //         }else{
//     //             len = cutpos_vec[i+1].pos - cutpos_vec[i].pos;
//     //         }
//     //         len = len >= min_insert && len <= max_insert? len : 0;
//     //         eff_len += len;
//     //     }
//     }else{                                  // whole genome
//         eff_len = ks->seq.l;
//         sum_score = eff_len;
//     }
//     tmp_len->chr_len = ks->seq.l;
//     tmp_len->eff_len = eff_len;
//     tmp_len->score   = sum_score;
// }

