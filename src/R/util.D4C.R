# this file submit jobs to data4cure to run analysis modules.
library(glue)

run_d4c_for_raw_count = function(cts.mat, sample.meta, compare_df, domain, user, pwd){
  col <- colnames(compare_df)
  row <- nrow(compare_df)
  for (i in 1:row) {
    tryCatch({
      # 1st part, extract control and test sample count and metadata
      ctrl = as.vector(compare_df[i,'control'])
      test = as.vector(compare_df[i,'test'])
      print(paste0('Process comparison ', ctrl, ' vs ', test))
      sub_cond_df = sample.meta[sample.meta$mergeCond %in% c(ctrl,test),]
      row.names(sub_cond_df) = sub_cond_df$Sample
      sub_expr_df = cts.mat[,row.names(sub_cond_df)]
      
      # 2nd part, create sub-folder and output the table and meta data into files
      comparison = str_replace_all(paste(test, ctrl, sep="_VS_"), '[ :,>]','_')
      path = file.path(d4c_out_dir, comparison)
      dir.create(path)
      sub_count_fn = glue('{path}/{comparison}.tsv')
      write.table(sub_expr_df, sub_count_fn, sep='\t',quote=F)
      sub_meta_fn = glue('{path}/{comparison}_meta.tsv')
      write.table(sub_cond_df, sub_meta_fn,sep='\t',quote=F, row.names=F)
      sub_compare_fn = glue('{path}/compare.tsv')
      write.table(compare_df[i,], sub_compare_fn, sep='\t',quote=F, row.names=F)
      
      # 3rd part, run the d4c api
      code = file.path(d4c_code_path,'run_data4cure_for_raw_RNASeq.py')
      cmd = glue("python3 '{code}' -i '{sub_count_fn}' -m '{sub_meta_fn}' -c '{sub_compare_fn}' -v mergeCond -d count --domain '{domain}' --user '{user}' --pw '{pwd}'")
      cmd = gsub('\\\\','/',cmd)
      system(cmd)
    },
    error = function(err) {
      print(paste("Error:", err))
    },
    finally = function(f) {
    })    
  }
}


run_d4c_for_deseq2_results = function(compare_df, domain, user, pwd){
  col <- colnames(compare_df)
  row <- nrow(compare_df)
  for (i in 1:row) {
    tryCatch({
      # 1st part, extract control and test sample metadata
      ctrl = as.vector(compare_df[i,'control'])
      test = as.vector(compare_df[i,'test'])
      print(paste0('Process comparison ', ctrl, ' vs ', test))
      sub_cond_df = sample.meta[sample.meta$mergeCond %in% c(ctrl,test),]
      row.names(sub_cond_df) = sub_cond_df$Sample
      
      # 2nd part, create sub-folder and copy the DE results
      comparison = str_replace_all(paste(test, ctrl, sep="_VS_"), '[ :,>]','_')
      path = glue('{d4c_out_dir}/{comparison}_DESeq2_result')
      dir.create(path)
      deseq2_fn = glue('{outdir}/{comparison}/{comparison}.result.csv')
      cp_deseq2_fn = glue('{path}/{comparison}.result.csv')
      file.copy(deseq2_fn,cp_deseq2_fn)
      
      # 3rd part, run the d4c api
      code = file.path(d4c_code_path,'run_data4cure_for_DESeq2_results.py')
      deseq2_tbl = glue('{path}/{comparison}.result.tbl')
      cmd = glue("python3 '{code}' -i '{cp_deseq2_fn}' -o '{deseq2_tbl}' --domain '{domain}' --user {user} --pw {pwd}")
      cmd = gsub('\\\\','/',cmd)
      system(cmd)
      
    },
    error = function(err) {
      print(paste("Error:", err))
    },
    finally = function(f) {
    })    
  }
}
