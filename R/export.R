export <- function (x, dest_dir=getwd(), dest_dir_create=TRUE, delim=c("tsv","csv"), input.data=TRUE, summary.data=TRUE, output.data=TRUE)
{
    
    if ((!dir.exists(dest_dir))&&(dest_dir_create)) {
        dir.create(dest_dir)
    }
    
    delim = match.arg(delim)
    if (delim=="tsv") {
        sep="\t"
        ext="tsv"
    } else if (delim=="csv") {
        sep=","
        ext="csv"
    }
    
    if (input.data) {
        root_str = "inputData_"
        
        output = cbind(c("",x$instance),rbind(x$feature,as.matrix(x$predictor)))
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,"predictor.",ext)),sep=sep)
        if (x$family=="cox") {
            output = rbind(c("instance",colnames(x$response)),cbind(x$instance,x$response))
        } else {
            output = rbind(c("instance","response"),cbind(x$instance,as.character(x$response)))
        }
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,"response.",ext)),sep=sep)
        
        if (x$family=="gaussian") {
            args.desc = c("family","nlambda","nlambda.ext","seed","scaled","n_fold","n_run","n_perm_null","QF_label",
            "cor_method")
            args.value = list(x$family,x$nlambda,x$nlambda.ext,x$seed,x$scaled,x$n_fold,x$n_run,x$n_perm_null,x$QF_label,
            x$cor_method)
        } else if (x$family=="binomial") {
            args.desc = c("family","nlambda","nlambda.ext","seed","scaled","n_fold","n_run","n_perm_null","QF_label",
            "binom_method","binom_pos","fscore_beta","fold_distrib_fail.max")
            args.value = list(x$family,x$nlambda,x$nlambda.ext,x$seed,x$scaled,x$n_fold,x$n_run,x$n_perm_null,x$QF_label,
            x$binom_method,x$binom_pos,x$fscore_beta,x$fold_distrib_fail.max)
        } else if (x$family=="multinomial") {
            args.desc = c("family","nlambda","nlambda.ext","seed","scaled","n_fold","n_run","n_perm_null","QF_label",
            "multinom_method","fscore_beta","fold_distrib_fail.max")
            args.value = list(x$family,x$nlambda,x$nlambda.ext,x$seed,x$scaled,x$n_fold,x$n_run,x$n_perm_null,x$QF_label,
            x$multinom_method,x$fscore_beta,x$fold_distrib_fail.max)
        } else if (x$family=="cox") {
            args.desc = c("family","nlambda","nlambda.ext","seed","scaled","n_fold","n_run","n_perm_null","QF_label",
            "cox_index","logrank","survAUC","survAUC_time")
            args.value = list(x$family,x$nlambda,x$nlambda.ext,x$seed,x$scaled,x$n_fold,x$n_run,x$n_perm_null,x$QF_label,
            x$cox_index,x$logrank,x$survAUC,paste0(x$survAUC_time,collapse=","))
        }
        
        args.value[sapply(args.value, is.null)] <- NA
        output = rbind(c("argument","value"),cbind(args.desc,unlist(args.value)))
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,"arguments.",ext)),sep=sep)
    }
    
    if (summary.data) {
        root_str = "summaryData"
        
        summary = coef(summary(x))
        output = rbind(colnames(summary),summary)
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,".",ext)),sep=sep)
    }
    
    if (output.data) {
        root_str = "outputData_"
        n_alpha = length(x$alpha)
        
        if (x$save_lambda_QF_full) {
            for (i_alpha in 1:n_alpha) {
               output =  rbind(c("lambda",paste0("run_",1:x$n_run)),cbind(x$lambda_values[[i_alpha]],x$lambda_QF_est_full[[i_alpha]]))
               write(t(output),ncolumns=ncol(output),
               file=file.path(dest_dir,paste0(root_str,"lambda_QF_est_full_a",x$alpha[i_alpha],".",ext)),sep=sep)
            }
        }
        
        nlambda = max(sapply(x$lambda_values,length)) # to pad missing lambda values if needed
        for (i_alpha in 1:n_alpha) {
            pad_length = nlambda - length(x$lambda_values[[i_alpha]])
            if (pad_length>0) {
                x$lambda_values[[i_alpha]] = c(x$lambda_values[[i_alpha]],rep(NA,pad_length))
                x$lambda_QF_est[[i_alpha]] = c(x$lambda_QF_est[[i_alpha]],rep(NA,pad_length))
            }
        }

        output = rbind(paste0("a",x$alpha),matrix(unlist(x$lambda_values),ncol=n_alpha))
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,"lambda_values.",ext)),sep=sep)
        
        output = rbind(paste0("a",x$alpha),matrix(unlist(x$lambda_QF_est),ncol=n_alpha))
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,"lambda_QF_est.",ext)),sep=sep)
        
        if (x$family%in%c("gaussian","cox")) {
            pred_header = c("median","mad")
        } else if (x$family%in%c("binomial","multinomial")) {
            pred_header = levels(x$response)
        }
        header = "instance"
        for (i_alpha in 1:n_alpha) {
            header = c(header,paste0("a",x$alpha[i_alpha],"_",pred_header))
        }
        output = rbind(header,cbind(x$instance,matrix(unlist(x$predicted_values),
        ncol=n_alpha*length(pred_header))))
        write(t(output),ncolumns=ncol(output),
        file=file.path(dest_dir,paste0(root_str,"predicted_values.",ext)),sep=sep)
        
        if (x$family%in%c("gaussian","binomial","cox")) {
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_coef_wmean)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_coef_wmean.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_coef_wsd)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_coef_wsd.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_freq_mean)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_freq_mean.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_freq_sd)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_freq_sd.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_coef_wmean)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"null_feature_coef_wmean.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_coef_wsd)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"null_feature_coef_wsd.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_freq_mean)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"null_feature_freq_mean.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_freq_sd)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"null_feature_freq_sd.",ext)),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,
            as.matrix(x$feature_coef_model_vs_null_pval)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_coef_model_vs_null_pval.",ext))
            ,sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,
            as.matrix(x$feature_freq_model_vs_null_pval)))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_freq_model_vs_null_pval.",ext))
            ,sep=sep)
            
        } else if (x$family=="multinomial") {
            
            class = levels(x$response)
            n_class = length(class)
            header = "feature"
            for (i_class in 1:n_class) {
                header = c(header,paste0(class[i_class],"_a",x$alpha))
            }
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_coef_wmean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_coef_wmean.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_coef_wsd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"feature_coef_wsd.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_freq_mean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"feature_freq_mean.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_freq_sd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"feature_freq_sd.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_coef_wmean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"null_feature_coef_wmean.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_coef_wsd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"null_feature_coef_wsd.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_freq_mean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"null_feature_freq_mean.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_freq_sd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),file=file.path(dest_dir,paste0(root_str,"null_feature_freq_sd.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_coef_model_vs_null_pval[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_coef_model_vs_null_pval.",ext)),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_freq_model_vs_null_pval[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncolumns=ncol(output),
            file=file.path(dest_dir,paste0(root_str,"feature_freq_model_vs_null_pval.",ext)),sep=sep)
            
        }
        
        if (x$family=="cox") {
            if (x$logrank) {
                output = rbind(c("alpha","logrank_pval"),cbind(x$alpha,x$logrank_pval))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"logrank_pval.",ext)),sep=sep)
            }
            if (x$survAUC) {
                output = rbind(c("time",paste0("a",x$alpha)),cbind(x$survAUC_time,x$AUC_mean))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"survAUC_mean.",ext)),sep=sep)
                
                output = rbind(c("time",paste0("a",x$alpha)),cbind(x$survAUC_time,x$AUC_sd))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"survAUC_sd.",ext)),sep=sep)
                
                output = rbind(c("time",paste0("a",x$alpha)),cbind(x$survAUC_time,x$AUC_perc025))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"survAUC_perc025.",ext)),sep=sep)
                
                output = rbind(c("time",paste0("a",x$alpha)),cbind(x$survAUC_time,x$AUC_perc500))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"survAUC_perc500.",ext)),sep=sep)
                
                output = rbind(c("time",paste0("a",x$alpha)),cbind(x$survAUC_time,x$AUC_perc975))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"survAUC_perc975.",ext)),sep=sep)
                
                output = rbind(c("time",paste0("a",x$alpha)),cbind(x$survAUC_time,x$AUC_pval))
                write(t(output),ncolumns=ncol(output),
                file=file.path(dest_dir,paste0(root_str,"survAUC_pval.",ext)),sep=sep)
            }
        }
    }
}
