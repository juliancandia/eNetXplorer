eNetXplorerExport <- function (x, path=getwd(), delimiter="tab", input.data=T, output.data.summary=T, output.data=T)
{
    if (delimiter=="tab") {
        sep="\t"
        ext="txt"
    } else if (delimiter=="csv") {
        sep=","
        ext="csv"
    } else {
        stop("Delimiter must be \'tab\' or \'csv\'\n")
    }
    
    setwd(path)
    
    if (input.data) {
        root_str = "inputData_"
        
        output = cbind(c("",x$instance),rbind(x$feature,x$predictor))
        write(t(output),ncol=ncol(output),file=paste0(root_str,"predictor.",ext),sep=sep)
        
        output = rbind(c("instance","response"),cbind(x$instance,as.character(x$response)))
        write(t(output),ncol=ncol(output),file=paste0(root_str,"response.",ext),sep=sep)
        
        args.desc = c("family","nlambda","nlambda.ext","seed",
        "scaled","n_fold","n_run","n_perm_null","QF_label","cor_method")
        args.value = list(x$family,x$nlambda,
        x$nlambda.ext,x$seed,x$scaled,x$n_fold,x$n_run,x$n_perm_null,x$QF_label,x$cor_method)
        args.value[sapply(args.value, is.null)] <- NA
        output = rbind(c("argument","value"),cbind(args.desc,unlist(args.value)))
        write(t(output),ncol=ncol(output),file=paste0(root_str,"arguments.",ext),sep=sep)
    }
    
    if (output.data.summary) {
        root_str = "outputDataSummary"
        
        summary = coef(summary(x))
        output = rbind(colnames(summary),summary)
        write(t(output),ncol=ncol(output),file=paste0(root_str,".",ext),sep=sep)
    }
    
    if (output.data) {
        root_str = "outputData_"
        n_alpha = length(x$alpha)
        
        nlambda = max(sapply(x$lambda_values,length)) # to pad missing lambda values if needed
        for (i_alpha in 1:n_alpha) {
            pad_length = nlambda - length(x$lambda_values[[i_alpha]])
            if (pad_length>0) {
                x$lambda_values[[i_alpha]] = c(x$lambda_values[[i_alpha]],rep(NA,pad_length))
                x$lambda_QF_est[[i_alpha]] = c(x$lambda_QF_est[[i_alpha]],rep(NA,pad_length))
            }
        }

        output = rbind(paste0("a",x$alpha),matrix(unlist(x$lambda_values),ncol=n_alpha))
        write(t(output),ncol=ncol(output),file=paste0(root_str,"lambda_values.",ext),sep=sep)
        
        output = rbind(paste0("a",x$alpha),matrix(unlist(x$lambda_QF_est),ncol=n_alpha))
        write(t(output),ncol=ncol(output),file=paste0(root_str,"lambda_QF_est.",ext),sep=sep)
        
        if (x$family=="gaussian") {
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
        write(t(output),ncol=ncol(output),file=paste0(root_str,"predicted_values.",ext),sep=sep)
        
        if (x$family%in%c("gaussian","binomial")) {
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_coef_wmean)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_coef_wmean.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_coef_wsd)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_coef_wsd.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_freq_mean)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_freq_mean.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$feature_freq_sd)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_freq_sd.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_coef_wmean)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_coef_wmean.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_coef_wsd)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_coef_wsd.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_freq_mean)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_freq_mean.",ext),sep=sep)
            
            output = rbind(c("feature",paste0("a",x$alpha)),cbind(x$feature,as.matrix(x$null_feature_freq_sd)))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_freq_sd.",ext),sep=sep)
            
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
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_coef_wmean.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_coef_wsd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_coef_wsd.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_freq_mean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_freq_mean.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$feature_freq_sd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"feature_freq_sd.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_coef_wmean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_coef_wmean.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_coef_wsd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_coef_wsd.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_freq_mean[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_freq_mean.",ext),sep=sep)
            
            tmp = NULL
            for (i_class in 1:n_class) {
                tmp = cbind(tmp,as.matrix(x$null_feature_freq_sd[[i_class]]))
            }
            output = rbind(header,cbind(x$feature,tmp))
            write(t(output),ncol=ncol(output),file=paste0(root_str,"null_feature_freq_sd.",ext),sep=sep)
            
        }
    }
}
