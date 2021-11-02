# Cox model
eNetXplorerCox <- function(x, y, family, alpha, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, save_lambda_QF_full, QF.FUN, QF_label, cox_index, logrank, survAUC, survAUC_time, ...)
{
    
    n_instance = nrow(x)
    n_feature = ncol(x)
    instance = rownames(x)
    if (is.null(instance)) {
        instance = paste0("Inst.",1:nrow(x))
    }
    feature = colnames(x)
    if (is.null(feature)) {
        feature = paste0("Feat.",1:ncol(x))
    }
    prediction_type = "response"
    
    if (is.null(QF.FUN)) {
        if (cox_index=="concordance") {
            QF = function(predicted,response) {
                return (concordance.index(x=predicted,surv.time=response[,"time"],surv.event=response[,"status"])$c.index)
            }
            QF_label = "concordance index"
        } else if (cox_index=="D-index") {
            QF = function(predicted,response) {
                return (D.index(x=predicted,surv.time=response[,"time"],surv.event=response[,"status"])$d.index)
            }
            QF_label = "D-index"
        }
    } else {
        QF = QF.FUN
        if (is.null(QF_label)) {
            QF_label = "user-provided"
        }
    }
    
    if (!is.null(seed)) {
        set.seed(seed) # set the seed for reproducibility
    }
    
    if (scaled) { # data are z-score transformed
        x_col_mean = apply(x,2,mean)
        x_col_sd = apply(x,2,sd)
        for (i in 1:ncol(x)) {
            x[,i] = (x[,i]-x_col_mean[i])/x_col_sd[i]
        }
    }
    
    foldid = NULL # we determine the fold template
    fold_size = floor(n_instance/n_fold)
    for (i_fold in 1:n_fold) {
        foldid = c(foldid,rep(i_fold,fold_size))
    }
    fold_rest = n_instance%%n_fold
    if (fold_rest>0) {
        for (i_fold in 1:fold_rest) {
            foldid = c(foldid,i_fold)
        }
    }
    
    n_alpha = length(alpha)
    n_alpha_eff = 0
    best_lambda = rep(NA,n_alpha)
    model_QF_est = rep(NA,n_alpha)
    
    lambda_values = vector("list",n_alpha)
    lambda_QF_est = vector("list",n_alpha)
    lambda_QF_est_full = vector("list",n_alpha)
    
    foldid_per_run = vector("list",n_alpha)
    predicted_values = vector("list",n_alpha)
    
    feature_coef_wmean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_coef_wsd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_mean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_sd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    null_feature_coef_wmean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_coef_wsd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_freq_mean = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    null_feature_freq_sd = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    QF_model_vs_null_pval = rep(NA,n_alpha)
    feature_coef_model_vs_null_pval = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    feature_freq_model_vs_null_pval = matrix(rep(NA,n_feature*n_alpha),ncol=n_alpha)
    
    # logrank
    if (logrank) {
        risk_group_template = rep(1,n_instance)
        logrank_pval = rep(NA,n_alpha)
    }
    
    # AUC
    if (survAUC) {
        if (!"package:survival"%in%search()) {
            attachNamespace("survival") # to solve some issues with timeROC
        }
        if (is.null(survAUC_time)) {
            stop("Error: survAUC_time must be provided")
        }
        n_survAUC_time = length(survAUC_time)
        AUC_mean = matrix(rep(NA,n_survAUC_time*n_alpha),ncol=n_alpha)
        AUC_sd = matrix(rep(NA,n_survAUC_time*n_alpha),ncol=n_alpha)
        AUC_perc025 = matrix(rep(NA,n_survAUC_time*n_alpha),ncol=n_alpha)
        AUC_perc500 = matrix(rep(NA,n_survAUC_time*n_alpha),ncol=n_alpha)
        AUC_perc975 = matrix(rep(NA,n_survAUC_time*n_alpha),ncol=n_alpha)
        AUC_pval = matrix(rep(NA,n_survAUC_time*n_alpha),ncol=n_alpha)
    }
    
    glmnet.control(factory=T) # resets internal glmnet parameters
    glmnet.control(...) # to allow changes of factory default parameters in glmnet
    
    # return object
    return_obj <- function(n_alpha_eff) {
        rownames(x) = instance
        colnames(x) = feature
        rownames(y) = instance
        alpha_label = paste0("a",alpha[1:n_alpha_eff])
        best_lambda = best_lambda[1:n_alpha_eff]
        names(best_lambda) = alpha_label
        model_QF_est = model_QF_est[1:n_alpha_eff]
        names(model_QF_est) = alpha_label
        QF_model_vs_null_pval = QF_model_vs_null_pval[1:n_alpha_eff]
        names(QF_model_vs_null_pval) = alpha_label
        lambda_values = lambda_values[1:n_alpha_eff]
        names(lambda_values) = alpha_label
        lambda_QF_est = lambda_QF_est[1:n_alpha_eff]
        names(lambda_QF_est) = alpha_label
        lambda_QF_est_full = lambda_QF_est_full[1:n_alpha_eff]
        names(lambda_QF_est_full) = alpha_label
        feature_coef_wmean = feature_coef_wmean[,1:n_alpha_eff,drop=F]
        rownames(feature_coef_wmean) = feature
        colnames(feature_coef_wmean) = alpha_label
        feature_coef_wsd = feature_coef_wsd[,1:n_alpha_eff,drop=F]
        rownames(feature_coef_wsd) = feature
        colnames(feature_coef_wsd) = alpha_label
        feature_freq_mean = feature_freq_mean[,1:n_alpha_eff,drop=F]
        rownames(feature_freq_mean) = feature
        colnames(feature_freq_mean) = alpha_label
        feature_freq_sd = feature_freq_sd[,1:n_alpha_eff,drop=F]
        rownames(feature_freq_sd) = feature
        colnames(feature_freq_sd) = alpha_label
        null_feature_coef_wmean = null_feature_coef_wmean[,1:n_alpha_eff,drop=F]
        rownames(null_feature_coef_wmean) = feature
        colnames(null_feature_coef_wmean) = alpha_label
        null_feature_coef_wsd = null_feature_coef_wsd[,1:n_alpha_eff,drop=F]
        rownames(null_feature_coef_wsd) = feature
        colnames(null_feature_coef_wsd) = alpha_label
        null_feature_freq_mean = null_feature_freq_mean[,1:n_alpha_eff,drop=F]
        rownames(null_feature_freq_mean) = feature
        colnames(null_feature_freq_mean) = alpha_label
        null_feature_freq_sd = null_feature_freq_sd[,1:n_alpha_eff,drop=F]
        rownames(null_feature_freq_sd) = feature
        colnames(null_feature_freq_sd) = alpha_label
        feature_coef_model_vs_null_pval = feature_coef_model_vs_null_pval[,1:n_alpha_eff,drop=F]
        rownames(feature_coef_model_vs_null_pval) = feature
        colnames(feature_coef_model_vs_null_pval) = alpha_label
        feature_freq_model_vs_null_pval = feature_freq_model_vs_null_pval[,1:n_alpha_eff,drop=F]
        rownames(feature_freq_model_vs_null_pval) = feature
        colnames(feature_freq_model_vs_null_pval) = alpha_label
        predicted_values = predicted_values[1:n_alpha_eff]
        names(predicted_values) = alpha_label
        for (i_alpha in 1:n_alpha_eff) {
            colnames(predicted_values[[i_alpha]]) = c("median","mad")
            rownames(predicted_values[[i_alpha]]) = instance
        }
        if (logrank) {
            logrank_pval = logrank_pval[1:n_alpha_eff]
            names(logrank_pval) = alpha_label
        }
        if (survAUC) {
            AUC_mean = AUC_mean[,1:n_alpha_eff,drop=F]
            rownames(AUC_mean) = survAUC_time
            colnames(AUC_mean) = alpha_label
            AUC_sd = AUC_sd[,1:n_alpha_eff,drop=F]
            rownames(AUC_sd) = survAUC_time
            colnames(AUC_sd) = alpha_label
            AUC_perc025 = AUC_perc025[,1:n_alpha_eff,drop=F]
            rownames(AUC_perc025) = survAUC_time
            colnames(AUC_perc025) = alpha_label
            AUC_perc500 = AUC_perc500[,1:n_alpha_eff,drop=F]
            rownames(AUC_perc500) = survAUC_time
            colnames(AUC_perc500) = alpha_label
            AUC_perc975 = AUC_perc975[,1:n_alpha_eff,drop=F]
            rownames(AUC_perc975) = survAUC_time
            colnames(AUC_perc975) = alpha_label
            AUC_pval = AUC_pval[,1:n_alpha_eff,drop=F]
            rownames(AUC_pval) = survAUC_time
            colnames(AUC_pval) = alpha_label
        }
        
        output =
        list(
        # input data and parameters
        predictor = as(x,"CsparseMatrix"), response = y, alpha = alpha[1:n_alpha_eff], family = family, nlambda = nlambda,
        nlambda.ext = nlambda.ext, seed = seed, scaled = scaled, n_fold = n_fold, n_run = n_run,
        n_perm_null = n_perm_null, save_lambda_QF_full = save_lambda_QF_full,
        QF_label = QF_label, cox_index = cox_index, logrank = logrank, survAUC = survAUC,
        survAUC_time = survAUC_time,
        instance = instance, feature = feature, glmnet_params = glmnet.control(),
        # summary results
        best_lambda = best_lambda, model_QF_est = model_QF_est, QF_model_vs_null_pval = QF_model_vs_null_pval,
        # detailed results for plots and downstream analysis
        lambda_values = lambda_values, lambda_QF_est = lambda_QF_est, lambda_QF_est_full = lambda_QF_est_full,
        predicted_values = predicted_values,
        feature_coef_wmean = as(feature_coef_wmean,"CsparseMatrix"), feature_coef_wsd = as(feature_coef_wsd,"CsparseMatrix"),
        feature_freq_mean = as(feature_freq_mean,"CsparseMatrix"), feature_freq_sd = as(feature_freq_sd,"CsparseMatrix"),
        null_feature_coef_wmean = as(null_feature_coef_wmean,"CsparseMatrix"), null_feature_coef_wsd = as(null_feature_coef_wsd,"CsparseMatrix"),
        null_feature_freq_mean = as(null_feature_freq_mean,"CsparseMatrix"), null_feature_freq_sd = as(null_feature_freq_sd,"CsparseMatrix"),
        feature_coef_model_vs_null_pval = as(feature_coef_model_vs_null_pval,"CsparseMatrix"),
        feature_freq_model_vs_null_pval = as(feature_freq_model_vs_null_pval,"CsparseMatrix")
        )
        if (logrank) {
            output = append(output, list(logrank_pval = logrank_pval))
        }
        if (survAUC) {
            output = append(output, list(
            AUC_mean = AUC_mean, AUC_sd = AUC_sd, AUC_perc025 = AUC_perc025, AUC_perc500 = AUC_perc500, AUC_perc975 = AUC_perc975, AUC_pval = AUC_pval
            ))
        }
        output
    }
    
    tryCatch({
        for (i_alpha in 1:n_alpha) { # beginning of alpha loop
            fit = glmnet(x,y,alpha=alpha[i_alpha],family=family,nlambda=nlambda)
            # We extend the range of lambda values (if nlambda.ext is set)
            if (!is.null(nlambda.ext)&&nlambda.ext>nlambda) {
                lambda_max = fit$lambda[1]*sqrt(nlambda.ext/length(fit$lambda))
                lambda_min = fit$lambda[length(fit$lambda)]/sqrt(nlambda.ext/length(fit$lambda))
                lambda_values[[i_alpha]] = lambda_max*(lambda_min/lambda_max)**(((1:nlambda.ext)-1)/(nlambda.ext-1))
            } else {
                lambda_values[[i_alpha]] = fit$lambda
            }
            n_lambda = length(lambda_values[[i_alpha]])
            predicted_values_all_lambda = matrix(rep(NA,n_instance*n_run*n_lambda),ncol=n_lambda)
            feature_coef_per_lambda = vector("list",n_lambda)
            feature_freq_per_lambda = vector("list",n_lambda)
            for (i_lambda in 1:n_lambda) {
                feature_coef_per_lambda[[i_lambda]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
                feature_freq_per_lambda[[i_lambda]] = matrix(rep(NA,n_feature*n_run),ncol=n_run)
            }
            foldid_per_run[[i_alpha]] = matrix(rep(NA,n_instance*n_run),ncol=n_run)
            pb <- progress_bar$new(
            format = paste0("  running MODEL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
            total = n_run, clear = T, width= 60)
            for (i_run in 1:n_run) {
                foldid_per_run[[i_alpha]][,i_run] = sample(foldid)
                feature_coef_per_run = vector("list",n_lambda)
                for (i_lambda in 1:n_lambda) {
                    feature_coef_per_run[[i_lambda]] = matrix(rep(NA,n_feature*n_fold),ncol=n_fold)
                }
                for (i_fold in 1:n_fold) {
                    instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                    fit = glmnet(x[instance_in_bag,],y[instance_in_bag,],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]])
                    n_lambda_eff = length(fit$lambda) # n_lambda_eff may be smaller than n_lambda
                    predicted_values_all_lambda[n_instance*(i_run-1)+which(!instance_in_bag),1:n_lambda_eff] = predict(fit, x[!instance_in_bag,], type=prediction_type)
                    for (i_lambda in 1:n_lambda_eff) {
                        #feature_coef_per_run[[i_lambda]][,i_fold] = coef(fit)[-1,][,i_lambda] # we remove the intercept
                        feature_coef_per_run[[i_lambda]][,i_fold] = coef(fit)[,i_lambda]
                    }
                }
                for (i_lambda in 1:n_lambda) {
                    is.na(feature_coef_per_run[[i_lambda]]) <- feature_coef_per_run[[i_lambda]] == 0
                    feature_coef_per_lambda[[i_lambda]][,i_run] = rowMeans(feature_coef_per_run[[i_lambda]], na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                    feature_freq_per_lambda[[i_lambda]][,i_run] = rowSums(!is.na(feature_coef_per_run[[i_lambda]]))/n_fold
                }
                pb$tick()
            }
            pb <- progress_bar$new(
            format = paste0("  MODEL predictions for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
            total = n_lambda, clear = T, width= 60)
            lambda_QF_est_all_runs = matrix(rep(NA,n_lambda*n_run),ncol=n_run)
            for (i_lambda in 1:n_lambda) {
                for (i_run in 1:n_run) {
                    predicted_OOB = predicted_values_all_lambda[(1:n_instance)+(i_run-1)*n_instance,i_lambda]
                    if (sum(is.na(predicted_OOB))==0) {
                        lambda_QF_est_all_runs[i_lambda,i_run] = QF(predicted_OOB,y)
                    }
                }
                pb$tick()
            }
            if (save_lambda_QF_full) {
                lambda_QF_est_full[[i_alpha]] = lambda_QF_est_all_runs
            }
            lambda_QF_est[[i_alpha]] = apply(lambda_QF_est_all_runs,1,median,na.rm=T)
            best_lambda_index = which.max(lambda_QF_est[[i_alpha]])
            
            best_lambda[i_alpha] = lambda_values[[i_alpha]][best_lambda_index]
            model_QF_est[i_alpha] = lambda_QF_est[[i_alpha]][best_lambda_index]
            
            predicted_values[[i_alpha]] = matrix(rep(NA,n_instance*2),ncol=2)
            for (i_instance in 1:n_instance) {
                predicted_values[[i_alpha]][i_instance,1] = median(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
                best_lambda_index],na.rm=T)
                predicted_values[[i_alpha]][i_instance,2] = mad(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
                best_lambda_index],na.rm=T)
            }
            
            weight_tot = rowSums(feature_freq_per_lambda[[best_lambda_index]])
            feature_coef_wmean[,i_alpha] = rowSums(feature_coef_per_lambda[[best_lambda_index]]*feature_freq_per_lambda[[best_lambda_index]],na.rm=T)/weight_tot
            for (i_feature in 1:n_feature) {
                feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(feature_freq_per_lambda[[best_lambda_index]][i_feature,]/weight_tot[i_feature] * (feature_coef_per_lambda[[best_lambda_index]][i_feature,] - feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
            }
            feature_freq_mean[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,mean)
            feature_freq_sd[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,sd)
            
            if (logrank) {
                pb <- progress_bar$new(
                format = paste0(" MODEL logrank for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
                total = n_run, clear = T, width= 60)
                logrank_chisq = rep(NA,n_run)
                for (i_run in 1:n_run) {
                    risk_group = risk_group_template
                    risk_group[order(predicted_values_all_lambda[n_instance*(i_run-1)+(1:n_instance),best_lambda_index])[round(n_instance/2+1):n_instance]]=2
                    logrank_chisq[i_run] = survdiff(Surv(time,status)~risk_group,data=data.frame(y,as.factor(risk_group)))$chisq
                    pb$tick()
                }
            }
            
            if (survAUC) {
                pb <- progress_bar$new(
                format = paste0("  MODEL AUC for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
                total = n_run*n_survAUC_time, clear = T, width= 60)
                AUC = vector("list",n_survAUC_time)
                for (i_survAUC_time in 1:n_survAUC_time) {
                    AUC[[i_survAUC_time]] = rep(NA,n_run)
                    for (i_run in 1:n_run) {
                        AUC[[i_survAUC_time]][i_run] = as.numeric(timeROC(T=y[,"time"],delta=y[,"status"],cause=1,marker=predicted_values_all_lambda[n_instance*(i_run-1)+(1:n_instance),best_lambda_index],times=survAUC_time[i_survAUC_time],ROC=FALSE)$AUC[2])
                        pb$tick()
                    }
                    AUC_mean[i_survAUC_time,i_alpha] = mean(AUC[[i_survAUC_time]])
                    AUC_sd[i_survAUC_time,i_alpha] = sd(AUC[[i_survAUC_time]])
                    AUC_perc025[i_survAUC_time,i_alpha] = quantile(AUC[[i_survAUC_time]],probs=0.025)
                    AUC_perc500[i_survAUC_time,i_alpha] = quantile(AUC[[i_survAUC_time]],probs=0.50)
                    AUC_perc975[i_survAUC_time,i_alpha] = quantile(AUC[[i_survAUC_time]],probs=0.975)
                }
            }
            
            # permutation-based null
            null_QF_est_all_runs = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
            null_feature_coef = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
            null_feature_freq = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
            pb <- progress_bar$new(
            format = paste0("  running NULL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
            total = n_run, clear = T, width= 60)
            if (logrank) {
                logrank_chisq_null = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
            }
            if (survAUC) {
                AUC_null = vector("list",n_survAUC_time)
                for (i_survAUC_time in 1:n_survAUC_time) {
                    AUC_null[[i_survAUC_time]] = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
                }
            }
            for (i_run in 1:n_run) {
                for (i_perm_null in 1:n_perm_null) {   # permutation loop
                    y_RDM = y[sample(n_instance),] # we randomly permute the response vector
                    null_predicted_values = rep(NA,n_instance)
                    null_feature_coef_per_run = matrix(rep(NA,n_feature*n_fold),ncol=n_fold)
                    for (i_fold in 1:n_fold) {
                        instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                        fit = glmnet(x[instance_in_bag,],y_RDM[instance_in_bag,],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]][best_lambda_index])
                        null_predicted_values[which(!instance_in_bag)] = predict(fit, x[!instance_in_bag,], s=lambda_values[[i_alpha]][best_lambda_index], type=prediction_type)
                        #null_feature_coef_per_run[,i_fold] = coef(fit)[-1,] # we remove the intercept
                        null_feature_coef_per_run[,i_fold] = as.numeric(coef(fit))
                    }
                    null_QF_est_all_runs[i_run,i_perm_null] = QF(null_predicted_values,y_RDM)
                    is.na(null_feature_coef_per_run) <- null_feature_coef_per_run == 0
                    null_feature_coef[,(i_run-1)*n_perm_null+i_perm_null] = rowMeans(null_feature_coef_per_run, na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                    null_feature_freq[,(i_run-1)*n_perm_null+i_perm_null] = rowSums(!is.na(null_feature_coef_per_run))/n_fold
                    if (logrank) {
                        risk_group = risk_group_template
                        risk_group[order(null_predicted_values)[round(n_instance/2+1):n_instance]]=2
                        logrank_chisq_null[i_run,i_perm_null] = survdiff(Surv(time,status)~risk_group,data=data.frame(y_RDM,as.factor(risk_group)))$chisq
                    }
                    if (survAUC) {
                        for (i_survAUC_time in 1:n_survAUC_time) {
                            AUC_null[[i_survAUC_time]][i_run,i_perm_null] =
                            as.numeric(timeROC(T=y_RDM[,"time"],delta=y_RDM[,"status"],cause=1,marker=null_predicted_values,times=survAUC_time[i_survAUC_time],ROC=FALSE)$AUC[2])
                        }
                    }
                }
                pb$tick()
            }
            weight_tot = rowSums(null_feature_freq)
            null_feature_coef_wmean[,i_alpha] = rowSums(null_feature_coef*null_feature_freq,na.rm=T)/weight_tot
            for (i_feature in 1:n_feature) {
                null_feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(null_feature_freq[i_feature,]/weight_tot[i_feature] * (null_feature_coef[i_feature,] - null_feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
            }
            null_feature_freq_mean[,i_alpha] = apply(null_feature_freq,1,mean)
            null_feature_freq_sd[,i_alpha] = apply(null_feature_freq,1,sd)
            
            # Comparisons of model vs null
            n_tail = 0
            n_tot = 0
            for (i_run in 1:n_run) {
                null_above_model = null_QF_est_all_runs[i_run,]>lambda_QF_est_all_runs[best_lambda_index,i_run]
                n_tail = n_tail + sum(null_above_model,na.rm=T)
                n_tot = n_tot + sum(!is.na(null_above_model))
            }
            QF_model_vs_null_pval[i_alpha] = (n_tail+1)/(n_tot+1) # correction based on Phipson & Smyth (2010)
            
            if (logrank) {
                n_tail = 0
                n_tot = 0
                for (i_run in 1:n_run) {
                    null_above_model = logrank_chisq_null[i_run,]>logrank_chisq[i_run]
                    n_tail = n_tail + sum(null_above_model,na.rm=T)
                    n_tot = n_tot + sum(!is.na(null_above_model))
                }
                logrank_pval[i_alpha] = (n_tail+1)/(n_tot+1)
            }
            
            if (survAUC) {
                for (i_survAUC_time in 1:n_survAUC_time) {
                    n_tail = 0
                    n_tot = 0
                    for (i_run in 1:n_run) {
                        null_above_model = AUC_null[[i_survAUC_time]][i_run,]>AUC[[i_survAUC_time]][i_run]
                        n_tail = n_tail + sum(null_above_model,na.rm=T)
                        n_tot = n_tot + sum(!is.na(null_above_model))
                    }
                    AUC_pval[i_survAUC_time,i_alpha] = (n_tail+1)/(n_tot+1)
                }
            }
            
            for (i_feature in 1:n_feature) {
                n_coef = 0
                n_tail_coef = 0
                n_tail_freq = 0
                for (i_run in 1:n_run) {
                    coef_null_vs_model = abs(null_feature_coef[i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)])>=abs(feature_coef_per_lambda[[best_lambda_index]][i_feature,i_run])
                    n_coef = n_coef + sum(!is.na(coef_null_vs_model))
                    n_tail_coef = n_tail_coef + sum(coef_null_vs_model,na.rm=T)
                    n_tail_freq = n_tail_freq + sum(null_feature_freq[i_feature,(i_run-1)*n_perm_null+(1:n_perm_null)]>=feature_freq_per_lambda[[best_lambda_index]][i_feature,i_run])
                }
                feature_coef_model_vs_null_pval[i_feature,i_alpha] = (n_tail_coef+1)/(n_coef+1)
                feature_freq_model_vs_null_pval[i_feature,i_alpha] = (n_tail_freq+1)/(n_run*n_perm_null+1)
            }
            
            n_alpha_eff = i_alpha # to keep track of fully completed cycles
            cat("alpha = ",alpha[i_alpha]," completed.\n")
        } # end of alpha loop
        return_obj(n_alpha_eff)
    }, error = function(e) {
        cat("Error during execution of alpha = ",alpha[i_alpha],"\n")
        if (n_alpha_eff>0) {
            return_obj(n_alpha_eff)
        }
    }) # end of tryCatch
}

