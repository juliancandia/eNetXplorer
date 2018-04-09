# Gaussian model
eNetXplorerGaussian <- function(x, y, family, alpha, nlambda, nlambda.ext, seed, scaled,
n_fold, n_run, n_perm_null, QF.FUN, QF_label, cor_method, ...)
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
    prediction_type = "link"
    
    if (is.null(QF.FUN)) {
        QF = function(predicted,response) {
            cor_test = cor.test(predicted,response,method=cor_method)
            return (cor_test$estimate)
        }
        QF_label = paste0("correlation (",cor_method,")")
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
    best_lambda = rep(NA,n_alpha)
    model_QF_est = rep(NA,n_alpha)
    
    lambda_values = vector("list",n_alpha)
    lambda_QF_est = vector("list",n_alpha)
    
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
    
    glmnet.control(factory=T) # resets internal glmnet parameters
    glmnet.control(...) # to allow changes of factory default parameters in glmnet
    
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
            for (i_fold in 1:n_fold) {
                instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                fit = glmnet(x[instance_in_bag,],y[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]])
                predicted_values_all_lambda[n_instance*(i_run-1)+which(instance_out_of_bag),] = predict(fit, x[instance_out_of_bag,], type=prediction_type)
                for (i_lambda in 1:n_lambda) {
                    feature_coef_per_run[[i_lambda]] = cbind(feature_coef_per_run[[i_lambda]],coef(fit)[-1,][,i_lambda]) # we remove the intercept
                }
            }
            for (i_lambda in 1:n_lambda) {
                is.na(feature_coef_per_run[[i_lambda]]) <- feature_coef_per_run[[i_lambda]] == 0
                feature_coef_per_lambda[[i_lambda]][,i_run] = rowMeans(feature_coef_per_run[[i_lambda]], na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                feature_freq_per_lambda[[i_lambda]][,i_run] = rowSums(!is.na(feature_coef_per_run[[i_lambda]]))/n_fold
            }
            pb$tick()
        }
        
        lambda_QF_est_all_runs = matrix(rep(NA,n_lambda*n_run),ncol=n_run)
        for (i_lambda in 1:n_lambda) {
            for (i_run in 1:n_run) {
                lambda_QF_est_all_runs[i_lambda,i_run] = QF(predicted_values_all_lambda[(1:n_instance)+(i_run-1)*n_instance,i_lambda],y)
            }
        }
        lambda_QF_est[[i_alpha]] = apply(lambda_QF_est_all_runs,1,median)
        best_lambda_index = which.max(lambda_QF_est[[i_alpha]])
        
        best_lambda[i_alpha] = lambda_values[[i_alpha]][best_lambda_index]
        model_QF_est[i_alpha] = lambda_QF_est[[i_alpha]][best_lambda_index]
        
        predicted_values[[i_alpha]] = matrix(rep(NA,n_instance*2),ncol=2)
        for (i_instance in 1:n_instance) {
            predicted_values[[i_alpha]][i_instance,1] = median(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
            best_lambda_index])
            predicted_values[[i_alpha]][i_instance,2] = mad(predicted_values_all_lambda[seq(i_instance,n_instance*n_run,by=n_instance),
            best_lambda_index])
        }
        
        weight_tot = rowSums(feature_freq_per_lambda[[best_lambda_index]])
        feature_coef_wmean[,i_alpha] = rowSums(feature_coef_per_lambda[[best_lambda_index]]*feature_freq_per_lambda[[best_lambda_index]],na.rm=T)/weight_tot
        for (i_feature in 1:n_feature) {
            feature_coef_wsd[i_feature,i_alpha] = sqrt(sum(feature_freq_per_lambda[[best_lambda_index]][i_feature,]/weight_tot[i_feature] * (feature_coef_per_lambda[[best_lambda_index]][i_feature,] - feature_coef_wmean[i_feature,i_alpha])^2,na.rm=T))
        }
        feature_freq_mean[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,mean)
        feature_freq_sd[,i_alpha] = apply(feature_freq_per_lambda[[best_lambda_index]],1,sd)
        
        # permutation-based null
        null_QF_est_all_runs = matrix(rep(NA,n_run*n_perm_null),ncol=n_perm_null)
        null_feature_coef = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        null_feature_freq = matrix(rep(NA,n_feature*n_run*n_perm_null),ncol=n_run*n_perm_null)
        pb <- progress_bar$new(
        format = paste0("  running NULL for alpha = ",alpha[i_alpha]," [:bar] :percent in :elapsed"),
        total = n_run, clear = T, width= 60)
        for (i_run in 1:n_run) {
            for (i_perm_null in 1:n_perm_null) {   # permutation loop
                y_RDM = sample(y) # we randomly permute the response vector
                null_predicted_values = rep(NA,n_instance)
                null_feature_coef_per_run = NULL
                for (i_fold in 1:n_fold) {
                    instance_in_bag = foldid_per_run[[i_alpha]][,i_run]!=i_fold
                    instance_out_of_bag = foldid_per_run[[i_alpha]][,i_run]==i_fold
                    fit = glmnet(x[instance_in_bag,],y_RDM[instance_in_bag],alpha=alpha[i_alpha],family=family,lambda=lambda_values[[i_alpha]][best_lambda_index])
                    null_predicted_values[which(instance_out_of_bag)] = predict(fit, x[instance_out_of_bag,], s=lambda_values[[i_alpha]][best_lambda_index], type=prediction_type)
                    null_feature_coef_per_run = cbind(null_feature_coef_per_run,coef(fit)[-1,]) # we remove the intercept
                }
                null_QF_est_all_runs[i_run,i_perm_null] = QF(null_predicted_values,y_RDM)
                is.na(null_feature_coef_per_run) <- null_feature_coef_per_run == 0
                null_feature_coef[,(i_run-1)*n_perm_null+i_perm_null] = rowMeans(null_feature_coef_per_run, na.rm=T) # we obtain the mean (nonzero) coefficients over folds.
                null_feature_freq[,(i_run-1)*n_perm_null+i_perm_null] = rowSums(!is.na(null_feature_coef_per_run))/n_fold
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
        for (i_run in 1:n_run) {
            n_tail = n_tail + sum(null_QF_est_all_runs[i_run,]>lambda_QF_est_all_runs[best_lambda_index,i_run])
        }
        QF_model_vs_null_pval[i_alpha] = (n_tail+1)/(n_run*n_perm_null+1) # correction based on Phipson & Smyth (2010)
        
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
        
        cat("alpha = ",alpha[i_alpha]," completed.\n")
    } # end of alpha loop
    
    # return object
    list(
    # input data and parameters
    predictor = x, response = y, alpha = alpha, family = family, nlambda = nlambda,
    nlambda.ext = nlambda.ext, seed = seed, scaled = scaled, n_fold = n_fold, n_run = n_run,
    n_perm_null = n_perm_null, QF_label = QF_label, cor_method = cor_method, instance = instance,
    feature = feature, glmnet_params = glmnet.control(),
    # summary results
    best_lambda = best_lambda, model_QF_est = model_QF_est, QF_model_vs_null_pval = QF_model_vs_null_pval,
    # detailed results for plots and downstream analysis
    lambda_values = lambda_values, lambda_QF_est = lambda_QF_est,
    predicted_values = predicted_values,
    feature_coef_wmean = as(feature_coef_wmean,"CsparseMatrix"), feature_coef_wsd = as(feature_coef_wsd,"CsparseMatrix"),
    feature_freq_mean = as(feature_freq_mean,"CsparseMatrix"), feature_freq_sd = as(feature_freq_sd,"CsparseMatrix"),
    null_feature_coef_wmean = as(null_feature_coef_wmean,"CsparseMatrix"), null_feature_coef_wsd = as(null_feature_coef_wsd,"CsparseMatrix"),
    null_feature_freq_mean = as(null_feature_freq_mean,"CsparseMatrix"), null_feature_freq_sd = as(null_feature_freq_sd,"CsparseMatrix"),
    feature_coef_model_vs_null_pval = as(feature_coef_model_vs_null_pval,"CsparseMatrix"),
    feature_freq_model_vs_null_pval = as(feature_freq_model_vs_null_pval,"CsparseMatrix")
    
    )
}

