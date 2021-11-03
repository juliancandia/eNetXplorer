mergeObj <- function(source_obj, source_dir=getwd(), dest_obj="eNet_merged.Robj", dest_dir=NULL)
{
    n_obj = length(source_obj)
    if (n_obj<2) {stop("Merge function requires two or more source objects")}
    if (!is.null(source_dir)) {
        source_obj = file.path(source_dir,source_obj)
    }
    if (is.null(dest_dir)) {
        if (!is.null(source_dir)) {
            dest_dir = source_dir
        } else {
            dest_dir = getwd()
        }
    }
    fit = vector("list",n_obj)
    for (i_obj in 1:n_obj) {
        load(source_obj[i_obj])
        fit[[i_obj]] = eNet
    }
    # 1. check that all runs are compatible.
    family = fit[[1]]$family
    if (family=="multinomial") {
        stop("Not supported yet") # This feature will be added later...
    }
    check_passed = T
    for (i_obj in 2:n_obj) {
        check_passed = check_passed&identical(fit[[i_obj]]$predictor,fit[[1]]$predictor)&
        identical(fit[[i_obj]]$response,fit[[1]]$response)&
        identical(fit[[i_obj]]$family,fit[[1]]$family)&
        identical(fit[[i_obj]]$nlambda,fit[[1]]$nlambda)&
        identical(fit[[i_obj]]$nlambda.ext,fit[[1]]$nlambda.ext)&
        identical(fit[[i_obj]]$scaled,fit[[1]]$scaled)&
        identical(fit[[i_obj]]$n_fold,fit[[1]]$n_fold)&
        identical(fit[[i_obj]]$n_run,fit[[1]]$n_run)&
        identical(fit[[i_obj]]$n_perm_null,fit[[1]]$n_perm_null)&
        identical(fit[[i_obj]]$QF_label,fit[[1]]$QF_label)&
        identical(fit[[i_obj]]$instance,fit[[1]]$instance)&
        identical(fit[[i_obj]]$feature,fit[[1]]$feature)&
        identical(fit[[i_obj]]$glmnet_params,fit[[1]]$glmnet_params)
        if (family=="gaussian") {
            check_passed = check_passed&identical(fit[[i_obj]]$cor_method,fit[[1]]$cor_method)
        }
        if (family=="binomial") {
            check_passed = check_passed&identical(fit[[i_obj]]$binom_method,fit[[1]]$binom_method)&
            identical(fit[[i_obj]]$binom_pos,fit[[1]]$binom_pos)&
            identical(fit[[i_obj]]$fscore_beta,fit[[1]]$fscore_beta)&
            identical(fit[[i_obj]]$fold_distrib_fail.max,fit[[1]]$fold_distrib_fail.max)
        }
        if (family=="cox") {
            check_passed = check_passed&identical(fit[[i_obj]]$cox_index,fit[[1]]$cox_index)&
            identical(fit[[i_obj]]$logrank,fit[[1]]$logrank)&
            identical(fit[[i_obj]]$survAUC,fit[[1]]$survAUC)&
            identical(fit[[i_obj]]$survAUC_time,fit[[1]]$survAUC_time)
        }
    }
    if (!check_passed) {stop("Source objects are incompatible")}

    predictor = fit[[1]]$predictor
    response = fit[[1]]$response
    nlambda = fit[[1]]$nlambda
    nlambda.ext = fit[[1]]$nlambda.ext
    scaled = fit[[1]]$scaled
    n_fold = fit[[1]]$n_fold
    n_run = fit[[1]]$n_run
    n_perm_null = fit[[1]]$n_perm_null
    QF_label = fit[[1]]$QF_label
    instance = fit[[1]]$instance
    feature = fit[[1]]$feature
    glmnet_params = fit[[1]]$glmnet_params
    n_feature = length(feature)

    # we merge the alphas from all objects
    alpha = NULL
    index_obj = NULL
    index_alpha = NULL
    for (i_obj in 1:n_obj) {
        alpha = c(alpha,fit[[i_obj]]$alpha)
        index_obj = c(index_obj,rep(i_obj,length(fit[[i_obj]]$alpha)))
        index_alpha = c(index_alpha,1:length(fit[[i_obj]]$alpha))
    }
    select = which(!duplicated(alpha))
    select = select[order(alpha[select])]
    alpha = alpha[select]
    n_alpha = length(alpha)
    index_obj = index_obj[select]
    index_alpha = index_alpha[select]

    # we merge objects as needed (these are all family-independent objects)
    best_lambda = rep(NA,n_alpha)
    names(best_lambda) = alpha
    model_QF_est = rep(NA,n_alpha)
    names(model_QF_est) = alpha
    QF_model_vs_null_pval = rep(NA,n_alpha)
    names(QF_model_vs_null_pval) = alpha
    lambda_values = rep(NA,n_alpha)
    names(lambda_values) = alpha
    lambda_QF_est = rep(NA,n_alpha)
    names(lambda_QF_est) = alpha
    predicted_values = vector("list",n_alpha)
    names(predicted_values) = alpha
    feature_coef_wmean = NULL
    feature_coef_wsd = NULL
    feature_freq_mean = NULL
    feature_freq_sd = NULL
    null_feature_coef_wmean = NULL
    null_feature_coef_wsd = NULL
    null_feature_freq_mean = NULL
    null_feature_freq_sd = NULL
    feature_coef_model_vs_null_pval = NULL
    feature_freq_model_vs_null_pval = NULL
    for (i_alpha in 1:n_alpha) {
        best_lambda[i_alpha] = fit[[index_obj[i_alpha]]]$best_lambda[index_alpha[i_alpha]]
        model_QF_est[i_alpha] = fit[[index_obj[i_alpha]]]$model_QF_est[index_alpha[i_alpha]]
        QF_model_vs_null_pval[i_alpha] = fit[[index_obj[i_alpha]]]$QF_model_vs_null_pval[index_alpha[i_alpha]]
        lambda_values[i_alpha] = fit[[index_obj[i_alpha]]]$lambda_values[index_alpha[i_alpha]]
        lambda_QF_est[i_alpha] = fit[[index_obj[i_alpha]]]$lambda_QF_est[index_alpha[i_alpha]]
        predicted_values[[i_alpha]] = (fit[[index_obj[i_alpha]]]$predicted_values)[[index_alpha[i_alpha]]]
        feature_coef_wmean = cbind(feature_coef_wmean,as.matrix(fit[[index_obj[i_alpha]]]$feature_coef_wmean)[,index_alpha[i_alpha]])
        feature_coef_wsd = cbind(feature_coef_wsd,as.matrix(fit[[index_obj[i_alpha]]]$feature_coef_wsd)[,index_alpha[i_alpha]])
        feature_freq_mean = cbind(feature_freq_mean,as.matrix(fit[[index_obj[i_alpha]]]$feature_freq_mean)[,index_alpha[i_alpha]])
        feature_freq_sd = cbind(feature_freq_sd,as.matrix(fit[[index_obj[i_alpha]]]$feature_freq_sd)[,index_alpha[i_alpha]])
        null_feature_coef_wmean = cbind(null_feature_coef_wmean,as.matrix(fit[[index_obj[i_alpha]]]$null_feature_coef_wmean)[,index_alpha[i_alpha]])
        null_feature_coef_wsd = cbind(null_feature_coef_wsd,as.matrix(fit[[index_obj[i_alpha]]]$null_feature_coef_wsd)[,index_alpha[i_alpha]])
        null_feature_freq_mean = cbind(null_feature_freq_mean,as.matrix(fit[[index_obj[i_alpha]]]$null_feature_freq_mean)[,index_alpha[i_alpha]])
        null_feature_freq_sd = cbind(null_feature_freq_sd,as.matrix(fit[[index_obj[i_alpha]]]$null_feature_freq_sd)[,index_alpha[i_alpha]])
        feature_coef_model_vs_null_pval = cbind(feature_coef_model_vs_null_pval,as.matrix(fit[[index_obj[i_alpha]]]$feature_coef_model_vs_null_pval)[,index_alpha[i_alpha]])
        feature_freq_model_vs_null_pval = cbind(feature_freq_model_vs_null_pval,as.matrix(fit[[index_obj[i_alpha]]]$feature_freq_model_vs_null_pval)[,index_alpha[i_alpha]])
    }
    rownames(feature_coef_wmean) = feature
    colnames(feature_coef_wmean) = alpha
    rownames(feature_coef_wsd) = feature
    colnames(feature_coef_wsd) = alpha
    rownames(feature_freq_mean) = feature
    colnames(feature_freq_mean) = alpha
    rownames(feature_freq_sd) = feature
    colnames(feature_freq_sd) = alpha
    rownames(null_feature_coef_wmean) = feature
    colnames(null_feature_coef_wmean) = alpha
    rownames(null_feature_coef_wsd) = feature
    colnames(null_feature_coef_wsd) = alpha
    rownames(null_feature_freq_mean) = feature
    colnames(null_feature_freq_mean) = alpha
    rownames(null_feature_freq_sd) = feature
    colnames(null_feature_freq_sd) = alpha
    rownames(feature_coef_model_vs_null_pval) = feature
    colnames(feature_coef_model_vs_null_pval) = alpha
    rownames(feature_freq_model_vs_null_pval) = feature
    colnames(feature_freq_model_vs_null_pval) = alpha

    eNet = list(
    # input data and parameters
    predictor = predictor, response = response, alpha = alpha, family = family, nlambda = nlambda,
    nlambda.ext = nlambda.ext, seed = NA, scaled = scaled, n_fold = n_fold, n_run = n_run,
    n_perm_null = n_perm_null,
    QF_label = QF_label, instance = instance,
    feature = feature, glmnet_params = glmnet_params,
    # the list below needs to be merged... Different families will require different outputs!!!
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

    # Family-dependent objects are added here
    if (family=="gaussian") {
        eNet = append(eNet, list(cor_method = fit[[1]]$cor_method))
    }
    if (family=="binomial") {
        eNet = append(eNet, list(binom_method = fit[[1]]$binom_method, binom_pos = fit[[1]]$binom_pos, fscore_beta = fit[[1]]$fscore_beta, fold_distrib_fail.max = fit[[1]]$fold_distrib_fail.max))
    }
    if (family=="cox") {
        eNet = append(eNet, list(cox_index = fit[[1]]$cox_index, logrank = fit[[1]]$logrank, survAUC = fit[[1]]$survAUC,
        survAUC_time = fit[[1]]$survAUC_time))
        if (fit[[1]]$logrank) {
            logrank_pval = rep(NA,n_alpha)
            for (i_alpha in 1:n_alpha) {
                logrank_pval[i_alpha] = fit[[index_obj[i_alpha]]]$logrank_pval[index_alpha[i_alpha]]
            }
            names(logrank_pval) = alpha
            eNet = append(eNet, list(logrank_pval = logrank_pval))
        }
        if (fit[[1]]$survAUC) {
            AUC_mean = NULL
            AUC_sd = NULL
            AUC_perc025 = NULL
            AUC_perc500 = NULL
            AUC_perc975 = NULL
            AUC_pval = NULL
            for (i_alpha in 1:n_alpha) {
                AUC_mean = cbind(AUC_mean,fit[[index_obj[i_alpha]]]$AUC_mean[,index_alpha[i_alpha]])
                AUC_sd = cbind(AUC_sd,fit[[index_obj[i_alpha]]]$AUC_sd[,index_alpha[i_alpha]])
                AUC_perc025 = cbind(AUC_perc025,fit[[index_obj[i_alpha]]]$AUC_perc025[,index_alpha[i_alpha]])
                AUC_perc500 = cbind(AUC_perc500,fit[[index_obj[i_alpha]]]$AUC_perc500[,index_alpha[i_alpha]])
                AUC_perc975 = cbind(AUC_perc975,fit[[index_obj[i_alpha]]]$AUC_perc975[,index_alpha[i_alpha]])
                AUC_pval = cbind(AUC_pval,fit[[index_obj[i_alpha]]]$AUC_pval[,index_alpha[i_alpha]])
            }
            rownames(AUC_mean) = fit[[1]]$survAUC_time
            colnames(AUC_mean) = alpha
            rownames(AUC_sd) = fit[[1]]$survAUC_time
            colnames(AUC_sd) = alpha
            rownames(AUC_perc025) = fit[[1]]$survAUC_time
            colnames(AUC_perc025) = alpha
            rownames(AUC_perc500) = fit[[1]]$survAUC_time
            colnames(AUC_perc500) = alpha
            rownames(AUC_perc975) = fit[[1]]$survAUC_time
            colnames(AUC_perc975) = alpha
            rownames(AUC_pval) = fit[[1]]$survAUC_time
            colnames(AUC_pval) = alpha
            eNet = append(eNet, list(
            AUC_mean = AUC_mean, AUC_sd = AUC_sd, AUC_perc025 = AUC_perc025, AUC_perc500 = AUC_perc500, AUC_perc975 = AUC_perc975, AUC_pval = AUC_pval
            ))
        }
    }
    class(eNet) <- "eNetXplorer"
    save(eNet,file=file.path(dest_dir,dest_obj))
    eNet
}

