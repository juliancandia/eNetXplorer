plot.feature.heatmap.eNetXplorer <- function (
x, alpha.index, stat, feature.all = F, feature.pval.thres = NULL, feature.set = NULL,
feature.top.n = 25, signif.code = T, xlab=NULL, ylab=NULL, main=NULL, col.main="black", cex.main=0.95,
line=1, col=NULL, breaks=NULL, scale="none", Rowv=F, Colv=F, na.color=NULL,
cexRow = NULL, srtRow=0, cexCol=0.75, srtCol=45, margins=c(5,5), key=T, key.title=NA, dendogram="none",
trace="none", notecol.freq = "black", notecol.coef = "white", notecex = 1, subtitle1=NULL, col.subtitle1="black",
line.subtitle1=-1,cex.subtitle1=0.65, subtitle2=NULL, col.subtitle2="darkgray",line.subtitle2=-2,cex.subtitle2=0.55,
...) {
    if ((x$family=="multinomial")&&(sum(levels(x$response)>2))) {
        plot.feature.heatmap.multinom.eNetXplorer(
        x=x, alpha.index=alpha.index, stat=stat, feature.all=feature.all,
        feature.pval.thres=feature.pval.thres, feature.set=feature.set,
        feature.top.n=feature.top.n, signif.code=signif.code,
        xlab=xlab, ylab=ylab, main=main, col.main=col.main, cex.main=cex.main, line=line,
        col=col, breaks=breaks, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color,
        cexRow=cexRow, srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, margins=margins,
        key=key, key.title=key.title, dendogram=dendogram, trace=trace,
        notecol.freq=notecol.freq, notecol.coef=notecol.coef, notecex=notecex, subtitle1=subtitle1,
        col.subtitle1=col.subtitle1, line.subtitle1=line.subtitle1, cex.subtitle1=cex.subtitle1,
        subtitle2=subtitle2, col.subtitle2=col.subtitle2, line.subtitle2=line.subtitle2, cex.subtitle2=cex.subtitle2,
        ...)
    } else {
        n_alpha = length(x$alpha)
        for (i_alpha in alpha.index) {
            if (stat=="freq") {
                feature_heatmap = as.matrix(x$feature_freq_mean)
                pval = x$feature_freq_model_vs_null_pval
                main.def = "frequencies"
                notecol = notecol.freq
                if (x$family=="binomial") {
                    class = levels(x$response)
                }
            } else if (stat=="coef") {
                feature_heatmap = as.matrix(x$feature_coef_wmean)
                pval = x$feature_coef_model_vs_null_pval
                main.def = "coefficients"
                notecol = notecol.coef
                if (x$family=="binomial") {
                    class = levels(x$response)
                    main.def = paste0(main.def," for class ",class[2])
                }
            }
            n_feature = length(x$feature)
            o = order(pval[,i_alpha])
            if (!feature.all) {
                if (!is.null(feature.pval.thres)) {
                    n_feature = sum(pval[o,i_alpha]<feature.pval.thres)
                } else if (!is.null(feature.set)) {
                    feature_select = which(x$feature%in%feature.set)
                    n_feature = length(feature_select)
                    o = feature_select[order(pval[feature_select,i_alpha])]
                } else {
                    n_feature = min(n_feature,feature.top.n)
                }
            }
            if (n_feature<length(x$feature)) {
                main.def = paste("Selected feature",main.def)
            } else {
                main.def = paste("Feature",main.def)
            }
            if (n_feature>1) {
                o = o[1:n_feature]
                feature_heatmap = feature_heatmap[o,]
                colnames(feature_heatmap) = paste0("alpha = ",x$alpha)
                rownames(feature_heatmap) = x$feature[o]
                if (is.null(main)) {
                    main.title = main.def
                } else {
                    main.title = main
                }
                if (is.null(breaks)) {
                    if (stat=="freq") {
                        breaks = seq(0,1,by=0.2)
                    } else if (stat=="coef") {
                        scale_max = quantile(abs(feature_heatmap),probs=0.95,na.rm=T)
                        scale_int = scale_max/5
                        scale_breaks = seq(scale_int,scale_max,by=scale_int)
                        breaks = c(-rev(scale_breaks),0,scale_breaks)
                    }
                }
                if (is.null(col)) {
                    if (stat=="freq") {
                        col = colorRampPalette(brewer.pal(length(breaks)-1,"Blues"))
                    } else if (stat=="coef") {
                        col = redgreen
                    }
                }
                if (is.null(na.color)) {
                    if (stat=="freq") {
                        na.color = "black"
                    } else if (stat=="coef") {
                        na.color = "white"
                    }
                }
                if (is.null(cexRow)) {
                    cexRow = min(1.6/log(n_feature),0.5)
                }
                if (signif.code) {
                    x$feature_coef_model_vs_null_pval[o,]
                    annot = matrix(rep("",n_feature*n_alpha),ncol=n_alpha)
                    annot[which(pval[o,]<0.1,arr.ind=T)] = "."
                    annot[which(pval[o,]<0.05,arr.ind=T)] = "*"
                    annot[which(pval[o,]<0.01,arr.ind=T)] = "**"
                    annot[which(pval[o,]<0.001,arr.ind=T)] = "***"
                    heatmap.2(feature_heatmap, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color, col=col,
                    breaks = breaks,dendrogram = dendogram, margins=margins, cexRow=cexRow,
                    srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, key=key, trace=trace, key.title=key.title,
                    xlab= xlab, ylab = ylab,
                    cellnote = annot, notecol = notecol, notecex = notecex, ...)
                    title(main=main.title, col.main=col.main, line=line, cex.main=cex.main)
                    if (is.null(subtitle1)) {
                        main.subtitle1 = paste0("Features ranked by p-value for alpha=",x$alpha[i_alpha])
                    } else {
                        main.subtitle1 = subtitle1
                    }
                    title(main=main.subtitle1,col.main=col.subtitle1, line=line.subtitle1, cex.main=cex.subtitle1)
                    if (is.null(subtitle2)) {
                        main.subtitle2 = "P-value significance codes:  <0.001 (***), <0.01 (**), <0.05 (*), <0.1 (.)"
                    } else {
                        main.subtitle2 = subtitle2
                    }
                    title(main=main.subtitle2,col.main=col.subtitle2, line=line.subtitle2, cex.main=cex.subtitle2)
                } else {
                    heatmap.2(feature_heatmap, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color, col=col,
                    breaks = breaks,dendrogram = dendogram, margins=margins, cexRow=cexRow,
                    srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, key=key, trace=trace,
                    xlab= xlab, ylab = ylab,
                    ...)
                    title(main=main.title, col.main=col.main, cex.main=cex.main)
                }
            } else {
                warning("Number of features selected is < 2")
            }
        }
    }
}
