plot.contingency.eNetXplorer <- function (x, alpha.index,
xlab="class (true)", ylab="class (predicted)", cex.lab=0.95, main=NULL, col.main = "black", cex.main=0.85,
cex.axis=1, symbol.size.inches=0.5, bg.color="steelblue2", fg.color=NULL, margin=0.2,
frequency.label=T, frequency.label.cex=1, frequency.label.offset=0,
...)
{
    for (i_alpha in alpha.index) {
        class = levels(x$response)
        n_class = length(class)
        contingency = matrix(rep(NA,n_class**2),ncol=n_class) # rows=true, cols=predicted
        for (i_class in 1:n_class) {
            instance_select = x$response==class[i_class]
            contingency[i_class,] = apply(x$predicted_values[[i_alpha]][instance_select,],2,sum)
        }
        
        freq = NULL
        for (i in 1:nrow(contingency)) {
            for (j in 1:ncol(contingency)) {
                freq = rbind(freq,c(i,j,contingency[i,j]))
            }
        }
        
        x.min=1-(n_class-1)*margin
        x.max=n_class+(n_class-1)*margin
        symbols(x=freq[,1],y=freq[,2],circles=sqrt(freq[,3]/pi), inches=symbol.size.inches, bg=bg.color, fg=fg.color,xlab=xlab,ylab=ylab,cex.lab=cex.lab,xaxt="n",yaxt="n",panel.first = abline(h=1:n_class,v=1:n_class,lty=3),xlim=c(x.min,x.max),ylim=c(x.min,x.max))
        axis(1,at=1:n_class,labels=as.character(class),cex.axis=cex.axis)
        axis(2,at=1:n_class,labels=as.character(class),cex.axis=cex.axis)
        if (is.null(main)) {
            main.title = paste0("contingency matrix for alpha=",x$alpha[i_alpha])
        } else {
            main.title = main
        }
        title(main=main.title, col.main=col.main, cex.main=cex.main)
        if (frequency.label) {
            textxy(freq[,1],freq[,2],freq[,3],cex=frequency.label.cex,offset=frequency.label.offset)
        }
    }
}
