plot.mob <- function(obj, terminal_panel = NULL, ...) {
  if(is.null(terminal_panel))
    terminal_panel <- if(is.factor(response(obj)[[1]])) node_spineplot else node_scatterplot
  plot.BinaryTree(obj, terminal_panel = terminal_panel, ...)
}

## FIXME plot.BinaryTree: grapcon_generator

node_scatterplot <- function(ctreeobj, which = NULL, col = "black", linecol = "red",
  cex = 0.5, jitter = FALSE, xscale = NULL, yscale = NULL, ylines = 3, id = TRUE, labels = FALSE)
{
    ## extract dependent variable
    y <- response(ctreeobj)
    ynam <- names(y)[1]
    y <- y[[1]]
    y <- if(is.factor(y)) as.numeric(y) - 1 else as.numeric(y)
    if(jitter) y <- jitter(y)

    ## extract regressor matrix
    x <- model.matrix(ctreeobj@tree$model)
    if(is.null(which)) { which <- if(NCOL(x) > 1) 2:NCOL(x) else 1 }
    xnam <- colnames(x)[which]
    x <- x[,which,drop=FALSE]
    k <- NCOL(x)
    
    if(is.null(xscale)) xscale <- apply(x, 2, function(xi) range(xi) + c(-0.1, 0.1) * diff(range(xi)))
        else xscale <- matrix(xscale)
    if(is.null(yscale)) yscale <- range(y) + c(-0.1, 0.1) * diff(range(y))
         
    ## panel function for scatter plots in nodes
    rval <- function(node) {
    
        ## dependent variable setup
	y <- rep.int(y, node$weights)
	yhat <- rep.int(fitted(node$model), node$weights)

        ## viewport setup
        top_vp <- viewport(layout = grid.layout(nrow = 2*k, ncol = 3,
                           widths = unit(c(ylines, 1, 1), c("lines", "null", "lines")),  
			   heights = unit(rep.int(c(2.5, 1), k) - c(1.5, rep.int(0, 2*k-1)),
			                  rep.int(c("lines", "null"), k))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste("node_scatterplot", node$nodeID, sep = ""))
        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = "white", col = NULL))

        ## main title
        top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
        pushViewport(top)
	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")
        grid.text(mainlab)
        popViewport()
	
	for(i in 1:k) {
            plot_vpi <- viewport(layout.pos.col = 2, layout.pos.row = 2*i, xscale = xscale[,i],
	        yscale = yscale, name = paste("node_scatterplot", i, node$nodeID, "plot", sep = ""))
            pushViewport(plot_vpi)
	
            ## regressor setup
	    xi <- rep.int(x[,i], node$weights)
	    oi <- order(xi)

            ## scatterplot
	    grid.points(xi, y, gp = gpar(col = col, cex = cex))
            grid.lines(xi[oi], yhat[oi], default.units = "native", gp = gpar(col = linecol))

            grid.xaxis(at = c(ceiling(xscale[1,i]*10), floor(xscale[2,i]*10))/10)
            grid.yaxis(at = c(ceiling(yscale[1]), floor(yscale[2])))

            if(labels) {
                grid.text(xnam[i], x = unit(0.5, "npc"), y = unit(-2, "lines"))
                grid.text(ynam, y = unit(0.5, "npc"), x = unit(-3, "lines"), rot = 90)
	    }

            grid.rect()
            upViewport()
        }
	upViewport()
    }
	    
    return(rval)
}

node_spineplot <- function(ctreeobj, which = NULL, id = TRUE, ...)
{
    ## re-use spine() from vcd package
    if(!require("vcd")) stop("the `vcd' package is required for spine plots")

    y <- response(ctreeobj)[[1]]
    if(is.factor(y)) y <- factor(y, levels = rev(levels(y)))

    x <- ctreeobj@data@get("input")
    if(is.null(which)) which <- 1
    which <- which[1] ## currently only a single regressor supported
    x <- x[,which]
    
    ### panel function for spine plots in nodes
    rval <- function(node) {
    
        ## parameter setup
	x <- rep(x, node$weights)
	y <- rep(y, node$weights)

	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")

        grid.rect(width = unit(1, "npc") - unit(2.5, "lines"),
	          height = unit(1, "npc") - unit(0.5, "lines"),
	          gp = gpar(fill = "white", col = "white"))
        spine(x, y, xlab = "", ylab = "", 
	      name = paste("node_spineplot", node$nodeID, sep = ""), newpage = FALSE,
	      margins = rep(1.5, 4), ...)
        grid.text(mainlab, y = unit(1, "npc") - unit(1, "lines"))
    }
    
    return(rval)
}

node_cdplot <- function(ctreeobj, which = NULL, id = TRUE, ...)
{
    ## re-use cd_plot() from vcd package
    if(!require("vcd")) stop("the `vcd' package is required for CD plots")

    y <- response(ctreeobj)[[1]]
    if(is.factor(y)) y <- factor(y, levels = rev(levels(y)))

    x <- ctreeobj@data@get("input")
    if(is.null(which)) which <- 1
    which <- which[1] ## currently only a single regressor supported
    x <- as.numeric(x[,which])
    
    ### panel function for spine plots in nodes
    rval <- function(node) {
    
        ## parameter setup
	x <- rep(x, node$weights)
	y <- rep(y, node$weights)

	mainlab <- paste(ifelse(id, paste("Node", node$nodeID, "(n = "), ""),
	                 sum(node$weights), ifelse(id, ")", ""), sep = "")

        grid.rect(width = unit(1, "npc") - unit(2.5, "lines"),
	          height = unit(1, "npc") - unit(0.5, "lines"),
	          gp = gpar(fill = "white", col = "white"))
        cd_plot(x, y, xlab = "", ylab = "", 
	        name = paste("node_spineplot", node$nodeID, sep = ""), newpage = FALSE,
	        margins = rep(1.5, 4), ...)
        grid.text(mainlab, y = unit(1, "npc") - unit(1, "lines"))
    }
    
    return(rval)
}
