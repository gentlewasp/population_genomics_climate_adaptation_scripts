draw.tree <- function(tree_file,n,summary_tau_file = NA,lwd.edge = 2.5,col.node = "white",edges = TRUE,leafs = FALSE,nodes = TRUE) {
	if (leafs & nodes) {
		stop("edges vs. nodes: choose between one of the two options")
	}
	if (!file.exists(tree_file)) {
		stop(paste("The file ",tree_file," does not exist. Check the input file name...",sep = ""))
	}
	tree <- scan(tree_file)
	ancstors <- as.vector(tree,mode = "numeric")
	if (ancstors[length(ancstors)] != 0) {
		stop("The last node in the sequence is not the root node (i.e., the last number in the sequence should be zero)... Check the input file")
	}
	ancstors <- ancstors[ancstors > 0]
	nbr.nodes <- length(ancstors)
	descdnts <- seq(1,nbr.nodes)
	nbr.leafs <- min(ancstors) - 1
	if (length(ancstors) > (2 * n - 1)) {
		stop(paste("The total number of nodes in the sequence cannot exceed (2n - 1), where n = ",toString(n)," is the number of sampled populations... Check the input file",sep = ""))
	}
	test <- which(ancstors > (2 * n - 1))
	if (length(test) > 0) {
		stop(paste("The node ",toString(test[1])," in the sequence cannot be larger than (2n - 1), where n = ",toString(n)," is the number of sampled populations... Check the input file",sep = ""))
	}
	test <- which(ancstors <= n)
	if (length(test) > 0) {
		stop(paste("The node ",toString(test[1])," in the sequence cannot be lower than or equal to n, where n = ",toString(n)," is the number of sampled populations... Check the input file",sep = ""))
	}
	
	if (nbr.leafs != n) {
		stop("The number of terminal nodes in the tree does not match the number of populations in the data file... Check the input file")
	}
	test <- which(ancstors == descdnts)
	if (length(test) > 0) {
		stop(paste("The node ",toString(test[1])," in the sequence cannot be its own ancestor... Check the input file",sep = ""))
	}
	test <- which(ancstors < descdnts)
	if (length(test) > 0) {
		stop(paste("The ancestor of node ",toString(test[1])," in the sequence must be larger than ",toString(test[1]),"... Check the input file",sep = ""))
	}
	ancst.nodes <- sort(unique(ancstors))
	descd.nodes <- list()
	idx <- 0
	X <- rep(0,(nbr.nodes + 1))
	for (i in 1:length(ancst.nodes)) {
		totl.dghtrs <- sort(descdnts[ancstors == ancst.nodes[i]])
 		descd.nodes[[ ancst.nodes[i] ]] <- totl.dghtrs 		
 		leaf.dghtrs <- totl.dghtrs[totl.dghtrs %in% 1:nbr.leafs]
  		if (length(leaf.dghtrs) > 0) {
  			X[leaf.dghtrs] <- (idx + 1):(idx + length(leaf.dghtrs))
  			idx <- idx + length(leaf.dghtrs)
  		}
	}
	for (i in ancst.nodes) {
		X[i] <- mean(X[descd.nodes[[ i ]]])	
	}
	if (is.na(summary_tau_file)) {
		tau <- rep(1,nbr.nodes)
	} else {
 		tau <- read.table(summary_tau_file,header = TRUE)$mean
	}
	dstnce.root <- c(tau,0)
	for (i in rev(ancst.nodes)) {
		dstnce.root[descd.nodes[[ i ]]] <- dstnce.root[descd.nodes[[ i ]]] + dstnce.root[i]
	}
	Y <- max(dstnce.root) - dstnce.root
	if (is.na(summary_tau_file)) {
		Y[1:nbr.leafs] <- 0
	}
	if (is.na(summary_tau_file)) {
		plot(X,Y,frame.plot = FALSE,xaxt = "n",yaxt = "n",xlab ="",ylab = "",type = "n",ylim = c(((min(Y) - max(Y)) / 20),max(Y)))
	} else {
		plot(X,Y,frame.plot = FALSE,xaxt = "n",xlab ="",ylab = "Branch length",type = "n",cex.lab = 1.5,ylim = c(((min(Y) - max(Y)) / 20),max(Y)))
	}
	for (i in ancst.nodes) {
		for (j in descd.nodes[[ i ]]) {
			lines(c(X[j],X[j]),c(Y[i],Y[j]),lwd = lwd.edge)
			lines(c(X[i],X[j]),c(Y[i],Y[i]),lwd = lwd.edge)
			if (edges) {
				points(X[j],mean(c(Y[i],Y[j])),col ="white",pch = 16,cex = 3)
				text(X[j],mean(c(Y[i],Y[j])),j)
			}
		}
	}
	if (leafs) {
		text(X[1:nbr.leafs],Y[1:nbr.leafs],(1:nbr.leafs),pos = 1,cex = 1.5)
	}
	if (nodes) {
		points(X,Y,col = col.node,pch = 16,cex = 4.5)
		points(X,Y,col = "black",pch = 1,cex = 4.5)
		text(X,Y,(1:(nbr.nodes + 1)))
	}
}




