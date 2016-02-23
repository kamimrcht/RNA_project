## 1 v 1
# t <- read.table("~/RNA_project/out.txt")
# names(t) <- c("windowsInTarget", "windowsInQuery", "score")

# i <- 0
# while (length(t[which(t$windowsInQuery==i), 1]) == 0) {
# 	i <- i+1
# }

# if (i != max(t$windowsInTarget)+1 ){
# 	firstVec <- t[which(t$windowsInQuery==i), 1]
# }

# for (i in 0:max(t$windowsInQuery)){
# 	if (i == 0){
# 		if (length(t[which(t$windowsInQuery==0), 1])!=0) {
# 			vec <- data.frame(t[which(t$windowsInQuery==0), 1])
# 			colnames(vec) <- c("windowT")
# 			vec["windowQ0"] <- t[which(t$windowsInQuery==0), 3]
# 		} else {
# 			vec <- data.frame(firstVec)
# 			colnames(vec) <- c("windowT")
# 			vec["windowQ0"] <- rep(0,max(t$windowsInTarget)+1)
# 		}
#     } else {
#     	name <- paste("windowQ",i,sep="")
#     	if (length(t[which(t$windowsInQuery==i), 1])!=0) {
# 	    	vec[name] <- t[which(t$windowsInQuery==i),3]
# 	    } else {
# 	    	vec[name] <- rep(0,max(t$windowsInTarget)+1)
# 	    }
#     }
# }
# rownames(vec) <- paste(rep("windowT",length(vec$windowT)),vec$windowT, sep ="")
# vec <- vec[order(vec$windowT),]
# vv <- vec[,2:dim(vec)[2]]
# vvv <- data.matrix(vv)
# vvv_heatmap <- heatmap(vvv, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

## 1 v all
t <- read.table("out.txt")
names(t) <- c("windowsInTarget", "windowsInQuery", "score")
i <- 0
while (length(t[which(t$windowsInQuery==i), 1]) == 0) {
	i <- i+1
}
if (i != max(t$windowsInTarget)+1 ){
	firstVec <- t[which(t$windowsInQuery==i), -2]
	firstVec <- aggregate(firstVec$score, by = list(firstVec$windowsInTarget), mean)
	firstVec <- firstVec[,1]

}
for (i in 0:max(t$windowsInQuery)){
	if (i == 0){
		if (length(t[which(t$windowsInQuery==0), 1])!=0) {
			firstVec <- t[which(t$windowsInQuery==0), -2]
			firstVec <- aggregate(firstVec$score, by = list(firstVec$windowsInTarget), mean)
			vec <- data.frame(firstVec[,1])
			colnames(vec) <- c("windowT")
			vec["windowQ0"] <- firstVec[,2]
		} else {
			vec <- data.frame(firstVec)
			colnames(vec) <- c("windowT")
			vec["windowQ0"] <- rep(0,length(firstVec))
		}
    } else {
    	name <- paste("windowQ",i,sep="")
    	if (length(t[which(t$windowsInQuery==i), 1])!=0) {
	    	addVec <- t[which(t$windowsInQuery==i), -2]
			addVec <- aggregate(addVec$score, by = list(addVec$windowsInTarget), mean)
			vec[name] <- addVec[,2]
	    } else {
	    	vec[name] <- rep(0,length(firstVec))
	    }
    }
}
rownames(vec) <- paste(rep("windowT",length(vec$windowT)),vec$windowT, sep ="")
vec <- vec[order(vec$windowT),]
vv <- vec[,2:dim(vec)[2]]
vvv <- data.matrix(vv)
vvv_heatmap <- heatmap(vvv, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))



