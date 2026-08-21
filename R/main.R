
set.seed(999)

calmu_int = function(Z_int, label, beta,Num_celltype,num){
  Z_int_ts = apply(Z_int[label==1,], 2, mean)
  for (i in 2:length(unique(label))){
    Z_int_ts = rbind(Z_int_ts, apply(Z_int[label==i,], 2, mean))
  }
  Z_int_ts = t(Z_int_ts)
  
  ## initial mu
  mu_int = matrix(0, num, Num_celltype)
  for (i in 1:num){
    mu_int[i,] = t(solve(t(beta)%*%beta + diag(rep(0.00001,Num_celltype)), t(beta)%*%t(t(Z_int_ts[i,]))))
  }
  return(mu_int)
}

col_match =function(A,B){
# Load the clue package

# Calculate the L2 distance matrix
distance=matrix(rep(0,ncol(A)*ncol(B)),ncol=ncol(A), nrow=ncol(B))
for (i in 1:ncol(A)){
    for (j in 1:ncol(B)){
        distance[i, j] = sqrt(sum((A[, i]- B[, j]) ** 2))
        }
    }

# Solve the Linear Sum Assignment Problem (LSAP)
assignment <- solve_LSAP(distance)

reorder=c()
for (i in 1:length(assignment)){
    reorder=append(reorder, which(assignment==i))
    }

# Reorder the columns of A based on the assignment
A_reordered <- A[, reorder]
    
return (A_reordered)
}

#' @title 执行 QR-SIDE 解卷积
#' @description 主函数，用于空间转录组数据的细胞类型解卷积分析
#' @param sp_counts 空间表达矩阵
#' @param pos 位置坐标
#' @param markerframe 标记基因框架
#' @param Num_topic 主题数
#' @param Num_HVG 高变基因数
#' @param dim_embed 嵌入维度
#' @param max_marker 最大标记基因数
#' @param markerflag 标记基因标志
#' @return 解卷积结果
#' @export
run_QRSIDE <- function(sp_counts, pos, markerframe, Num_topic, Num_HVG,dim_embed,max_marker,markerflag) {

time_spot0 <- Sys.time()
   
sce = Createsce(sp_counts,nrow(sp_counts))
sce = spatialPreprocess(sce, n.HVGs=Num_HVG)
PCA = reducedDim(sce, "PCA")
logX = t(logcounts(sce)[rowData(sce)[["is.HVG"]],])
sum(apply(logX, 1, sum))
X2 =  as.matrix(t(counts(sce)[rowData(sce)[["is.HVG"]],]))
    

dt = 3   
############################################################
if (markerflag== FALSE){
    
CT_=unique(markerframe$label)
print('filtered celltype')
print(CT_)
    
Num_celltype=length(CT_)
markermat=list()
for (c in 1:length(CT_)){
    markermat[[c]]=markerframe$gene[markerframe$label==CT_[c]]
    }

gene_length=0
for (m in markermat){
gene_length=gene_length+length(m)
    } 
    
delta <- matrix(0, nrow = gene_length, ncol = Num_celltype)
rho <- matrix(0, nrow = gene_length, ncol = Num_celltype)

for (j in 1:Num_celltype){

    if (j==1){
    m_len=length(unique(intersect(markermat[[j]],colnames(X2))))
    markerList=unique(intersect(markermat[[j]],colnames(X2)))[1:min(m_len,max_marker)]
    print(markerList)
    if (m_len == 0) {
      stop("Error: This marker is empty.")
        }

    row_start <- 1
    row_end <- length(markerList) 

    rho[row_start:row_end,j] <- 1
    delta[row_start:row_end,j] <- dt

    }else{
    m_len=length(setdiff(intersect(markermat[[j]],colnames(X2)),markerList))
    marker_tmp=setdiff(intersect(markermat[[j]],colnames(X2)),markerList)[1:min(m_len,max_marker)]
    markerList=c(markerList, marker_tmp)
    print(marker_tmp)
    if (m_len == 0) {
      stop("Error: This marker is empty.")
        }

    row_start <- row_end+1
    row_end <- row_start+length(marker_tmp)-1

    rho[row_start:row_end,j] <- 1
    delta[row_start:row_end,j] <- dt

        }
    }

rho <- rho[rowSums(rho[, 1:ncol(rho)]) > 0, ]
delta <- delta[rowSums(delta[, 1:ncol(delta)]) > 0, ]

}
############################################################   
if (markerflag== TRUE){    
    
   Num_celltype <- length(markerframe)
   markerList=unlist(markerframe) 
    
   for (i in seq_along(markerframe)){

    if (i==1){
        t_m=matrix(0, nrow = length(markerframe[[i]]), ncol = Num_celltype)
        
        t_m[,1]=1
        
        }
    else{
        t_tmp=matrix(0, nrow = length(markerframe[[i]]), ncol = Num_celltype)
        t_tmp[,i]=1
        t_m=rbind(t_m,t_tmp)
        
        }

    }
    
   delta <- t_m*dt 
   rho <- t_m
      
}
############################################################      
print('Num of Marker')
print(length(markerList))
    
X = round(sp_counts)
n = nrow(X)
p = Num_HVG
p_m = length(markerList)
q = dim_embed

## STdeconvolve
cd <- t(X)
## remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)
## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
      
threshold <- 1
test_data <- t(as.matrix(corpus))[rowSums(t(as.matrix(corpus))) > threshold, ]
test_data <- test_data[, colSums(test_data) > threshold]

## categorize X_marker and others
idx = match(markerList, colnames(X2))
X_u = X2[,-idx]
X_m = X2[,idx]
    
#add marker genes
test_data=cbind(test_data,X_m[rownames(test_data),setdiff(colnames(X_m),intersect(colnames(X_m), colnames(test_data)))])

ldas <- fitLDA(test_data, Ks = Num_celltype)
    
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
result <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)

deconProp <- result$theta
deconGexp <- result$beta

## rho for STdeconvolve rho_STdeconvolve may have NA
markerexpr = as.matrix(t(result$beta[,match(markerList, colnames(result$beta))]))
rho_STdeconvolve = matrix(0, p_m, Num_celltype)
for (i in 1:p_m){
rho_STdeconvolve[i,which.max(markerexpr[i,])] = 1
}
########################################################################################################################      

## generate beta
gmm = Mclust(deconProp, G = Num_topic)
# adjustedRandIndex(gmm$classification, true_label)
beta = matrix(0, Num_topic, Num_celltype)
for (i in 1:Num_topic){
print(apply(deconProp[gmm$classification==i,], 2, mean))
tmp = apply(deconProp[gmm$classification==i,], 2, mean)
beta[i,] = tmp/sum(tmp)
}

## generate Z_int and W_int
fit_factor = factpc(as.matrix(logX[,-idx]), q)
Z_int = fit_factor$scores
W_int = fit_factor$loadings


## generate initial pi
pi_initial = matrix(0.000001, n, Num_topic)
for (i in 1:Num_topic){
pi_initial[which(gmm$classification == i),i]=0.999998
}

## generate initial gene expression mu_int from Z and beta
mu_int = calmu_int(Z_int, gmm$classification, beta,Num_celltype,q)


## dimension reduction
Lambda = matrix(1, p-p_m, 1)
alpha = matrix(0, n, 1)
epsilon2 = matrix(1, n, q)
pi0_initial = as.matrix(rep(1/Num_topic,Num_topic))  
    
Initial = CountDeconvolution(X_u, Z_int, W_int, K = Num_celltype, T = Num_topic, beta = beta, mu = t(mu_int), pi = pi_initial, pi0 = pi0_initial,
                               Lambda = Lambda, alpha = alpha, epsilon2 = epsilon2,
                               max_iter = 50, eps = 1e-5, verbose = TRUE,
                               fixed_loglambda_mu = FALSE, fixed_loglambda_s2 = FALSE,
                               fixed_z_mu = FALSE, fixed_z_s2 = FALSE,
                               fixed_m = TRUE, fixed_b = TRUE,
                               fixed_pi = TRUE,
                               fixed_W = FALSE, fixed_alpha = FALSE, fixed_Lambda = FALSE,
                               fixed_pi0 = TRUE, fixed_epsilon2 = TRUE,
                               fixed_mu = TRUE, fixed_sigma2 = TRUE, fixed_beta2 = TRUE, alpha_spot_specific = TRUE, cutoff = 1000000)
    
########################################################################################################################         
print('finish Initial')
    
index_A_to_B <- sapply(1:ncol(col_match(rho,rho_STdeconvolve)), function(i) {
  which(apply(rho, 2, identical, col_match(rho,rho_STdeconvolve)[, i]))
})
    
rho = col_match(rho,rho_STdeconvolve)
delta = col_match(delta,rho_STdeconvolve)
print('finish match')
## generate initial gene expression mu_int from Z and beta
mu_int = calmu_int(Initial$z_mu, gmm$classification, beta,Num_celltype,q)

Adj = SC.MEB::getneighborhood_fast(as.matrix(pos), cutoff = 1.3)
# summary(colSums(Adj))
y_int = apply(Initial$pi, 1, which.max)
epsilon2 = matrix(1, n, q)

beta0_y = matrix(0, Num_topic, 1)*(1/Num_topic)
beta_y_grid = matrix(seq(0,4,0.2),21 ,1)

Lambda_m <- rep(1, p_m)
m_int = matrix(0, p_m, Num_celltype)
print('finish init') 
## to obtain the loglambda_mu
fit <- DeconvMarker(X_m, rho, T = Num_topic, beta, delta, pi_initial, pi0_initial, Lambda_m, m_int, p, 1e-5, TRUE); 
#colnames(fit$TERM) <- c("cd", "pi_it", "nu_ijk", "m_jk", "delta_jk", "beta_tk", "pi_t", "nu_k", "Lambda_j")
print('finish fit')
########################################################################################################################   
    
mu_int = calmu_int(Initial$z_mu, gmm$classification, beta,Num_celltype,q)
m_int = calmu_int(fit$loglambda_mu, gmm$classification, beta,Num_celltype,p_m)
#iter 50 100 50
out2 = JointSTCountDecon(X_m, rho, X_u, y_int, Adj, Initial$z_mu, Initial$W, K = Num_celltype, T = Num_topic, beta = beta, mu = t(mu_int), pi = pi_initial, pi0 = pi0_initial,
                      Lambda = Initial$Lambda, alpha = Initial$alpha,epsilon2 = epsilon2,
                      beta0_y = beta0_y, beta_y_grid = beta_y_grid, Lambda_m = Lambda_m, delta = delta, m = m_int, 
                      max_iter = 50, maxIter_ICM = 50, eps = 1e-4, verbose = TRUE,
                      fixed_loglambda_mu = FALSE, fixed_loglambda_s2 = FALSE,
                      fixed_z_mu = FALSE, fixed_z_s2 = FALSE,
                      fixed_m = FALSE, fixed_b = FALSE,
                      fixed_pi = FALSE,
                      fixed_W = FALSE, fixed_alpha = FALSE, fixed_Lambda = FALSE,
                      fixed_pi0 = FALSE, fixed_epsilon2 = FALSE,
                      fixed_mu = FALSE, fixed_sigma2 = FALSE, fixed_beta2 = FALSE, fixed_beta_y = FALSE, alpha_spot_specific = TRUE, cutoff = 300)

time_spot1 <- Sys.time()
print('finish')
print(time_spot1-time_spot0)
    
# deconProp = matrix(0,n,Num_celltype)
# for (i in 1:n){
# idx = which.max(out2$pi[i,])
# deconProp[i,] = out2$beta[idx,]
# }
output_=list()
output_$pi=out2$pi
output_$beta=out2$beta
output_$order=index_A_to_B
# out2$STdecon=result 
output_$pos=pos
# out2$eta=fit$nu
########################################################################################################################   
    return (output_)
}
  

