cov_mat <- matrix(c(1, 0.8, 0.8, 1), byrow = T, ncol=2)
cov_mat
cov_mat2 <- cov_mat
diag(cov_mat2) <- diag(cov_mat2)*2
cov_mat3 <- cov_mat %*% matrix(c(2, 1, 1, 2), byrow = T, ncol=2)

n <- 6000
r1 <- mvtnorm::rmvnorm(n = n, mean = c(0,0), sigma = cov_mat)
r2 <- mvtnorm::rmvnorm(n = n, mean = c(0,0), sigma = cov_mat2)
r3 <- mvtnorm::rmvnorm(n = n, mean = c(0,0), sigma = cov_mat3)
par(mfrow=c(1,3))
plot(r1,
     xlim=c(-4, 4), ylim=c(-4, 4))
plot(r2,
     xlim=c(-4, 4), ylim=c(-4, 4))
plot(r3,
     xlim=c(-4, 4), ylim=c(-4, 4))


cor(r1[,1], r1[,2]) ## really different correlation (factor 2)
cor(r2[,1], r2[,2]) ## really different correlation (factor 2)
cor(r3[,1], r3[,2])


cov(r1[,1], r1[,2]) ## shared covariance
cov(r2[,1], r2[,2]) ## shared covariance
cov(r3[,1], r3[,2])
