dir = 'C:/Users/aga/Dropbox/MPHIL/VU/MitISEM/MyMit/include/Debug'
setwd(dir)

theta <-  read.csv("theta.csv", header = FALSE,  sep="," )
mu <-  read.csv("mu.csv", header = FALSE,  sep="," )
Sigma <-  read.csv("Sigma.csv", header = FALSE,  sep="," )
p <- read.csv("p.csv", header = FALSE,  sep="," )
df <- read.csv("df.csv", header = FALSE,  sep="," )

theta = matrix(unlist(theta),nrow=10000,ncol=2,byrow=FALSE)
mu = matrix(unlist(mu),nrow=2,ncol=2,byrow=FALSE)
Sigma = matrix(unlist(Sigma),nrow=2,ncol=4,byrow=FALSE)
df = as.vector(unlist(df))
p = as.vector(unlist(p))

mit = list(p=p, mu=mu, Sigma=Sigma, df=df)
isMit(mit)

d=dmvgt(theta=theta,mit=mit,log=TRUE)
write.csv(d, "d_R.csv",row.names=FALSE)


H=length(mit$p)
k = ncol(mit$mu)