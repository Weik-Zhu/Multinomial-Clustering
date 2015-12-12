
H<- matrix(readBin("histograms.bin","double",640000),40000,16)

MultinomialEM<- function(H,k,tau){
  c<- rep(1,k)/k
  #initial centroids
  random<- sample(1:nrow(H),k)
  centroids<- list()
  for (i in 1:k){
    centroids[[i]]<- (H[random[i],]+1)/11
  }
  A_new<-E_step(H,c,centroids,k)
  d<-norm(A_new,"O")
  
  while(d>tau){
    print(d)
    A_old<-A_new
    list<-M_step(A_new,H,k)
    c<-list[[1]]
    centroids<-list[[2]]
    A_new<-E_step(H,c,centroids,k)
    A<-A_new-A_old
    d<-norm(A,"O")
  }
  m<-c()
  print("yes successful!!!!")
  for (i in 1:40000){
    m[i]<-which.max(A_new[i,])
  }
  return(m)
}

#E_step
E_step<-function(H,c,centroids,k){
  phy<- matrix(nrow =40000, ncol = k)
  for (i in 1:40000){
    for (l in 1:k){
      count<- 0
      for (j in 1:16){
        count<- count+H[i,j]*log2(centroids[[l]][j])
      }
      phy[i,l]<- exp(count)
    }
  }
  A<- matrix(nrow =40000, ncol = k)
  for (i in 1:40000){
    denominator<- sum(c*phy[i,])
    for(l in 1:k){
      A[i,l]<-c[l]*phy[i,l]/denominator
    }
  }
  return(A)
}
#M_step
M_step<-function(A,H,k){
  c<-c()
  b<-list()
  centroids<-list()
  
  for (i in 1:k){
    c[i]<- sum(A[,i])/40000
    count<- 0
    for(j in 1:40000){
      count<- count + A[j,i]*H[j,]
    }
    b[[i]]<- count
  }
  
  for (i in 1:k){
    centroids[[i]]<- b[[i]]/sum(b[[i]])
  }
  
  list<-list()
  list[[1]]<-c
  list[[2]]<-centroids
  return(list)
}

#get m
m<- MultinomialEM(H,k,tau)

#get image
mm<-matrix(m,ncol = 200)
z<-matrix(ncol = 200,nrow = 200)

for (i in 1:200){
  c<-mm[i,]
  d<-c()
  for (j in 1:200){
    d[j]<-c[201-j]
  }
  z[i,]<-d
}

x<-1:200
y<-1:200
image(x,y,z)