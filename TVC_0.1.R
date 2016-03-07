library(tsDyn)
library(data.table)
data("barry")

# R port of the gauss code by Luis Filipe Martins
y1 <- log(barry[,"dolcan"], base = exp(1))	
y2 <- log(barry[,"cpiUSA"], base = exp(1))	
y3 <- log(barry[,"cpiCAN"], base = exp(1))


y <- cbind(y1, y2, y3)

y <- as.data.table(y)

n <- y[,.N]
k <- ncol(y)

# Function to create lagged data.tables
lag.data.table <- function(data, max.lag = 1){
        data[, shift(data, max.lag, give.names = TRUE)]
}
# Output conformable to the the output of varlags
cbind(y[-(1:lags),] ,na.exclude(lag.data.table(y, max.lag = 1)))

# Function diffmatrix to calculate the first differences
diffmatrix <- function(x,max.diff = 1,max.lag = 1){
        #Add if condition to make it possible to differentiate between matrix and vector                  
        if(is.vector(x) == TRUE ){
                myx <- embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
                colnames(myx) <- paste("v1.d",max.diff, sep=".")
                return(myx)
        }
        
        else if(is.matrix(x) == TRUE){
                myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
                mycolnames <- colnames(x)
                colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
                return(myx)
        }
        #if matrix the result should be 0, if only a vector it should be 1
        else if(as.integer(is.null(ncol(x))) == 0 ){
                myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
                mycolnames <- colnames(x)
                colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
                return(myx)
                
        }
}

as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1))
# Removal of na.s
na.exclude(as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1)))


# In gauss a preceding . means that the following calculations are meant to be carried out element-wise
#Function to calculate the time-varying cointegration
tvcoint <- function(y,p,m){
#local ystar_1,ysub,y_1,dy,dylags,T,betau,resu,betav,resv,S00,S01,S10,S11,S00inv,S11inv,A,eval,evec,evs;
                ystar_1 <- ycheb(y,p,m);		
                
                dy <- na.exclude(as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1)))
                dylags <-  na.exclude(lag.data.table(dy, max.lag = p))
                dy <- dy[-(1:lags)]
                myT <- nrow(dylags);             
                dylags <- cbind(rep(1,myT), dylags)  
                betau <- dy/dylags;          
                resu <- dy-dylags*betau;     
                betav <- ystar_1/dylags;     
                resv <- ystar_1-dylags*betav;
                S00 <- (crossprod(resu))/ myT;        
                S01 <- (crossprod(resu,resv))/myT;        
                S10 <- t(S01);                
                S11 <- (crossprod(resv))/myT;       
                S00inv <- solve(S00);       
                S11inv <- solve(S11);       
                A <- S11inv*S10*S00inv*S01; 
                list(eval,evec) <- eigen(A);     
                evs <- rev(sortc(eval~t(evec),1)); 
                evec <- t(evs[.,2:ncol(evs)])	
                return(evs[.,1],evec,det(S00))
}

# function to calculate the chebycheff polynomials              
ycheb <- function(mydata,varord,chebdim){
            #  local i,yst,k,n,nn,ind;
              k <- ncol(mydata);
              nn <- nrow(mydata)
              yst <- matrix(nrow = nn-varord-1, ncol = (chebdim+1)*k)
              yst[, 1:k] <- as.matrix(mydata[(varord+1):(nn-1), ])
              if(chebdim == 0){
              return(yst);
                } else {
              n <- length(yst[, 1]);
              ind <- seq(from = varord+2, by = 1, length.out = n);
              i = 1;
              while(i <= chebdim){
              yst[, (i*k+1):((i+1)*k)] <- sqrt(2)*cos(i*pi*(ind-0.5)/n)*yst[,1:k]
              i = i+1
              }
            #  return(yst);
                }
        }
# calculate the lagged time series              
varlags <- function(var,lags){
              xlags = t(shiftrt(t(matrix(data = 1, nrow = 1, ncol = lags)) * y,
                              seq(from = 1-lags,by = 1, length.out = lags) * matrix(data = 1, nrow = ncol(y), ncol = 1),
                           miss(0,0)));
              
              # PURPOSE Trims rows from the top and/or bottom of a matrix.
              # FORMAT y = trimr(x,t,b);
              # INPUT x N×K matrix from which rows are to be trimmed.
              # t scalar containing the number of rows which are to be removed from
              # the top of x.
              # b scalar containing the number of rows which are to be removed from
              # the bottom of x.    
        return(trimr(var,lags,0),trimr(xlags,0,lags));
             #First returned data matrix
              var[-(1:lags),]
             #Second returned data matrix
              na.exclude(lag.data.table(y, max.lag = 1))
}



mmax <- round(n/10)	# maximum dimension of Chebishev Polynomials
p <- 1;			# VAR order

#Setting up the matrices to hold the values

lrtvc <- matrix(nrow = mmax, ncol = k);	     # TVC Stat (row) m=1,...,mmax and (col) r=1,..,k
lrtvcpv <- matrix(nrow = mmax, ncol = k);	     # P-values using the Asymp Distr (chisquare(mrk))
betat <- matrix(nrow = n-p-1, ncol = k*mmax)  # b_t for r=1.  (row) i= observation; Cols 1 up to k= b1~...~bk (m=1); Cols k+1 up to 2k=b1~...~bk (m=2); ...
lnlikm <- matrix(nrow = mmax, ncol = k);	     # log likelihood for different m and r 
aic <- matrix(nrow = mmax, ncol = k);	     # AIC Akaike model selection criterion 
bic <- matrix(nrow = mmax, ncol = k);	     # BIC Schwarz msc 
hann <- matrix(nrow = mmax, ncol = k);	     # Hannan-Quinn msc 



list(ev0,evec0,det0)=tvcoint(y,p,0); #    /* standard cointegration */
ll0=ln(1-ev0);
#Function call to all previously defined functions
m=1;
while(m<=mmax){

        list(eval,evect,detm) <- tvcoint(y,p,m); #/* OUT: Lambda; Eigenvectors q1...qr...q(m+1)k; det(S00) */
                llm <- ln(1-eval);	
             #   /* k=3 and r=1 */		
                        ind <- seq(from = p+2, by = 1, length.out = n-p-1); 
                        beta1sum <- matrix(nrow = m, ncol = nrow(ind));	
                        beta2sum <- matrix(nrow = m, ncol = nrow(ind));	
                        beta3sum <- matrix(nrow = m, ncol = nrow(ind));	
                        
                        mm=1;
                        while(mm <= m){
                        beta1sum[mm,] <- evect[k*mm+1,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
			beta2sum[mm,] <- evect[k*mm+2,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
                        beta3sum[mm,] <- evect[k*mm+3,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
                                mm=mm+1;
                                    }
         betat[.,k*(m-1)+1:k*m] <- t(evect[1:k,1])+cbind(colSums(beta1sum),colSums(beta2sum),colSums(beta3sum));		
                                                  
                        r=1;
                        while(r<=k){
                        lrtvc[m,r] <- (n-p-1)*colSums(ll0[1:r, ]) - (n-p-1)*colSums(llm[1:r,.]);	
                        lrtvcpv[m,r] <- cdfchic(lrtvc[m,r],m*r*k);		
                        lnlikm[m,r] <- (ln(r)-k-k*ln(2*pi))*(n-p-1)/2 - colSums(llm[1:r,.])*(n-p-1)/2 - (ln(detm))*(n-p-1)/2;				npar = (m+1)*k*r+r*k+k^2+(k+(p-1)*k^2);		
                        aic[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*2/(n-p-1);					
                        bic[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*(ln(n-p-1))/(n-p-1);			
                        hann[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*(ln(ln(n-p-1)))*2/(n-p-1);	
                        r <- r+1;
                                  }
                        m <- m+1;
                }




save.image("myTVPData.RData", compress = "xz")
