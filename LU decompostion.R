
##  LU DECOMPOSITION (Roll: bs2016, B.Stat(Hons.) 1st year)


LU = function(A, b){
  
  A1 = A; # creates a duplicate copy of the original matrix
  
  # LU decomposition only exists for invertible matrices
  stopifnot(ncol(A) == nrow(A));
  # the process stops here if A is not a square matrix
  
  n = nrow(A)
  # we don't have to store the number of columns separately because
  # we have already ensured that we will only proceed if A is a square matrix
  
  # LU decomposition will only exist if A is non-singular
  # so A must have full row and column rank
  stopifnot(qr(A)$rank == n);
  
  # qr() is a built-in function which performs the QR decomposition of a matrix
  
  
  L = function(j) { # a function that builds a lower triangular matrix
    tot = rep(0,n-j+1)
    if(j>1) {
      for(k in 1:(j-1)) 
        tot = tot + A[j:n,k] * A[k, j]
    }
    A[j:n,j] <<- A[j:n,j] - tot
    A
  }
  
  U = function(i) { # a function that an upper triangular matrix
    tot = rep(0,n-i)
    if(i>1) {
      for(k in 1:(i-1)) 
        tot = tot + A[i,k] * A[k, (i+1):n]
    }
    A[i, (i+1):n] <<- (A[i, (i+1):n] - tot)/A[i, i]
    A
  }
  # the functions L and U are same as those described in the class web page
  
  
  # n = nrow(A1);
  print("Matrix calculated at each step of Efficient Crout's Decomposition: ")
  for(i in 1 : n)
  {
    print(L(i));
    
    if(i != n){
      print(U(i));
    }
  }
  
  # we find the decomposition of A into LU
  # n = ncol(A1);
  L1 = matrix(0, n, n); 
  U1 = matrix(0, n, n);
  
  for(i in 1 : n){
    for(j in 1 : n){
      
      #assigning the diagonal elements
      if(i == j){
        
        L1[i, j] = A[i, j];
        U1[i, j] = 1; #as diagonal elements of U in Crout's Factorization 
      }
      
      #assigning the elements above the main diagonal
      if(i < j){
        
        L1[i, j] = 0; #L is lower triangular 
        U1[i, j] = A[i, j];
      }
      
      #assigning the elements below the main diagonal
      if(i > j){
        
        L1[i, j] = A[i, j];
        U1[i, j] = 0; #U is upper triangular
      }
    }
  }
  print("Original matrix A : ");
  print(A1);
  print("The lower triangular matrix L :");
  print(L1);
  print("The upper triangular matrix U : ");
  print(U1);
  
  # 1. Let y = Ux, then
  # 2. Solving Ax=B, is equivalent to solving Lx = y and Ux = y
  # 3. y=Ux can be solved by forward substitution
  # 4. Ly=b can be solved by backwarrd substitution
  
  
  # the two systems are consisitent simultaneously only if 
  # L and U are both non-singular (they have full column and row rank)

  rL = qr(L1)$rank;
  rU = qr(U1)$rank;
  if( rL != n | rU != n){
    print("Error. L or U is singular.")
  }
  else{
    # solving Ly = b
    y = rep(0, n);
      
      for(i in 1 : n){
        if(i > 1){
          
          # from 1 to i-1, all the unknowns are being calculated  
          for(j in 1 : (i-1)){
            
            b[i] = b[i] - y[j]*L1[i, j]; # all the columns have been traversed
          }
        }
        y[i] = b[i]/L1[i, i]; # we are solving this by forward substitution
      }
      
    # solving y = Ux  
    x = rep(0, n);
      for(i in seq(n, 1, by = -1)){
        
        if(i != n){
          
          # from i+1 to n, all the unknowns are being calculated  
          for(j in (i+1) : n){
            
            y[i] = y[i] - x[j]*U1[i, j]; # similar algorithm as before
          }
        }
        x[i] = y[i]/U1[i, i]; # by backward substitution
        }
    print("The solution of the system Ax = b : ");
    print(x);
  }
}