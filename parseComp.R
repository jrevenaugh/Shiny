parseComp <- function( comp_str = "" ) {
  a <- unlist( strsplit( comp_str, "," ) )
  n <- length( a )
  if ( n == 0 ) return( NA )
  
  for ( i in 1:n ) {
    b <- as.numeric( unlist( strsplit( a[i], ":" ) ) )
    if ( length( b ) == 2 ) {
      if ( i == 1 ) 
        l <- list( b[1]:b[2] )
      else
        l[[i]] <- b[1]:b[2]
    }
    else {
      if ( i == 1 ) 
        l <- list( b[1] )
      else
        l[[i]] <- b[1]
    }
  }
  return( l )
}