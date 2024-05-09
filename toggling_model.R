#Dowloading packages
library(dplyr) #For data wrangling
library(adana) #For manipulations between binary and decimal systems
library(ggplot2)
library(stringr) #For editing strings
library(deSolve) #For creating ODE model
library(scales)

#Defning functions


n <- 3

viral.escape <- function(column1, column2){#insert two vectors - first host, second virus
  B.comp <- numeric(length = length(column1)) #Create an empty vector to store the information for B stage virus
  for (x in 1:length(column1)) { #start a for loop for each object in the vector
    if (column1[x] == 1 & column2[x] == 0 | column1[x] == 1 & column2[x] == 1 | column1[x] == 0 & column2[x] == 1) ( #check whether virus escapes
      B.comp[x] <- 1 #produce escape in the virus in the B compartment if eligible, or maintains escape mutation if eligible
    )
  }
  return((B.comp)) #return B compartment
}

to.B.comp <- function(n){ #where n is a number of loci in the model
  B.matrix <- matrix(NA,    # Create empty data frame, host type in rows, virus type in columns, the matrix is filled out with viral type for each host after escape
                     nrow = 2^n, #Number of host and virus types is equal to 2^n
                     ncol = 2^n,
                     dimnames = list(c(0:(2^n-1)), c(0:(2^n-1)))) #hosts start from 0 not from 1 so adjusting for it here
  for (x in c(1:(2^n))) { #for each row
    for (y in c(1:(2^n))) { #for each column
      B.matrix[x,y] <- bin2int(viral.escape(int2bin(x-1,n), int2bin(y-1,n))) #For each of the host type predicts the viral type after escape, additional functions which translate a binary number to a decimal number are used
    }
  }
  return(B.matrix) #Return ready matrix
}

matrix <- to.B.comp(n)

generate.A.sum.for.B <- function(n, i, j) { #Function for the loci matrix, i is host, j is virus type in the B stage of infection
  if (j %in% unique(matrix[i+1,])) {#check whether i host type can produce a j type virus in the B stage of infection, aka check if a certain Bij individual actually exists
    vector1 <- c() #Empty vector to store data
    for (x in 1:ncol(matrix)) { #for each column
      if (matrix[i+1,x] == j) { #check whether the virus type j can be produced in the B stage in the host type i
        vector1 <- c(vector1, paste("A", as.character(i), ".", as.character(x-1), sep = "", collapse = "")) #If yes, store in the vector - function goes through each identified host type, upon finding a match (virus type in the B stage) it extracts the virus type at the A stage of infection
      }
    }   
  }
  return(vector1) #return vector with found A stage individuals
}

generate.A.B.C.sum.for.A <- function(n, j) { #function to generate sum of all A and C individuals with a j virus, matrix represents of the "escape matrixes"
  sum <- c() #Create empty object
  for (x in 1:nrow(matrix)) { #Run the code through each row of the matrix
    sum <- c(sum, paste("A", as.character(x-1),".",as.character(j), sep = "", collapse = "")) #Extract individual sum compartments - each host with a virus type j
  }
  for (i in 1:nrow(matrix)) {#start for loop for each row
    if (j %in% matrix[i,]) {#Check whether a j type virus is obtained from any i type hosts after escape
      sum <- c(sum, paste("B", as.character(i-1), ".", as.character(j), sep = "", collapse = ""))
    } #store the individual in the vector
  }
  sum <- c(sum, paste("C", as.character(j), ".", as.character(j), sep = "", collapse = ""))
  return(sum)
}

generate.B.sum.for.C <- function(n, i) { #Generate a sum to describe each compartment C where i is the C host type, matrix is the escape matrix generated with the to.B.comp function
  vector <- c() #Create an empty vector
  for (j in unique(matrix[i+1,])) { #Run a for loop - for each row (host type) find unique virus types and link them together in a sum
    vector <- c(vector, paste("B", as.character(i), ".", as.character(j), sep = "", collapse = ""))
  }
  return(vector)
}

generate.all.A <- function(n){ #Generates sum of all A compartments for n loci
  vector <- c() #Create an empty vector
  for (i in 0:(nrow(matrix)-1)) {#Extract each row
    for (j in 0:(nrow(matrix)-1)) {#Extract each column
      vector <- c(vector, paste("A", as.character(i), ".", as.character(j), sep = "", collapse = ""))
    }} #Paste the Aij individual inside the vector
  return(vector)
}

generate.all.B <- function(n){ #Generates sum of all B compartments for n loci
  vector <- c() #Create an empty vector
  for (i in 0:(nrow(matrix))-1) { #for each row (host type)
    for (j in unique(matrix[i+1,])) { #check which virus types are created after escape
      vector <- c(vector, paste("B", as.character(i), ".", as.character(j), sep = "", collapse = ""))
    }} #Paste the Bij individual
  return(vector)
}

generate.all.C <- function(n){ #Generates sum of all C compartements for n loci
  vector <- c() #Create an empty vector
  for (x in 0:(nrow(matrix)-1)){ #for each row (host type = virus type for C)
    vector <- c(vector, paste("C", as.character(x), ".", as.character(x), sep = "", collapse = ""))
  } #Return C (C is only in a form of Ci.j where i=j)
  return(vector)
}

generate.all.N <- function(n) { #n is number of loci
  vector <- c("S", generate.all.A(n), generate.all.B(n), generate.all.C(n))
  return(vector)
}

generate.A.matrix.row <- function(matrix, n){
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      type <- rownames(matrix)[i][1]
      virus_type <- sub(".*\\.", "", type)
      vector <- generate.A.B.C.sum.for.A(n, virus_type)
      if (colnames(matrix)[j] %in% vector) {
        substring <- sub("\\..*", "", type)
        host_type <- substr(substring, 2, nchar(substring))
        proportion <- paste("p", as.character(host_type), sep = "")
        px <- get(proportion, envir = parent.frame())
        sum2 <- beta * px
        matrix[i,j] <- sum2
      } else {
        matrix[i,j] <- 0
      }
    }
  }
  return(matrix)
}

generate.B.matrix.row <- function(matrix, n){
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      type <- rownames(matrix)[i][1]
      virus_type <- sub(".*\\.", "", type)
      host_type <- regmatches(type, regexpr("(?<=B)(.*?)(?=\\.)", type, perl=TRUE))[[1]]
      vector <- generate.A.sum.for.B(n, as.numeric(host_type), as.numeric(virus_type))
      if (rownames(matrix)[i] == colnames(matrix)[j]) {
        sum1 <- -(rho+D)
        matrix[i,j] <- sum1
      } else if (colnames(matrix)[j] %in% vector) {
        matrix[i,j] <- epsilon
      } else {
        matrix[i,j] <- 0
      }
    }
  }
  return(matrix)
}

generate.C.matrix.row <- function(matrix, n){
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      type <- rownames(matrix)[i][1]
      host_type <- sub(".*\\.", "", type)
      vector <- generate.B.sum.for.C(n, as.numeric(host_type))
      if (rownames(matrix)[i] == colnames(matrix)[j]) {
        matrix[i,j] <- -d -D
      } else if (colnames(matrix)[j] %in% vector) {
        matrix[i,j] <- rho
      } else {
        matrix[i,j] <- 0
      }
    }
  }
  return(matrix)
}

generate.BC.matrix <- function(n) {
  allA <- generate.all.A(n)
  allB <- generate.all.B(n)
  allC <- generate.all.C(n)
  allI <- c(allA, allB, allC)
  transition_matrix <- matrix(data = NA, nrow = length(allI), ncol = length(allI), dimnames = list(allI, allI))
  transition_matrix[1:(length(allA)),] <- 0
  matrixBB <- transition_matrix[(length(allA)+1):(length(allA)+length(allB)),]
  transition_matrix[(length(allA)+1):(length(allA)+length(allB)),] <- generate.B.matrix.row(matrixBB,n)
  matrixCC <- transition_matrix[(1+length(allA)+length(allB)):(nrow(transition_matrix)),]
  transition_matrix[(1+length(allA)+length(allB)):(nrow(transition_matrix)),] <- generate.C.matrix.row(matrixCC, n)
  return(transition_matrix)
}

generate.A.matrix <- function(n){
  allA <- generate.all.A(n)
  allB <- generate.all.B(n)
  allC <- generate.all.C(n)
  allI <- c(allA, allB, allC)
  transition_matrix <- matrix(data = NA, nrow = length(allI), ncol = length(allI), dimnames = list(allI, allI))
  matrixAA <- transition_matrix[1:(length(allA)),]
  transition_matrix[1:(length(allA)),] <- generate.A.matrix.row(matrixAA, n)
  transition_matrix[(length(allA)+1):(nrow(transition_matrix)),] <- 0
  return(transition_matrix)
}

generate.A.diagonal.matrix <- function(n){
  allA <- generate.all.A(n)
  allB <- generate.all.B(n)
  allC <- generate.all.C(n)
  allI <- c(allA, allB, allC)
  transition_matrix <- matrix(data = 0, nrow = length(allI), ncol = length(allI), dimnames = list(allI, allI))
  diag(transition_matrix[1:(length(allA)),]) <- - epsilon - D
  return(transition_matrix)
  
}

escape_revert_model <- function(n){
  transition_matrixA <- generate.A.matrix(n)
  transition_matrix_diagonalA <- generate.A.diagonal.matrix(n)
  transition_matrixBC <- generate.BC.matrix(n)
  SI <- function(time, state, params) {
    with(as.list(c(time, state, params)), {
      total_infected <- as.numeric(sum(state[2:length(state)]))
      total_population <- as.numeric(sum(state[1:length(state)]))
      S <- state[1]
      dS <- -beta*S*total_infected/total_population - D*S + Q
      Xij <- state[2:length(state)]
      dXij <- (S * (1/total_population) * transition_matrixA%*%Xij + transition_matrix_diagonalA%*%Xij+ transition_matrixBC%*%Xij)
      return(list(c(dS, dXij)))
    })
  }
  inits <- c(1000, 1, rep(0, times = ncol(transition_matrixA)-1))
  
  params = params1
  
  times = seq(0, generations, by = 1)
  
  rootfun <- function(time, state, params) {
    dstate <- unlist(SI(time, state, params))
    sum(abs(dstate)) - 1e-4
  }
  
  out.root.finalmodel <- as.data.frame(lsodar(func = SI, y = inits, parms = params,
                                              times = times, rootfun = rootfun))
  return(out.root.finalmodel)
}

return_data_frame <- function(n){
  loci_model <- escape_revert_model(n) 
  colnames(loci_model)[2:ncol(loci_model)] <- c(generate.all.N(n))
  return(loci_model)
}

generate.matched.equations <- function(n){#For how many loci
  vector <- c() #Create an empty vector
  for (i in 0:(nrow(matrix)-1)) { #for each row (aka the number of all host/virus types)
    vector <- c(vector, paste("A", as.character(i), ".", as.character(i), sep = ""), paste("B", as.character(i), ".", as.character(i), sep = ""), paste("C", as.character(i), ".", as.character(i), sep = ""))
  }
  return(vector)
}

count.escapes <- function(x, y, n) {#Where x and y are the virus type before (x) and after (y) escape and n is the number of loci
  table <- data.frame(int2bin(x,n), int2bin(y,n)) #Create a dataframe to compare two viral types
  colnames(table) <- c("Before", "After") #Naming the columns 
  summy <- sum(table$Before == 0 & table$After == 1) #CSumming how many loci escape
  return(summy)
}

group.i.escapes <- function(i, n) { #Where i is the number of escapes we want to investigate and n is the number of loci
  vector <- c() #Create empty vector to store data
  for (k in 1:nrow(matrix)) { #For each host type
    for (l in 1:ncol(matrix)) { #For each virus
      if (count.escapes(l-1,matrix[k,l],n)==i){ #Check whether the certain pair produces i number of mutations
        vector <- c(vector, paste("A", as.character(k-1), ".", as.character(l-1), sep = "")) #Save the A stage individual
      }
    }
    
  }
  return(vector)
}

count.reversions <- function(x, y, n) {#Where x and y are the virus type before (x) and after (y) escape and n is the number of loci
  table <- data.frame(int2bin(x,n), int2bin(y,n)) #Create a dataframe to compare two viral types
  colnames(table) <- c("Before", "After") #Naming the columns 
  summy <- sum(table$Before == 1 & table$After == 0) #CSumming how many loci escape
  return(summy)
}

group.i.reversions <- function(i, n) { #Number of reversion is i, and number of loci is n
  vector <- c() #Create empty vector
  for (k in 1:nrow(matrix)) { #For each row (host type)
    for (l in unique(matrix[k,])) { #For each virus type for the host type in B compartment
      if (count.reversions(l, k-1, n)==i){ #Check whether such individual produces i reversions
        vector <- c(vector, paste("B", as.character(k-1), ".", as.character(l), sep = "")) #If yes, save the individual
      }
      
    }
    
  }
  return(vector)
}

find.0.loci <- function(n, x) {
  vector <- c()
  for (i in 0:(2^n-1)) {
    binary <- int2bin(i, n)
    if (binary[x] == 0){
      vector <- c(vector, paste("rPV", as.character(i), sep = ""))
    }
  }
  vector
}

generate_data_for_cluster <- function(matrix, n){
  test_model_cluster <- dplyr::mutate(matrix, 
                                      N = rowSums(matrix[2:ncol(matrix)]), 
                                      I = rowSums(matrix[3:ncol(matrix)]),
                                      matched = rowSums(matrix[,generate.matched.equations(n), drop = FALSE]),
                                      mismatched = I - matched,
                                      A = rowSums(matrix[,generate.all.A(n), drop = FALSE]),
                                      B = rowSums(matrix[,generate.all.B(n), drop = FALSE]),
                                      C = rowSums(matrix[,generate.all.C(n), drop = FALSE]))
  for (i in 1:n) {
    column <- paste0("leave", as.character(i), "m", sep = "")
    test_model_cluster <- dplyr::mutate(test_model_cluster,
                                        !!column := i * (rowSums(test_model_cluster[,group.i.escapes(i,n), drop = FALSE]) - rowSums(test_model_cluster[,group.i.escapes(i,n), drop = FALSE]) * (exp(1)^(-epsilon*1)) + rowSums(test_model_cluster[,group.i.reversions(i,n), drop = FALSE]) - rowSums(test_model_cluster[,group.i.reversions(i,n), drop = FALSE]) * (exp(1)^(-rho*1))))
  }
  test_model_cluster <- dplyr::mutate(test_model_cluster,
                                      wh.er = rowSums(dplyr::select(test_model_cluster, (ncol(test_model_cluster) - n + 1):ncol(test_model_cluster)))/I)
  for (j in 0:(2^n-1)) {
    column <- paste0("rPV", as.character(j), sep = "")
    test_model_cluster <- dplyr::mutate(test_model_cluster,
                                        !!column := rowSums(test_model_cluster[,generate.A.B.C.sum.for.A(n,j), drop = FALSE])/I)
  }
  for (m in 1:n) {
    column <- paste0("rPL", as.character(m), sep = "")
    test_model_cluster <- dplyr::mutate(test_model_cluster,
                                        !!column := rowSums(test_model_cluster[,find.0.loci(n,m), drop = FALSE]))
  }
  dummy_df <- test_model_cluster
  for (k in 1:n) {
    column <- paste0("rPLdif", as.character(k), sep = "")
    test_model_cluster <- dplyr::mutate(test_model_cluster,
                                        !!column := abs((dplyr::select(dummy_df, (ncol(dummy_df) - n + k)))-(lag((dplyr::select(dummy_df, (ncol(dummy_df) - n + k)))))))
  }
  test_model_cluster <- dplyr::mutate(test_model_cluster,
                                      pb.er = rowSums(dplyr::select(test_model_cluster, (ncol(test_model_cluster) - n + 1):ncol(test_model_cluster))))
  
  return(test_model_cluster)
}

summary_statistics <- function(df) {
  df_last <- tail(df, n = 1)
  matched <- df_last$matched/df_last$I
  S <- df_last$S
  I <- df_last$I
  A <- df_last$A
  B <- df_last$B
  C <- df_last$C
  wh.er <- df_last$wh.er
  pb.er <- df_last$pb.er
  wh.erC <- df_last$wh.er + Constant
  pb.erC <- df_last$pb.er + Constant
  ratio <- wh.er/pb.er
  ratioC <- wh.erC/pb.erC
  generations <- nrow(df)
  Crequired <- wh.er/3 
  summary_vector <- c(n, matched,S,I,A,B,C,wh.er,pb.er,ratio,wh.erC,pb.erC,ratioC,generations,Crequired)
  names_vector <- c('n', 'matched', 'S', 'I', 'A', 'B', 'C', 'wh.er', 'pb.er', 'ratio', 'wh.erC', 'pb.erC', 'ratioC', 'generations', 'Crequired')
  summary_df <- data.frame(matrix(summary_vector, nrow = 1, dimnames = list(NULL, names_vector)))
  return(summary_df)
}

beta=0.3/365
epsilon=1/365 #Escape rate
rho=0.25/365 #Reversion rate
Q=50/365 #Birth Rate
d=0.2/365 #Infection-related death
D=0.02/365 #Natural Death Rate
generations=200*365
Constant = 0.00025



p_vector_values <- c()
for(i in 0:(2^n - 1)) {
  nam <- paste("p", i, sep = "")
  integeri <- int2bin(i, n)    
  sumi <- sum(integeri == 1)
  object_value <- 0.1^sumi*0.9^(n - sumi)
  p_vector_values <- c(p_vector_values, object_value)
  assign(nam, object_value)
}


pvector_names <- paste0("p", 0:(2^n - 1))
presult_vector <- setNames(p_vector_values, pvector_names)


params1 = c(
  Q = Q, #Birth
  D = D, #natural death
  d = d, #virus related death
  epsilon = epsilon, #Escape rate
  beta = beta, #Transmission Rate
  rho = rho,
  presult_vector#Reversion rate
  
)

loci_data <- return_data_frame(n)
loci_data_full <- generate_data_for_cluster(loci_data, n)


plot_NSI <- ggplot(loci_data_full, aes(time)) +
  geom_line(aes(y=N, colour = "N")) +
  geom_line(aes(y=I, colour = "I")) +
  geom_line(aes(y=S, colour = "S")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept = 30000, linetype="dashed", 
             color = "grey", size= 0.4) +
  geom_vline(xintercept = 12000, linetype="dashed", 
             color = "grey", size=0.4) +
  annotate(geom="text", x=6000, y=2500, label="Exponential Growth",
           color="black") +
  annotate(geom="text", x=21000, y=2500, label="Saturation",
           color="black") +
  annotate(geom="text", x=37000, y=2500, label="Epidemic \nEquilibrium",
           color="black") +
  xlab(label = "Time (days)") +
  ylab(label = "Number of \nindividuals") +
  labs(colour = "Infection \nStage") +
  annotation_logticks() +
  scale_color_manual(breaks = c("N", "I", "S"), labels = c("All hosts", "Infected", "Susceptible"), values = c("black", "red", "blue")) +
  theme_classic() +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(size =12),
        axis.text.y = element_text(size =12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14))


plot_NSI

plot_matched <- ggplot(loci_data_full, aes(time)) +
  geom_line(aes(y=matched/I, colour = "Matched")) +
  geom_line(aes(y=mismatched/I, colour = "Mismatched")) +
  geom_vline(xintercept = 30000, linetype="dashed", 
             color = "grey", size=1) +
  geom_vline(xintercept = 12000, linetype="dashed", 
             color = "grey", size=1) +
  annotate(geom="text", x=6000, y=1.1, label="Exponential Growth",
           color="black") +
  annotate(geom="text", x=21000, y=1.1, label="Saturation",
           color="black") +
  annotate(geom="text", x=37000, y=1.1, label="Epidemic \nEquilibrium",
           color="black") +
  labs(colour = "Host and virus \ncombination") +
  scale_color_manual(values = c("blue", "red")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  ylab(label = "Prevalence") +
  xlab(label = "Time (days)") +
  theme_classic() +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(size =12),
        axis.text.y = element_text(size =12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14))

plot_matched

plot_virustypes <- ggplot(loci_data_full, aes(time)) +
  geom_line(aes(y=rPV0, colour = "Virus 0")) + 
  geom_line(aes(y=rPV1, colour = "Virus 1")) + 
  geom_line(aes(y=rPV2, colour = "Virus 2")) + 
  geom_line(aes(y=rPV3, colour = "Virus 3")) + 
  geom_line(aes(y=rPV4, colour = "Virus 4")) + 
  geom_line(aes(y=rPV5, colour = "Virus 5")) + 
  geom_line(aes(y=rPV6, colour = "Virus 6")) + 
  geom_line(aes(y=rPV7, colour = "Virus 7")) + 
  labs(colour = "Virus type") +
  geom_vline(xintercept = 30000, linetype="dashed", 
             color = "grey", size=0.4) +
  geom_vline(xintercept = 12000, linetype="dashed", 
             color = "grey", size=0.4) +
  annotate(geom="text", x=6000, y=1.1, label="Exponential Growth",
           color="black") +
  annotate(geom="text", x=21000, y=1.1, label="Saturation",
           color="black") +
  annotate(geom="text", x=36000, y=1.1, label="Epidemic \nEquilibrium",
           color="black") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_viridis_d() +
  ylab(label = "Prevalence") +
  xlab(label = "Time (days)") +
  theme_classic() +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(size =12),
        axis.text.y = element_text(size =12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14))

plot_virustypes

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


plot_ers <- ggplot(loci_data_full, aes(time)) +
  geom_line(aes(y=wh.er, colour = "Within-host")) + 
  geom_line(aes(y=pb.er, colour = "Between-host")) + 
  labs(colour = "Evolutionary rate") +
  scale_color_manual(values = c("red", "blue")) +
  geom_vline(xintercept = 30000, linetype="dashed", 
             color = "grey", size=1) +
  geom_vline(xintercept = 12000, linetype="dashed", 
             color = "grey", size=1) +
  annotate(geom="text", x=6000, y=0.0004, label="Exponential \nGrowth",
           color="black") +
  annotate(geom="text", x=21000, y=0.0004, label="Saturation",
           color="black") +
  annotate(geom="text", x=36000, y=0.0004, label="Epidemic \nEquilibrium",
           color="black") +
  ylab(label = "Evolutionary Rate \n(substitutions/genome/day)") +
  xlab(label = "Time (days)") +
  scale_y_continuous(label=scientific_10) +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=0.1) +
  theme_classic() +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(size =12),
        axis.text.y = element_text(size =12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14))


plot_ers


plot_bhers <- ggplot(loci_data_full, aes(time)) +
  geom_line(aes(y=pb.er, colour = "Between-host")) + 
  labs(colour = "Evolutionary rate") +
  scale_y_continuous(label=scientific_10, limits = c(NA, 0.0000011)) +
  ylab(label = "Evolutionary Rate") +
  xlab(label = "Time") +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "grey", size=0.1) +
  scale_x_continuous(breaks = seq(30000, 46000, by = 2000), limits = c(30000, NA)) +
  theme_classic()

plot_bhers


plot_ratio <- ggplot(loci_data_full, aes(time)) +
  geom_line(aes(y=wh.er/pb.er, colour = "Mis-match Ratio")) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(label = "Time (days)") +
  ylab(label = "Mis-match ratio") +
  geom_hline(yintercept = 2, linetype="dashed", 
             color = "grey", size=0.25) +
  geom_hline(yintercept = 6, linetype="dashed", 
             color = "grey", size=0.25) +
  annotate(geom="text", x=37500, y=6, label="6x",
           color="dark grey") +
  annotate(geom="text", x=37500, y=2, label="2x",
           color="dark grey") +
  theme_classic() +
  xlim(10,38304) +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text.x = element_text(size =12),
        axis.text.y = element_text(size =12),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

plot_ratio


write.csv(summary_statistics(loci_data_full), paste0("outputs/summary.vector.n.", as.character(n), ".csv"))


