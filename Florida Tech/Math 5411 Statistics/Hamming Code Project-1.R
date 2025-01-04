# Hamming Code Project
##########################################################################################################

#Consider the problem 10 in Chapter 2, now, by adding Hamming code of size 4xm, up to m-bit changes can
#be corrected. For example, for Hamming code of size 4, 4 bits will be added, and single-bit change can be
#corrected, using a Hamming code of size 8, up to 2-bit changes (1-bit and 2-bit) can be corrected, and so
#on. Find the number of Humming code bits to be added to an 8-bit word for maximizing the probability
#that an 8-bit word is communicated correctly.


x = seq(8,8*20, by=1) 
rec = vector()
cord = vector()
max_prob = 0  
optimal_k = 0


for (var in x) {
  bit = 8+var
  #print(bit)
  corec = floor(bit/8)-1
  #print(corec)
  cord = append(cord, corec)
  
  p = 0.05

  cdf = pbinom(corec,bit,p)
  rec = append(rec, cdf) 
  #print(cdf)

  if (cdf > max_prob) {
    max_prob = cdf
    optimal_k = var  # Update optimal number of Hamming bits
  }
}

print(max_prob)
print(optimal_k)
plot(x, rec)
plot(x, rec/x)
plot(cord, rec)


########################################################################################
########################################################################################

x = seq(8,8*20, by=1) 
rec = vector()
cord = vector()
max_prob = 0  
optimal_k = 0

for (var in x) {
  bit = 8+var
  #print(bit)
  corec = floor(bit/8)-1
  #print(corec)
  cord = append(cord, corec)
  
  p = 0.05
  
  pmf = dbinom(x-8,bit,p)
  print
  
  i = 1
  prob = 0
  for (i in 1:(corec + 1)) {
    
    prob = prob + pmf[i]
    

  }
  #print(prob)
  rec = append(rec, prob)
  if (prob > max_prob) {
    max_prob = prob
    optimal_k = var  # Update optimal number of Hamming bits
  }
}
print(max_prob)
print(optimal_k)
plot(x, rec)
plot(x, rec/x)
plot(cord, rec)