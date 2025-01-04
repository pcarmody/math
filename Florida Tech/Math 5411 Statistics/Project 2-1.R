###########################################
# Poll Values from November 4th 2024

Xt = 875*0.516 # Number of People who voted for Trump
Xh = 875*0.465 # Number of People who voted for Harris
n = 875*0.516 + 875*0.465 # Number of people who voted for Trump or Harris

#ph = X/n

a = 0.02
z.a = qnorm(1-a)
z.a2 = qnorm(1-a/2)

###########################################
#Put a 98% upper bound Confidence Interval for the true proportion of the votes for 
#candidate with the lower proportion of the votes.

phh = Xh/n

psh = (phh+(z.a^2)/(2*n))/(1+(z.a^2)/n)
sigh = z.a*(sqrt(phh*(1-phh)/n + (z.a^2)/(4*n^2)))/(1+(z.a^2)/n)

pu = psh + sigh

print(paste("The Upper Bound CI for Harris is", pu))

###########################################
#Put a 98% lower bound Confidence Interval for the true proportion of the votes for the 
#candidate with the higher proportion of the votes.

pht = Xt/n

pst = (pht+(z.a^2)/(2*n))/(1+(z.a^2)/n)
sigt = z.a*(sqrt(pht*(1-pht)/n + (z.a^2)/(4*n^2)))/(1+(z.a^2)/n)

pl = pst - sigt

print(paste("The Lower Bound CI for Trump is", pl))


###########################################
#Put a 98% two-sided Confidence Interval on the difference of true proportions. Use the
#higher proportion minus the lower proportion as the estimated difference.

phd = (Xt-Xh)/n

#print(phd)

psd = (phd+(z.a2^2)/(2*n))/(1+(z.a2^2)/n)
sigd = z.a2*(sqrt(phd*(1-phd)/n + (z.a2^2)/(4*n^2)))/(1+(z.a2^2)/n)

pld = psd - sigd
pud = psd + sigd

print(paste("The Lower Bound CI for the difference is", pld, "The Upper Bound CI for the difference is", pud ))

###########################################
# True Values from the Election
vt = 1770242 # Number of votes Trump received
vh = 1582860 # Number of votes Harris received
vpr= vt + vh # Number of votes for Trump and Harris

pt = vt/vpr # Percentage of people who voted for Trump
ph = vh/vpr # Percentage of people who voted for Harris

dp = (vt-vh)/vpr # Difference in the percentage of people voting for each candidate


print(paste("The percentage of people who voted for Trump is", pt))
print(paste("The percentage of people who voted for Harris is", ph))
print(paste("The difference in the percentages is", dp))
###########################################
#Find the true error of the estimate using post-election proportions.


errt = abs(pt-pst) # Error in the proportions for Trump
errh = abs(ph-psh) # Error in the proportions for Harris

print(paste("The error in the proportions for Trump:", errt))
print(paste("The error in the proportions for Harris:", errh))

print("Does the difference in the true percentages fall in the CI predicted?")

if (pld<dp & dp<pud){
  print("Yes")
}else{
  print("No")
}