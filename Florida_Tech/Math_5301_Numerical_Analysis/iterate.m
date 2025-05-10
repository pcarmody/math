ext = 0;
for i=1:n-1
	m = (pdata(i)+pdata(i+1))/2;
	est = est + inc/6*(pdata(i) + 4*m + pdata(i+1));
end
