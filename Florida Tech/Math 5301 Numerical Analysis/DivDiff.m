function div = DivDiff(x,y)
difference = y;
for j=2:length(y)
    temp = difference;
    for i=j:length(y)
        difference(i) = (temp(i)-temp(i-1))/(x(i)-x(i-(j-1)));
    end
end
div = difference;
end
