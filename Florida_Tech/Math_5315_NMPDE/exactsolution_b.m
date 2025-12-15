
function fn = exactsolution(x)

fn = sin(pi*x)^40;

%if(x>0.4 && x < 0.6)
%  fn = 1.0; %(abs(x-0.25) < 0.1);
%else if(x <= 0.55 && x >= 0.05)
%   fn = cos(2*pi*(x-0.3))^2;
%else
%   fn = 0.0;
%end
end
