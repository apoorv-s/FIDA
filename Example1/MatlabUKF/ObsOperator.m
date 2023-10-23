function skPos = ObsOperator(xk)
global time

lmd = double(xk(1));
uL = 2.0;
uR = 1.0;
xR = double(xk(2));

beta = 2.0/lmd;
alpha = (uL-uR)/(xR);
tBreak = 1/(beta*alpha);

if time < tBreak
    skPos = (beta*uL*time + xR + beta*uR*time)/2;
else
    sigma = (uL + uR)/lmd;
    skPos = sigma*time + (xR/2);
end

end