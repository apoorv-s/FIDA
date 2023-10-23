function [currState, jumpPoints] = IBSolverIC1(icPars, time)
%IBSolver: Takes in lambda (originally 2, denominator in forcing
% function in Inviscid burger. Takes in the parameters for the initial
% condition and solved for the u values at time t.

if iscell(icPars(1))
    alpha = double(cell2mat(icPars(1)));
    beta = double(cell2mat(icPars(2)));
    gamma = double(cell2mat(icPars(3)));
    a = double(cell2mat(icPars(4)));
    b = double(cell2mat(icPars(5)));
    c = double(cell2mat(icPars(6)));
    d = double(cell2mat(icPars(7)));
    lambda = double(cell2mat(icPars(8)));
else
    alpha = icPars(1);
    beta = icPars(2);
    gamma = icPars(3);
    a = icPars(4);
    b = icPars(5);
    c = icPars(6);
    d = icPars(7);
    lambda = icPars(8);
end

F=@(u) double(u.^2/lambda);
Fp=@(u) double(2*u/lambda);
% L = @(u) double(u.^2/lambda);

icDomain = [a, b, c, d];

tmpVar = chebfun('tmpVar', icDomain);
g = chebfun({sin(alpha*tmpVar), beta, sin(gamma*tmpVar)}, ...
                    icDomain, 'splitting', 'on');

G = cumsum(g);
d = minandmax(g)';
a = g.ends(abs(jump(g, g.ends)) > ...
                   10*vscale(g)*eps);

u=@(x)IBsolverHelper(Fp,F,g,G,d,a,x,time);
currState = chebfun(@(x) u(x), domain(g), ...
                    'splitting', 'on');
jumpPoints = currState.ends;
jumpPoints = jumpPoints(2:length(jumpPoints) - 1);

function u = IBsolverHelper(Fp, F, g, G, d, a, x, t)%
    cFp = chebfun(Fp, d, 'splitting', 'on'); 
    q=[];
    for k = 1:length(a)
    hk = a(k) - x + t*cFp;
        q = [q, roots(hk)'];
    end
    dom = unique([d, q]);
    h= chebfun(@(p) g(x - t*Fp(p)) - p, dom, 'splitting', 'on'); 
    p= [roots(h); q(:)];
    
    J = @(p) (p.*Fp(p) - F(p)).*t + G(x-Fp(p).*t);
%     J = @(p) L(-Fp(p))*t + G(x-Fp(p)*t); 
    [~, idx] = min(J(p));
    %
    u = p(idx);
end

end