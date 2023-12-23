function Fval = RHSWB_phi1_perturb(t,T_mid, y, num, W, W1, gsyn, taus, Ip, Tp)
y = y(:);
Fval = zeros(size(y));

index = 1:num;
v = y(index); h = y(num+index); n = y(2*num+index); s = y(3*num+index);
ena=55.0; ek=-90.0; el=-65.0; esyn=-75.0;
gna=35.0; gk=9.0; gl=0.1; Iapp = 0.4; 
phi=1.0; 
ae0=4.0;
       

am = @(v)(-0.1*(v+35.0)./(exp(-0.1*(v+35.0))-1.0));
bm = @(v)(4.0*exp(-(v+60.0)/18.0));
minf = @(v)am(v)./(am(v)+bm(v));

ah = @(v)(0.07*exp(-(v+58.0)/20.0));
bh = @(v)(1.0./(exp(-0.1*(v+28.0))+1.0));
       
an = @(v)(-0.01*(v+34.0)./(exp(-0.1*(v+34.0))-1.0));
bn = @(v)(0.125*exp(-(v+44.0)/80.0));
 
ae = @(v)(ae0./(1+exp(-v./5)));

if t<=T_mid
    Fval(index) =  -gna*h.*(v-ena).*minf(v).^3 - gk*(v-ek).*n.^4 - gl*(v-el) ...
        + Iapp - gsyn*(W*s).*(v-esyn) + Ip*heaviside(t-Tp(1))*heaviside(Tp(2)-t);
else
    Fval(index) =  -gna*h.*(v-ena).*minf(v).^3 - gk*(v-ek).*n.^4 - gl*(v-el)...
        + Iapp - gsyn*(W1*s).*(v-esyn) + Ip*heaviside(t-Tp(1))*heaviside(Tp(2)-t);
end

Fval(num+index) = phi*(ah(v).*(1-h)-bh(v).*h); 
Fval(2*num+index) =   phi*(an(v).*(1-n)-bn(v).*n);
Fval(3*num+index) = -s/taus + ae(v).*(1-s);