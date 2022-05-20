using LinearAlgebra, Polynomials, SpecialPolynomials
using QuadGK, Calculus, CalculusWithJulia

Max = 25
Q = 0.7
z = 1
a1 = 1
a2 = 1 
a3 = 1
a4 = 1
alpha = 1
rc = 1
L = 1

function Zl(r)
   1 + (z-1)*exp(-a1*r)- r*(a3+a4*r)*exp(-a2*r)
end

function potential(r)
   -Zl(r)/r-alpha/(2*r^4)*(1-exp(-(r/rc)^6))
end

function potentialH(r)
   -1/r
end

function angular(l,r)
   l*(l+1)/(2*r^2)
end

function Sturmian(k,l,r)
   ö = zeros(1,Max)
   ö[k] = 1
   sqrt(factorial(big(k))/(2*factorial(big(2*l+1+k))))*exp(-Q*r)*(2*Q*r)^(l+1)*Laguerre{2*l+1}(vec(ö))(2*r*Q)
end

function DSturmian(k,l,r)
   if k <= 1
      return 0
   else
      ö = zeros(1,Max)
      ö[k] = 1
      öö = zeros(1,Max)
      öö[k-1] = 1
      return sqrt(factorial(big(k))/(2*factorial(big(2*l+1+k))))*Q*(-2*exp(-Q*r)*(2*Q*r)^(l+1)*Laguerre{2*l+2}(vec(öö))(2*r*Q)+2*(1+l)*exp(-Q*r)*(2*Q*r)^(l)*Laguerre{2*l+1}(vec(ö))(2*r*Q)-exp(-Q*r)*(2*Q*r)^(l+1)*Laguerre{2*l+1}(vec(ö))(2*r*Q))
   end
end

O = zeros((Max,Max))
Hamil1 = zeros((Max,Max))
Hamil2 = zeros((Max,Max))

for k1 = 1:Max
   for k2 = k1:Max
      O[k1,k2] = quadgk(r-> Sturmian(k1,L,r)*Sturmian(k2,L,r),0,100)[1]
      O[k2,k1] = O[k1,k2]
      Hamil1[k1,k2] = quadgk(r-> Sturmian(k1,L,r)*Sturmian(k2,L,r)*(potentialH(r)+angular(1,r)),0,100)[1]
      Hamil1[k2,k1] = Hamil1[k1,k2]
      Hamil2[k1,k2] = quadgk(r-> DSturmian(k1,L,r)*DSturmian(k2,L,r)/2,0,100)[1]
      Hamil2[k2,k1] = Hamil2[k1,k2]
      print(k1,k2)
   end
end

Hamil = Hamil1 + Hamil2
print(eigvals(Hamil,O))
