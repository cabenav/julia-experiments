using LinearAlgebra, Polynomials, SpecialPolynomials
using QuadGK, Calculus, CalculusWithJulia
using GSL
using Plots 

Max =10
Q = 0.7
z = 1
a1 = 1
a2 = 1 
a3 = 1
a4 = 1
alpha = 1
rc = 1
L = 0

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
   sqrt(factorial(k)/(2*factorial(2*l+1+k)))*exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k,2*l+1,2*r*Q)
end

function DSturmian(k,l,r)
   if k < 1
      return  sqrt(factorial(k)/(2*factorial(2*l+1+k)))*Q*(2*(1+l)*exp(-Q*r)*(2*Q*r)^(l)*sf_laguerre_n(k,2*l+1,2*r*Q)-exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k,2*l+1,2*r*Q))
   else
      return  sqrt(factorial(k)/(2*factorial(2*l+1+k)))*Q*(-2*exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k-1,2*l+2,2*r*Q)+2*(1+l)*exp(-Q*r)*(2*Q*r)^(l)*sf_laguerre_n(k,2*l+1,2*r*Q)-exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k,2*l+1,2*r*Q))
   end
end

function Kinetic(k1,k2,L,r)
   DSturmian(k1,L,r^2)*DSturmian(k2,L,r^2)/2*2*r
end

function energy(k1,k2,L,r)
   Sturmian(k1,L,r^2)*Sturmian(k2,L,r^2)*(potentialH(r^2)+angular(L,r^2))*2*r
end

function Overlap(k1,k2,L,r)
   Sturmian(k1,L,r^2)*Sturmian(k2,L,r^2)*2*r
end

function loop(Max,L)
   O = zeros((Max+1,Max+1))
   Hamil1 = zeros((Max+1,Max+1))
   Hamil2 = zeros((Max+1,Max+1)) 
   for k1 = 0:Max
      for k2 = k1:Max
         O[k1+1,k2+1] = quadgk(r-> Overlap(k1,k2,L,r),0,50,rtol=1e-4)[1]
         Hamil1[k1+1,k2+1] = quadgk(r-> energy(k1,k2,L,r),0,50,rtol=1e-4)[1]
         Hamil2[k1+1,k2+1] = quadgk(r-> Kinetic(k1,k2,L,r),0,50,rtol=1e-4)[1]
     end
  end
  for k1 = 0:Max
     for k2 = k1:Max
        O[k2+1,k1+1] = O[k1+1,k2+1]
        Hamil1[k2+1,k1+1] = Hamil1[k1+1,k2+1]
        Hamil2[k2+1,k1+1] = Hamil2[k1+1,k2+1]
     end
  end
  Hamil = Hamil1 + Hamil2
  print(eigvals(Hamil,O))
end


function loopNum(Max,L)
   O = zeros((Max+1,Max+1))
   Hamil1 = zeros((Max+1,Max+1))
   Hamil2 = zeros((Max+1,Max+1)) 
   x=0.01:0.01:50
   for k1 = 0:Max
      for k2 = k1:Max
         O[k1+1,k2+1] = sum(Overlap.(k1,k2,L,x))*0.01
         Hamil1[k1+1,k2+1] = sum(energy.(k1,k2,L,x))*0.01
         Hamil2[k1+1,k2+1] = sum(Kinetic.(k1,k2,L,x))*0.01
     end
  end
  for k1 = 0:Max
     for k2 = k1:Max
        O[k2+1,k1+1] = O[k1+1,k2+1]
        Hamil1[k2+1,k1+1] = Hamil1[k1+1,k2+1]
        Hamil2[k2+1,k1+1] = Hamil2[k1+1,k2+1]
     end
  end
  Hamil = Hamil1 + Hamil2
  print(eigvals(Hamil,O))
end

#loop(Max,L)
#print("VVVVVVVVVVVVV")
loopNum(Max,L)

