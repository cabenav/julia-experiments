using LinearAlgebra, Polynomials, SpecialPolynomials
using QuadGK, Calculus, CalculusWithJulia
using GSL
using Plots 
using DelimitedFiles, DataFrames #CSV

Max =10
Q = 1
z = 1
a1 = 1
a2 = 1 
a3 = 1
a4 = 1
alpha = 1
rc = 1
L = 0

#alkali = readdlm(“alkaliparameters.csv”,‘,’);
#alkalidata = DataFrame(Li=alkali[:,1],Na=alkali[:,2],K=alkali[:,3],Rb=alkali[:,4],
#l=[NaN,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2],
#parameters=[“α_c”,“a_1”,“a_2",“a_3”,“a_4",“r_c”,“a_1",“a_2”,“a_3",“a_4”,“r_c”,“a_1”,“a_2",“a_3”,“a_4",“r_c”])

abstract type AbstractAtom end

struct Atom <: AbstractAtom
	a::Array{Float64}
	l::Int
	alpha_c::Float64
	r_c::Float64
end


function Zl(r,atom::Atom)
	a1, a2, a3, a4 = atom.a
        1 + (z-1)*exp(-a1*r)- r*(a3+a4*r)*exp(-a2*r)
end

function potential(r,atom::Atom)
	rc, alpha = atom.alpha_c, atom.r_c 
   	-Zl(r,atom)/r-alpha/(2*r^4)*(1-exp(-(r/rc)^6))
end

function potentialH(r)
   	-1/r
end

function angular(r,l)
	l*(l+1)/(2*r^2)
end

function Sturmian(r,k,l,Q=0.7)
	sqrt(factorial(k)/(2*factorial(2*l+1+k)))*exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k,2*l+1,2*r*Q)
end

function DSturmian(r,k,l,Q=0.7)
	if k < 1
	   return  sqrt(factorial(k)/(2*factorial(2*l+1+k)))*Q*(2*(1+l)*exp(-Q*r)*(2*Q*r)^(l)*sf_laguerre_n(k,2*l+1,2*r*Q)-exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k,2*l+1,2*r*Q))
   	else
           return  sqrt(factorial(k)/(2*factorial(2*l+1+k)))*Q*(-2*exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k-1,2*l+2,2*r*Q)+2*(1+l)*exp(-Q*r)*(2*Q*r)^(l)*sf_laguerre_n(k,2*l+1,2*r*Q)-exp(-Q*r)*(2*Q*r)^(l+1)*sf_laguerre_n(k,2*l+1,2*r*Q))
   	end
end


function Overlap(Num,l,Rang=50,Grid=1000,Q=0.7)
	x = LinRange(sqrt(0.01),sqrt(Rang),Grid).^2
        O = zeros(size(x)[1],Num)
	for j = 0:Num-1
	   O[:,j+1] = Sturmian.(x,j,l,Q) 
	end
	Overlap = transpose(O)*(O .* 2x)
        print(Overlap)
end

Overlap(5,0)

function Kinetic(k1,k2,L,r)
   DSturmian(r^2,k,l,Q)*DSturmian(r^2,k,l,Q)/2*2*r
end

function energy(k1,k2,L,r)
   Sturmian(k1,L,r^2)*Sturmian(k2,L,r^2)*(potentialH(r^2)+angular(L,r^2))*2*r
end

function Overlap(k1,k2,L,r)
   Sturmian(r^2,k1,L)*Sturmian(r^2,k2,L)*2*r
end

#exact calculation
#=
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
=#

#Numerical calculations
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

