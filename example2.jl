using LinearAlgebra, Polynomials, SpecialPolynomials
using QuadGK, Calculus, CalculusWithJulia
using Plots

Max = 5
E = 1

Hamil = zeros((Max,Max))
Hamil[1,1] = 1
Hamil[Max,Max] = 1
x =1:Max

for u = 0:8 
   int = u/4
   Hamil[1,2] = int
   for k1 = 2:Max-1
      Hamil[k1,k1] = 1
      Hamil[k1,k1+1] = int
      Hamil[k1,k1-1] = int
   end
   Hamil[Max,Max-1]= int
   EIG = eigvecs(Hamil)
   for i = 1:Max
      EIG[:,i] = EIG[:,i].^4
      display(plot(x,EIG[:,i]))
   end
end 

#print(eigvecs(Hamil))
x = 1:10; y = rand(10); # These are the plotting data
plot(x,y)
display(plot(x, y))

readline()
