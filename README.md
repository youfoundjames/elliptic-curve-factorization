## Different implementations of Lenstra's elliptic curve factorization algorithm ("breaking" ECC)


Here is some SageMath code I wrote at my cryptography REU at University of Kansas. It entails several different implementations of Lenstra's elliptic curve factoring algorithm, namely for different choices of coefficients _k_ as in step 3 in https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization. Included in the code is the ability to flip between different generative methods for these coefficients and compare their accuracy given different constraints placed on the randomly generated elliptic curves.

Lenstra's elliptic curve algorithm is the third-fastest known algorithm for general-purpose integer factorization, such as one might use to try to decrypt messages encoded using modern methods, such as RSA and elliptic curves. Most current encryption methods depend on the fact that it is currently difficult to factor large numbers, which are in practice typically "semi-primes" (i.e. the product of two prime numbers) with hundreds of digits. 
