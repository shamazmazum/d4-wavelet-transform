D4 wavelet transform implementation for Common Lisp
====================

## What is a wavelet transform?

Wavelet transform is just like Fourier transform defined by
![equation](https://latex.codecogs.com/gif.latex?F(z)%20%3D%20%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7Df(x)e%5E%7Bixz%7Ddx).

This transform gives you a "frequency representation" of a function, but as you can see, you lose all time information (i.e. you cannot tell where exaclty a
signal of some certain frequency is localized in time). Wavelet transform saves both time and frequency information. It is defined so:

![equation](https://latex.codecogs.com/gif.latex?F(a%2Cb)%20%3D%20%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7Df(x)%5Cpsi(%7B%7Bx-b%7D%20%5Cover%20%7Ba%7D%7D)dx) 

with some function ψ Usually, either ψ(x) has compact support or Ψ(z) (its Fourier-transformed
version) does (they cannot do both). Its usually desirable for ψ(x) to have compact support and for
Ψ(z) to have some fast decay. As you can see now, F(a,b) retains both frequency (provided by scaling
ψ(x) by a) and time (provided by shifting ψ(x) by b) information. Here I will not tell what ψ(x) are
good for exact reconstruction of f(x), nor I will tell how to do it.

There is also a discrete wavelet transform, where you calculate a scalar product of your function
f(x) with a set of functions
![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)%20%3D%20a_%7B0%7D%5E%7B-m%2F2%7D%5Cpsi(a_%7B0%7D%5E%7B-m%7Dx%20-%20nb_%7B0%7D))
where m,n are integer numbers. With a proper ψ(x) you can then reconstruct f(x) (exactly or with
some desired precision). Daubechies D4 wavelet is such a function ψ(x) with compact support for
which ![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)) constitute an
orthonormal basis on ![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D).

## What is Multiresolution Analysis?

Consider ![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)) defined by
![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)%20%3D%202%7D%5E%7B-m%2F2%7D%5Cpsi(2%5E%7B-m%7Dx%20-%20n))
This is the same ![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)) defined
earlier with
![equation](https://latex.codecogs.com/gif.latex?a_%7B0%7D%20%3D%202%2C%20b_%7B0%7D%20%3D%201).

Multiresolution Analysis is a set of spaces
![equation](https://latex.codecogs.com/gif.latex?%5Cldots%2CV_%7B-2%7D%2CV_%7B-1%7D%2CV_%7B0%7D%2CV_%7B1%7D%2CV_%7B2%7D%2C%5Cldots)
all of them being some closed subspace of
![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D). They are called multiresolution
analysis if:

 - ![equation](https://latex.codecogs.com/gif.latex?V_%7Bn%2B1%7D%20%5Csubset%20V_%7Bn%7D)
 - ![equation](https://latex.codecogs.com/gif.latex?%5Cbigcup_%7Bn%7D%20V_%7Bn%7D) is dense in ![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D)
 - ![equation](https://latex.codecogs.com/gif.latex?%5Cbigcap_%7Bn%7D%20V_%7Bn%7D%20%3D%20%5Cleft%5C%7B0%5Cright%5C%7D)
 - ![equation](https://latex.codecogs.com/gif.latex?f%20%5Cin%20V_%7Bn%7D) if and only if ![equation](https://latex.codecogs.com/gif.latex?f(2%5E%7Bn%7D%5Ccdot)%20%5Cin%20V_%7B0%7D)
 - There exists a function
   ![equation](https://latex.codecogs.com/gif.latex?%5Cphi%20%5Cin%20V_%7B0%7D) such that
   ![equation](https://latex.codecogs.com/gif.latex?%5Cleft%5C%7B%5Cphi_%7B0%2Ck%7D%3A%20k%20%5Cin%20%5Cmathbb%7BZ%7D%5Cright%5C%7D)
   is an orthonormal basis in ![equation](https://latex.codecogs.com/gif.latex?V_%7B0%7D).

![equation](https://latex.codecogs.com/gif.latex?%5Cphi_%7Bm%2Cn%7D(x)) are defined by
![equation](https://latex.codecogs.com/gif.latex?%5Cphi_%7Bm%2Cn%7D(x)%20%3D%202%7D%5E%7B-m%2F2%7D%5Cphi(2%5E%7B-m%7Dx%20-%20n)). For
every ![equation](https://latex.codecogs.com/gif.latex?V_%7Bn%7D) there is its orthonormal
complement ![equation](https://latex.codecogs.com/gif.latex?W_%7Bn%7D) in
![equation](https://latex.codecogs.com/gif.latex?V_%7Bn-1%7D) so
![equation](https://latex.codecogs.com/gif.latex?V_%7Bn-1%7D%20%3D%20V_%7Bn%7D%20%5Coplus%20W_%7Bn%7D),
where ![equation](https://latex.codecogs.com/gif.latex?W_%7Bn%7D) are orthogonal to each other. The
whole functional space is then can be represented as
![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D%20%3D%20%5Cbigoplus_%7Bn%3D-%5Cinfty%7D%5E%7B%5Cinfty%7DW_%7Bn%7D). Usually
we do not care about all ![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D) and some
"fine scaled" ![equation](https://latex.codecogs.com/gif.latex?V_%7Bn%7D) which will give "close"
representation of out signal is OK for us. Then we can get
![equation](https://latex.codecogs.com/gif.latex?V_%7Bn%7D%20%3D%20V_%7BN%7D%20%5Coplus%20(%5Cbigoplus_%7Bk%3D0%7D%5E%7BN-n-1%7DW_%7BN-k%7D)).
Scaled versions of ψ(x), ![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)) are
in 
![equation](https://latex.codecogs.com/gif.latex?W_%7Bm%7D) respectively and shifted versions of
![equation](https://latex.codecogs.com/gif.latex?%5Cpsi_%7Bm%2Cn%7D(x)) with some fixed m form a
basis in that ![equation](https://latex.codecogs.com/gif.latex?W_%7Bm%7D). Together with ϕ(x) on
some "coarse" scale, these finer scaled and shifted versions of ψ(x) form a basis in some space
![equation](https://latex.codecogs.com/gif.latex?V_%7Bn%7D) which is "closer" to
![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D) when more finer scales of ψ(x) are
taken into consideration.

With all proof skipped I will state that if
![equation](https://latex.codecogs.com/gif.latex?%5Cphi%20(x)%20%3D%20%5Csum_%7Bn%3D-%5Cinfty%7D%5E%7B%5Cinfty%7Dh_%7Bn%7D%5Cphi(2x-n))
is in ![equation](https://latex.codecogs.com/gif.latex?V_%7B0%7D), then there is a multiresolution
analysis for it and
![equation](https://latex.codecogs.com/gif.latex?%5Cpsi%20(x)%20%3D%20%5Csum_%7Bn%3D-%5Cinfty%7D%5E%7B%5Cinfty%7Dd_%7Bn%7D%5Cphi(2x-n))
where
![equation](https://latex.codecogs.com/gif.latex?d_%7Bn%7D%20%3D%20(-1)%5E%7Bn-1%7Dh_%7B-n-1%7D) is
in ![equation](https://latex.codecogs.com/gif.latex?W_%7B0%7D).

The best example (because of its simplicity) is ![equation](https://latex.codecogs.com/gif.latex?%5Cphi(x)%20%3D%20%5Cchi%20%5B0%3B%201%5D) (equals one on [0;1], zero
otherwise). This particular
![equation](https://latex.codecogs.com/gif.latex?%5Cphi(x)%20%3D%20%5Cphi(2x)%20%2B%20%5Cphi(2x%20-%201)). Then,
together with
![equation](https://latex.codecogs.com/gif.latex?%5Cpsi(x)%20%3D%20%5Cphi(2x)%20-%20%5Cphi(2x%20-%201))
they constitute an orthonormal basis
![equation](https://latex.codecogs.com/gif.latex?%5Cphi(x)%2C%20%5Cpsi_%7Bm%2Cn%7D(x)) on
![equation](https://latex.codecogs.com/gif.latex?L%5E%7B2%7D%5B0%3B1%5D). This is called Haar
wavelet.

ϕ(x) is called a scaling function and ψ(x) is a mother wavelet.

## What is Daubechies D4 wavelet?

D4 wavelet is a compactly supported orthonormal wavelet. Because it is compactly supported, there
is only a finite number (four exactly) of non-zero coefficients in
![equation](https://latex.codecogs.com/gif.latex?%5Cphi%20(x)%20%3D%20%5Csum_%7Bn%3D-%5Cinfty%7D%5E%7B%5Cinfty%7Dh_%7Bn%7D%5Cphi(2x-n))
so in other words
![equation](https://latex.codecogs.com/gif.latex?%5Cphi%20(x)%20%3D%20h_%7B0%7D%5Cphi(2x)%20%2B%20h_%7B1%7D%5Cphi(2x%20%2B%201)%20%2B%20h_%7B2%7D%5Cphi(2x%20%2B%202)%20%2B%20h_%7B3%7D%5Cphi(2x%20%2B%203))

I'll skip all derivations and just say that
![equation](https://latex.codecogs.com/gif.latex?h_%7Bn%7D) are composed so that
![equation](https://latex.codecogs.com/gif.latex?%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7D%20%5Cpsi(x)%20dx%20%3D%200)
and
![equation](https://latex.codecogs.com/gif.latex?%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7D%20x%5Cpsi(x)%20dx%20%3D%200). It's
said, that ψ(x) has two vanishing moments. Generally speaking if your input signal is at least n
times differetiable then transform coefficients decay to zero faster when moving to a finer scale
with a wavelet having n vanishing moments compared to wavelet having less than n vanishing
moments. Also it has a property that
![equation](https://latex.codecogs.com/gif.latex?%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7D%20%5Cphi(x)%20dx%20%3D%201)

## So what is this all about?

Here, we are close to cause of our short intro to wavelets. Namely, we want to transform a signal to
its wavelet representation and back. As for any orthonormal basis normally we would calculate scalar
product of our signal with the functions constituting our basis. But for Daubechies wavelet we do
not even know neither ψ(x) nor ϕ(x). All we know is a recursive equation above, called a dilation
equation. As far as I know, there is no way to represent scaling function or mother wavelet in terms
of elementary functions. So what we can do with dilation equation? First of all, we have our signal
represented as an array of samples taken a short period of time one from another. If this period is
short enough we can consider our signal constant within it. Then
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x-n)%5Crangle%20%3D%20f_%7Bn%7D)
because of the property of ϕ(x) shown above. Consider ϕ(x) is a scaling function on the finest scale
possible. Then we need to find
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x/2-n)%5Crangle)
and
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cpsi(x/2-n)%5Crangle),
which are scalar products with the scaling function and the mother wavelet on some "coarser"
scale. Here, we use our dilation equation to see that we need to convolve our input signal with
coefficients ![equation](https://latex.codecogs.com/gif.latex?h_%7Bn%7D) and
![equation](https://latex.codecogs.com/gif.latex?d_%7Bn%7D) for that and then retain only even
samples in the resulting sequence. So for n samples we get n/2 scalar products with the scaling
function and n/2 scalar products with the mother wavelet. We repeat the process using
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x/2-n)%5Crangle) as an
input data for the next iteration until we get only one
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x%2F2%5E%7Bm%7D)%20%5Crangle)
representing scalar product with the scaling function on the most coarse scale.

Inverse transform is done in the exactly reverse order. We use
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x%2F2%5E%7Bm%7D)%20%5Crangle)
and
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cpsi(x%2F2%5E%7Bm%7D)%20%5Crangle)
to get two values
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x%2F2%5E%7Bm-1%7D-n)%20%5Crangle)
and so on until we restore
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cphi(x-n)%5Crangle%20%3D%20f_%7Bn%7D).

## Lisp implementation.

In this lisp library there is `d4-wavelet-transform` package which consists of `transform` and
`inverse-transform` functions. Also, remember, that any orthonormal basis weakly converges to zero,
so expect values for
![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20f(x)%2C%20%5Cpsi(x-n)%20%5Crangle)
(scalar product with the mother wavelet on the finest scale) to be almost zero. You can "lose some
precision" with `lose-precision`. This function zeroes all coefficients which are close enough to
zero and therefore imitates some compression of input signal. There is `generate-sin` function,
which can get you sinusoidal signal, on which you can test `transform`. Also there is `array2file`
function which can store your signal in the form which is understood by gnuplot.

## Useful reading
 - [A Linear Algebra View of the Wavelet Transform](http://www.bearcave.com/misl/misl_tech/wavelets/matrix/index.html)
 - Ten Lectures on Wavelets, Ingrid Daubechies
 - Introduction to Hilbert Spaces with Applications, Third Edition by Lokenath Debnath, Piotr
 Mikusinski
