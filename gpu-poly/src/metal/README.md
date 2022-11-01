## Adding support for your field

Right not only implementation for prime field arithmetic has been done although adding other field support is encouraged. Prime fields can easily be added. Take a look at the implementation of `FP270497897142230380135924736767050121217` and the comments on `Felt128`.

## GPU functions

### Number Theory Transform (FFT)

The FFT happens in-place and employs the [Cooley–Tukey algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm). As input it takes an array of values $a = [a_0,a_1,...,a_{n-1}]$ to be transformed and a constant array of twiddle factors

$$twiddles[i] = ω^{permute(i)}\quad\forall\ 0 \leq i < n/2$$

where $n \in 2^\mathbb{N}$ is the length of the input values $a$ and $ω$ is $n$-th root of unity such that $ω^n \equiv 1 \pmod{p}\  \land\  ω^i \not\equiv 1 \pmod{p}\ \ \forall\ 0 \lt i \lt n$. The function $permute(i)$ can be chosen depending on the order of the input values getting transformed. If the input values $a$ are in bit reversed order the $permute(i)$ function refers to the identity function. If the input values $a$ is in natural order the $permute(i)$ function performs a bit reversal of the first $\log_2(n)-1$ bits of input $i$.

Due to the nature of the Cooley-Tuckey algorithm the in-place transformed values are reshuffled (unmuffling results in non-negligible overhead). Either (1) bit-reversed order input values are transformed to natural order output values or (2) natural order input values are transformed to bit-reversed order output values. It's required that (1) has twiddle factors in natural order and (2) has twiddle factors in bit-reversed order.

There are two functions for performing the FFT:

- `FftMultiple` uses threadgroup memory (faster than global) to perform multiple FFT iterations
- `FftSingle` performs a single iteration without threadgroup memory

### Inverse Number Theory Transform (IFFT)

The IFFT happens in-place and employs the [Cooley–Tukey algorithm](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm) with the Gentleman-Sande butterfly structure. As input it takes an array of values $a^\prime = [a^\prime_0,a^\prime_1,...,a^\prime_{n-1}]$ to be transformed and a constant array of inverse twiddle factors

$$
inv\_twiddles[i] = \dfrac{1}{n}\ ω^{-permute(i)}\quad\forall\ 0 \leq i < n/2
$$

where $n\in 2^\mathbb{N}$ is the length of the input values $a^\prime$ and $ω$ is $n$-th root of unity such that $ω^n \equiv 1 \pmod{p}\  \land\  ω^i \not\equiv 1 \pmod{p}\ \ \forall\ 0 \lt i \lt n$. The function $permute(i)$ can be chosen depending on the order of the input values getting transformed. If the input values $a^\prime$ are in bit reversed order the $permute(i)$ function refers to the identity function. If the input values $a$ is in natural order the $permute(i)$ function performs a bit reversal of the first $\log_2(n)-1$ bits of input $i$.

Due to the nature of the Cooley-Tuckey algorithm the in-place transformed values are reshuffled (unmuffling results in non-negligible overhead). Either (1) bit-reversed order input values are transformed to natural order output values or (2) natural order input values are transformed to bit-reversed order output values. It's required that (1) has twiddle factors in natural order and (2) has twiddle factors in bit-reversed order.

There are two functions for performing the FFT:

- `IFftMultiple` uses threadgroup memory (faster than global) to perform multiple IFFT iterations
- `IFftSingle` performs a single IFFT iteration without threadgroup memory
