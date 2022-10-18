# GPU accelerated STARK engine for Apple devices

My time at Immutable left me fascinated with theoretical proof systems and their applications. I've developed an all- consuming amateur interest for scalable transparent arguments of knowledge (STARKs). They give the ability to prove to a 3rd party that some computation with $n$ steps ran correctly. Proofs can be verified in $O(log^{2}(n))$ vs naive $O(n)$. 

The library is primarily written in Rust but runs fast due to moving heavily parellisable polynomial arithmetic from the CPU to the GPU. Rust has poor support for coding GPUs so [Apple's metal shader language](https://developer.apple.com/metal/) was used instead.

# Things I don't like

- remembering what the longest table is
- all terminal and challence and column you have to remember the numerical index. Could implement trait and enums to mitigate this.

# TODO

- debugging memory table. Remove constraint for memory stay the same clock cycle increase
