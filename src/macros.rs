// Adapted from the `forward_ref_binop!` macro in the Rust standard library.
// Implements "&T op U", "T op &U" based on "&T op &U"
macro_rules! forward_ref_binop {
    (impl < F: GpuField > $imp:ident, $method:ident for $t:ty, $u:ty) => {
        impl<'a, F: GpuField> $imp<$u> for &'a $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, other: $u) -> <&'a $t as $imp<&'a $u>>::Output {
                $imp::$method(self, &other)
            }
        }

        impl<'a, F: GpuField> $imp<&'a $u> for $t {
            type Output = <&'a $t as $imp<&'a $u>>::Output;

            #[inline]
            fn $method(self, other: &$u) -> <&'a $t as $imp<&'a $u>>::Output {
                $imp::$method(&self, other)
            }
        }
    };
}
