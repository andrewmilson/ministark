// Adapted from the `forward_ref_binop!` macro in the Rust standard library.
// Implements "&T op U", "T op &U" based on "T op U"
#[macro_export]
macro_rules! forward_ref_binop {
    (
        impl
        $(< $($lt:tt $(: $clt:tt $(< $($elt:tt),+ >)? $(+ $dlt:tt $(< $($flt:tt),+ >)?)*)?),+ >)?
        $imp:ident,
        $method:ident for
        $t:ty,
        $u:ty
    ) => {
        impl< 'frb $(, $( $lt $( : $clt $(< $($elt),+ >)? $(+ $dlt $(< $($flt),+ >)? )* )? ),+ )? > $imp<$u> for &'frb $t {
            type Output = <$t as $imp<$u>>::Output;

            #[inline]
            fn $method(self, other: $u) -> <$t as $imp<$u>>::Output {
                $imp::$method(<$t as Clone>::clone(self), other)
            }
        }

        impl< 'frb $(, $( $lt $( : $clt $(< $($elt),+ >)? $(+ $dlt $(< $($flt),+ >)? )* )? ),+ )? > $imp<&'frb $u> for $t {
            type Output = <$t as $imp<$u>>::Output;
            // type Output = <&'frb $t as $imp<&'frb $u>>::Output;

            #[inline]
            fn $method(self, other: &$u) -> <$t as $imp<$u>>::Output {
                $imp::$method(self, <$u as Clone>::clone(other))
            }
        }
    };
}

// Adapted from the `forward_ref_op_assign!` macro in the Rust standard library.
// implements "T op= &U", based on "T op= U"
#[macro_export]
macro_rules! forward_ref_op_assign {
    (
        impl
        $(< $($lt:tt $(: $clt:tt $(< $($elt:tt),+ >)? $(+ $dlt:tt $(< $($flt:tt),+ >)?)*)?),+ >)?
        $imp:ident,
        $method:ident for
        $t:ty,
        $u:ty
    ) => {
        impl< 'frb $(, $( $lt $( : $clt $(< $($elt),+ >)? $(+ $dlt $(< $($flt),+ >)? )* )? ),+ )? > $imp<&'frb $u> for $t {
            #[inline]
            fn $method(&mut self, other: &'frb $u) {
                $imp::$method(self, other.clone());
            }
        }
    };
}
