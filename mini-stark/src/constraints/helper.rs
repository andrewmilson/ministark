use super::Constraint;
use fast_poly::GpuField;

/// Returns zero only when a == b.
pub fn are_eq<F: GpuField>(a: &Constraint<F>, b: &Constraint<F>) -> Constraint<F> {
    a - b
}

/// Returns zero only when a == zero.
pub fn is_zero<F: GpuField>(a: &Constraint<F>) -> Constraint<F> {
    a.clone()
}

/// Returns zero only when a = zero || a == one.
pub fn is_binary<E: GpuField>(a: E) -> E {
    a * a - a
}
