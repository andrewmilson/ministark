use crate::Constraint;
use crate::ProofOptions;
use crate::TraceInfo;
use ark_serialize::CanonicalSerialize;
use fast_poly::GpuField;

pub trait Air {
    type BaseField: GpuField;
    type PublicInputs: CanonicalSerialize;

    // TODO: could make this borrow info and options if so inclined
    fn new(info: TraceInfo, inputs: Self::PublicInputs, options: ProofOptions) -> Self;

    fn trace_info(&self) -> &TraceInfo;

    fn boundary_constraints(&self) -> Vec<Constraint<Self::BaseField>> {
        Vec::new()
    }

    fn transition_constraints(&self) -> Vec<Constraint<Self::BaseField>> {
        Vec::new()
    }

    fn terminal_constraints(&self) -> Vec<Constraint<Self::BaseField>> {
        Vec::new()
    }
}
