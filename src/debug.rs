//! Tools for debugging issues that may arrive with AIR or STARK
//! TODO:

use crate::challenges::Challenges;
use crate::hints::Hints;
use crate::stark::Stark;
use crate::Matrix;

/// Checks AIR constraints are valid
pub const fn default_validate_constraints<S: Stark>(
    _this: &S,
    _challenges: &Challenges<S::Fq>,
    _hints: &Hints<S::Fq>,
    _base_trace: &Matrix<S::Fp>,
    _extension_trace: Option<&Matrix<S::Fq>>,
) {
    // ```text
    // TODO: move constraint checking from air.rs into here
    // #[cfg(all(feature = "std", debug_assertions))]
    // fn validate_constraints(
    //     &self,
    //     challenges: &Challenges<C::Fq>,
    //     hints: &Hints<C::Fq>,
    //     base_trace: &crate::Matrix<C::Fp>,
    //     extension_trace: Option<&crate::Matrix<C::Fq>>,
    // ) {
    //     use AlgebraicItem::*;
    //     use Expr::*;

    //     let num_execution_trace_columns = C::NUM_BASE_COLUMNS +
    // C::NUM_EXTENSION_COLUMNS;     let mut col_indicies = vec![false;
    // num_execution_trace_columns];     let mut challenge_indicies =
    // vec![false; challenges.len()];     let mut hint_indicies =
    // vec![false; hints.len()];

    //     for constraint in &self.constraints {
    //         constraint.traverse(&mut |node| match node {
    //             Leaf(Challenge(i)) => challenge_indicies[*i] = true,
    //             Leaf(Trace(i, _)) => col_indicies[*i] = true,
    //             Leaf(Hint(i)) => hint_indicies[*i] = true,
    //             _ => {}
    //         })
    //     }

    //     for (index, exists) in col_indicies.into_iter().enumerate() {
    //         if !exists {
    //             // TODO: make assertion
    //             println!("WARN: no constraints for execution trace column
    // {index}");         }
    //     }

    //     for (index, exists) in challenge_indicies.into_iter().enumerate()
    // {         if !exists {
    //             // TODO: make assertion
    //             println!("WARN: challenge at index {index} never used");
    //         }
    //     }

    //     for (index, exists) in hint_indicies.into_iter().enumerate() {
    //         if !exists {
    //             // TODO: make assertion
    //             println!("WARN: hint at index {index} never used");
    //         }
    //     }

    //     let trace_domain = self.trace_domain();
    //     let base_column_range = Self::base_column_range();
    //     let extension_column_range = Self::extension_column_range();

    //     // helper function to get a value from the execution trace
    //     let get_trace_value = |row: usize, col: usize, offset: isize| {
    //         let pos = (row as isize +
    // offset).rem_euclid(trace_domain.size() as isize) as usize;
    // if base_column_range.contains(&col) {
    // FieldVariant::Fp(base_trace.0[col][pos])         } else if
    // extension_column_range.contains(&col) {             let col =
    // col - C::NUM_BASE_COLUMNS;
    // FieldVariant::Fq(extension_trace.unwrap().0[col][pos])
    //         } else {
    //             unreachable!("requested column {col} does not exist")
    //         }
    //     };

    //     for (c_idx, constraint) in
    // self.constraints().into_iter().enumerate() {         for
    // (row, x) in trace_domain.elements().enumerate() {
    // let is_valid = constraint                 .check(&mut |leaf|
    // match leaf {                     X => FieldVariant::Fp(x),
    //                     &Hint(i) => FieldVariant::Fq(hints[i]),
    //                     &Challenge(i) => FieldVariant::Fq(challenges[i]),
    //                     &Trace(col, offset) => get_trace_value(row, col,
    // offset),                     &Constant(c) => c,
    //                 })
    //                 .is_some();

    //             if !is_valid {
    //                 let mut vals = vec![format!("x = {x}")];
    //                 constraint.traverse(&mut |node| match *node {
    //                     // get a description of each leaf node
    //                     Leaf(Trace(col, offset)) => vals.push(format!(
    //                         "Trace(col={col:0>3}, offset={offset:0>3}) =
    // {}",                         get_trace_value(row, col,
    // offset)                     )),
    //                     Leaf(Challenge(i)) => {
    //                         vals.push(format!("Challenge({i}) = {}",
    // challenges[i]))                     }
    //                     Leaf(Hint(i)) => vals.push(format!("Hint({i}) =
    // {}", hints[i])),                     // skip tree nodes
    //                     _ => (),
    //                 });

    //                 vals.sort();
    //                 vals.dedup();

    //                 // TODO: display constraint? eprintln!("Constraint
    // is:\n{constraint}\n");                 #[cfg(feature = "std")]
    //                 eprint!("Constraint {c_idx} does not evaluate to a
    // low degree polynomial. ");                 #[cfg(feature =
    // "std")]                 eprintln!("Divide by zero occurs at
    // row {row}.\n");                 #[cfg(feature = "std")]
    //                 eprintln!("Expression values:\n{}", vals.join("\n"));
    //                 panic!();
    //             }
    //         }
    //     }
    // }
    // ```
}
