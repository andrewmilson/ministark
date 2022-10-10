pub trait BrainfuckColumn {
    const FIRST_TRACE_COL_INDEX: usize;
    const LAST_TRACE_COL_INDEX: usize;
    const NUM_TRACE_COLUMNS: usize = Self::LAST_TRACE_COL_INDEX - Self::FIRST_TRACE_COL_INDEX + 1;
}

pub enum Challenge {
    A,
    B,
    C,
    D,
    E,
    F,
    Alpha,
    Beta,
    Gamma,
    Delta,
    Eta,
}

impl mini_stark::constraint::Challenge for Challenge {
    fn index(&self) -> usize {
        match self {
            Self::A => 0,
            Self::B => 1,
            Self::C => 2,
            Self::D => 3,
            Self::E => 4,
            Self::F => 5,
            Self::Alpha => 6,
            Self::Beta => 7,
            Self::Gamma => 8,
            Self::Delta => 9,
            Self::Eta => 10,
        }
    }
}

pub enum ProcessorBaseColumn {
    Cycle,
    Ip,
    CurrInstr,
    NextInstr,
    Mp,
    MemVal,
    MemValInv,
}

pub enum ProcessorExtensionColumn {
    InstructionPermutation,
    MemoryPermutation,
    InputEvaluation,
    OutputEvaluation,
}

pub enum MemoryBaseColumn {
    Cycle,
    Mp,
    MemVal,
    Dummy,
}

pub enum MemoryExtensionColumn {
    Permutation,
}

pub enum InstructionBaseColumn {
    Ip,
    CurrInstr,
    NextInstr,
}

pub enum InstructionExtensionColumn {
    ProcessorPermutation,
    ProgramEvaluation,
}

pub enum InputBaseColumn {
    Value,
}

pub enum InputExtensionColumn {
    Evaluation,
}

pub enum OutputBaseColumn {
    Value,
}

pub enum OutputExtensionColumn {
    Evaluation,
}

impl BrainfuckColumn for ProcessorBaseColumn {
    const FIRST_TRACE_COL_INDEX: usize = ProcessorBaseColumn::Cycle as usize;
    const LAST_TRACE_COL_INDEX: usize = ProcessorBaseColumn::MemValInv as usize;
}

impl BrainfuckColumn for MemoryBaseColumn {
    const FIRST_TRACE_COL_INDEX: usize = ProcessorBaseColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::Dummy as usize;
}

impl BrainfuckColumn for InstructionBaseColumn {
    const FIRST_TRACE_COL_INDEX: usize = MemoryBaseColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::NextInstr as usize;
}

impl BrainfuckColumn for InputBaseColumn {
    const FIRST_TRACE_COL_INDEX: usize = InstructionBaseColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::Value as usize;
}

impl BrainfuckColumn for OutputBaseColumn {
    const FIRST_TRACE_COL_INDEX: usize = InputBaseColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::Value as usize;
}

impl BrainfuckColumn for ProcessorExtensionColumn {
    const FIRST_TRACE_COL_INDEX: usize = OutputBaseColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize =
        Self::FIRST_TRACE_COL_INDEX + Self::OutputEvaluation as usize;
}

impl BrainfuckColumn for MemoryExtensionColumn {
    const FIRST_TRACE_COL_INDEX: usize = ProcessorExtensionColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::Permutation as usize;
}

impl BrainfuckColumn for InstructionExtensionColumn {
    const FIRST_TRACE_COL_INDEX: usize = MemoryExtensionColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize =
        Self::FIRST_TRACE_COL_INDEX + Self::ProgramEvaluation as usize;
}

impl BrainfuckColumn for InputExtensionColumn {
    const FIRST_TRACE_COL_INDEX: usize = InstructionExtensionColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::Evaluation as usize;
}

impl BrainfuckColumn for OutputExtensionColumn {
    const FIRST_TRACE_COL_INDEX: usize = InputExtensionColumn::LAST_TRACE_COL_INDEX + 1;
    const LAST_TRACE_COL_INDEX: usize = Self::FIRST_TRACE_COL_INDEX + Self::Evaluation as usize;
}

macro_rules! impl_column {
    ($t:ty) => {
        impl mini_stark::constraint::Column for $t {
            fn index(&self) -> usize {
                Self::FIRST_TRACE_COL_INDEX + *self as usize
            }
        }
    };
}

impl_column!(ProcessorBaseColumn);
impl_column!(ProcessorExtensionColumn);

impl_column!(MemoryBaseColumn);
impl_column!(MemoryExtensionColumn);

impl_column!(InstructionBaseColumn);
impl_column!(InstructionExtensionColumn);

impl_column!(InputBaseColumn);
impl_column!(InputExtensionColumn);

impl_column!(OutputBaseColumn);
impl_column!(OutputExtensionColumn);
