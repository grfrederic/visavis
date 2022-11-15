// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use crate::molecule::Mol;
use crate::molecule::N_MOLECULE_SPECIES;

type MolArray = [u8; N_MOLECULE_SPECIES];

#[derive(Clone, Copy)]
pub struct Cell {
    pub alive: bool,
    pub molecules: MolArray,
}

macro_rules! level_bound {
    ($sm:ident, $mi:expr) => {
        Cell::$sm.molecules[$mi]
    };
}

impl Cell {
    // Array entries below correspond to molecules given in enum crate::molecule::Mol.

    pub const MIN: Cell = Cell {
        alive: true,
        molecules: [0, 0, 0, 0, 0, 0, 0],
    };

    pub const MAX: Cell = Cell {
        alive: true,
        molecules: [1, 3, 3, 3, 3, 3, 3],
    };

    pub const ACT: Cell = Cell {
        alive: true,
        molecules: [1, 3, 3, 3, 3, 1, 1],
    };

    #[inline]
    pub fn is_active(m: Mol, ms: &MolArray) -> bool {
        let i = m as usize;
        ms[i] >= level_bound!(ACT, i)
    }

    #[inline]
    pub fn can_increase(m: Mol, ms: &MolArray) -> bool {
        let i = m as usize;
        ms[i] < level_bound!(MAX, i)
    }

    #[inline]
    pub fn can_decrease(m: Mol, ms: &MolArray) -> bool {
        let i = m as usize;
        ms[i] > level_bound!(MIN, i)
    }
}
