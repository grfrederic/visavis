// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use crate::lattice::Lattice;
use crate::molecule::Mol::*;
use crate::rates::Rates;

#[derive(Debug, Copy, Clone)]
pub enum Event {
    VinfIncr,
    VrnaIncr,
    VprotIncr,
    Pirf3Incr,
    IfniIncr,
    PstatIncr,
    IsgIncr,
    VinfDecr,
    VrnaDecr,
    VprotDecr,
    Pirf3Decr,
    IfniDecr,
    PstatDecr,
    IsgDecr,
    Die,
}

impl Event {
    pub fn rate_coef(self, rates: &Rates) -> f64 {
        match self {
            Event::VinfIncr => rates.vinf_incr,
            Event::VinfDecr => rates.vinf_decr,
            Event::VrnaIncr => rates.vrna_incr,
            Event::VrnaDecr => rates.vrna_decr,
            Event::VprotIncr => rates.vprot_incr,
            Event::VprotDecr => rates.vprot_decr,
            Event::Pirf3Incr => rates.pirf3_incr,
            Event::Pirf3Decr => rates.pirf3_decr,
            Event::IfniIncr => rates.ifni_incr,
            Event::IfniDecr => rates.ifni_decr,
            Event::PstatIncr => rates.pstat_incr,
            Event::PstatDecr => rates.pstat_decr,
            Event::IsgIncr => rates.isg_incr,
            Event::IsgDecr => rates.isg_decr,
            Event::Die => rates.die,
        }
    }

    pub fn occur(event_i: usize, lattice: &mut Lattice, cell_i: usize) -> Vec<usize> {
        let cell = &mut lattice.cells[cell_i];
        let mols = &mut cell.molecules;
        let neighs = &lattice.neighborhoods[cell_i];
        macro_rules! increment {
            ($m:ident) => {
                mols[$m as usize] += 1;
            };
        }
        macro_rules! decrement {
            ($m:ident) => {
                mols[$m as usize] -= 1;
            };
        }
        macro_rules! current_cell {
            () => {
                vec![cell_i]
            };
        }
        macro_rules! current_cell_and_neighboring_cells {
            () => {
                vec![cell_i, neighs[0], neighs[1], neighs[2], neighs[3], neighs[4], neighs[5]]
            };
        }
        let event = Event::from_index(event_i);
        match event {
            Event::VinfIncr => {
                increment!(Vinf);
                current_cell!()
            }
            Event::VinfDecr => {
                decrement!(Vinf);
                current_cell!()
            }
            Event::VrnaIncr => {
                increment!(Vrna);
                current_cell!()
            }
            Event::VrnaDecr => {
                decrement!(Vrna);
                current_cell!()
            }
            Event::VprotIncr => {
                increment!(Vprot);
                current_cell_and_neighboring_cells!()
            }
            Event::VprotDecr => {
                decrement!(Vprot);
                current_cell_and_neighboring_cells!()
            }
            Event::Pirf3Incr => {
                increment!(Pirf3);
                current_cell!()
            }
            Event::Pirf3Decr => {
                decrement!(Pirf3);
                current_cell!()
            }
            Event::IfniIncr => {
                increment!(Ifni);
                current_cell!()
            }
            Event::IfniDecr => {
                decrement!(Ifni);
                current_cell!()
            }
            Event::PstatIncr => {
                increment!(Pstat);
                current_cell!()
            }
            Event::PstatDecr => {
                decrement!(Pstat);
                current_cell!()
            }
            Event::IsgIncr => {
                increment!(Isg);
                current_cell!()
            }
            Event::IsgDecr => {
                decrement!(Isg);
                current_cell!()
            }
            Event::Die => {
                cell.alive = false;
                mols.iter_mut().for_each(|x| *x = 0);
                current_cell_and_neighboring_cells!()
            }
        }
    }

    #[inline]
    fn from_index(event_i: usize) -> Event {
        match event_i {
            0 => Event::VinfIncr,
            1 => Event::VrnaIncr,
            2 => Event::VprotIncr,
            3 => Event::Pirf3Incr,
            4 => Event::IfniIncr,
            5 => Event::PstatIncr,
            6 => Event::IsgIncr,
            7 => Event::VinfDecr,
            8 => Event::VrnaDecr,
            9 => Event::VprotDecr,
            10 => Event::Pirf3Decr,
            11 => Event::IfniDecr,
            12 => Event::PstatDecr,
            13 => Event::IsgDecr,
            14 => Event::Die,
            _ => {
                panic!("☠ @ event_i={}", event_i)
            }
        }
    }

    #[inline]
    pub fn to_index(self) -> usize {
        self as usize
    }
}
