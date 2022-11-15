// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use crate::cell::Cell;
use crate::config::OUT_FILE_NAME_TIME_IN_MIN;
use crate::event::Event;
use crate::lattice::{CytokineArray, Lattice};
use crate::molecule::Mol::{Vinf, Vrna, Vprot, Pirf3, Ifni, Pstat, Isg};
use crate::molecule::N_MOLECULE_SPECIES;
use crate::rates::{Rates, TIMESTEP};
use crate::rates::transport::{K_IFNE_LL_DT, K_IFNE_LU_DT, K_IFNE_UL_DT, K_IFNE_UU_DT};
use crate::units::{HOUR, MIN};

use rand::{rngs::StdRng, Rng};
use std::io::Write; // for .flush()
use threadpool::ThreadPool;

#[inline]
const fn ceil_pow2(i: u32) -> u32 {
    let mut v = i - 1;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v + 1 + ((v == 0) as u32)
}

const PROPENS_TREE_SIZE: usize =
    ceil_pow2(Lattice::CAPACITY as u32) as usize + Lattice::CAPACITY - 1;
const PROPENS_TREE_CELL_INDEX_BASE: usize = PROPENS_TREE_SIZE - Lattice::CAPACITY;
const PROPENS_EVENTS_SIZE: usize = 2 * N_MOLECULE_SPECIES + 1; // +molecule,-molecule, and 1 for Die
type Propensities = [[f64; PROPENS_EVENTS_SIZE]; PROPENS_TREE_SIZE];

pub struct Simulation {}

impl Simulation {
    #[inline]
    fn unset_cell_event_prop(propens: &mut Propensities, cell_i: usize, event_i: usize) {
        let mut propens_i = PROPENS_TREE_CELL_INDEX_BASE + cell_i;
        let rate = propens[propens_i][event_i];
        debug_assert!(rate >= 0.);
        if rate > 0. {
            loop {
                propens[propens_i][event_i] -= rate;
                if propens_i == 0 { break; } else { propens_i = (propens_i - 1) / 2 }
            }
        }
    }

    #[inline]
    fn unset_cell_events_props(propens: &mut Propensities, cell_i: usize) {
        for event_i in 0..PROPENS_EVENTS_SIZE {
            Simulation::unset_cell_event_prop(propens, cell_i, event_i)
        }
    }

    #[inline]
    fn set_event_propensity(propens: &mut Propensities, cell_i: usize, event_i: usize, rate: f64) {
        let mut propens_i = PROPENS_TREE_CELL_INDEX_BASE + cell_i;
        loop {
            propens[propens_i][event_i] += rate;
            if propens_i == 0 { break; } else { propens_i = (propens_i - 1) / 2 }
        }
    }

    fn set_cell_events_props(
        propens: &mut Propensities,
        lattice: &Lattice,
        rates: &Rates,
        cell_i: usize,
        ifni_secretion: bool,
    ) {
        let &cell = &lattice.cells[cell_i];
        if !cell.alive {
            if cfg!(debug_assertions) {
                for event_i in 0..PROPENS_EVENTS_SIZE {
                    debug_assert!(
                        propens[PROPENS_TREE_CELL_INDEX_BASE + cell_i][event_i].abs() < 1.0e-6
                    );
                }
            }
            return;
        }

        let &ms = &cell.molecules;

        macro_rules! mol_count_   { ($m:ident) => { ms[$m as usize] }; }
        macro_rules! mol_count    { ($m:ident) => { ms[$m as usize] as f64 }; }
        macro_rules! is_active    { ($m:ident) => { Cell::is_active($m, &ms) }; }
        macro_rules! can_increase { ($m:ident) => { Cell::can_increase($m, &ms) }; }
        macro_rules! can_decrease { ($m:ident) => { Cell::can_decrease($m, &ms) }; }

        macro_rules! set_ev_prop {
            ($rxn:ident, $rate_mul:expr, $rate_add:expr) => {
                let r = Event::$rxn;
                let rate = r.rate_coef(rates) * $rate_mul + $rate_add;
                Simulation::set_event_propensity(propens, cell_i, r.to_index(), rate);
            };
            ($rxn:ident, $rate_mul:expr) => {
                set_ev_prop!($rxn, $rate_mul, 0.)
            };
            ($rxn:ident) => {
                set_ev_prop!($rxn, 1., 0.)
            };
        }

        // The wiring of regulary interactions is expressed in terms of rates in the section below.
        //------------------------------------------------------------------------------------------
        let (vprot, isg) = (mol_count!(Vprot), mol_count!(Isg)); // extract often used values

        // Vinf
        if can_increase!(Vinf) {
            for neigh_cell_i in lattice.neighborhoods[cell_i].iter() {
                let neigh_ms = lattice.cells[*neigh_cell_i].molecules;
                if Cell::is_active(Vprot, &neigh_ms) {
                    set_ev_prop!(VinfIncr);
                }
            }
        }
        if can_decrease!(Vinf) && (mol_count_!(Vrna) == 0) {
            set_ev_prop!(VinfDecr);
        }

        // Vrna
        if is_active!(Vinf) && can_increase!(Vrna) {
            set_ev_prop!(VrnaIncr, mol_count!(Vinf) / (isg * rates.isg_inh_vrna + 1.));
        }
        if can_decrease!(Vrna) {
            set_ev_prop!(VrnaDecr);
        }

        // Vprot
        if is_active!(Vrna) && can_increase!(Vprot) {
            set_ev_prop!(VprotIncr, 1. / (isg * rates.isg_inh_vprot + 1.));
        }
        if can_decrease!(Vprot) {
            set_ev_prop!(VprotDecr);
        }

        // Pirf3
        if is_active!(Vrna) && can_increase!(Pirf3) {
            set_ev_prop!(
                Pirf3Incr,
                (isg * rates.isg_pro_pirf3 + 1.) / (vprot * rates.vprot_inh_pirf3 + 1.)
            );
        }
        if can_decrease!(Pirf3) {
            set_ev_prop!(Pirf3Decr);
        }

        // Ifni
        if is_active!(Pirf3) && can_increase!(Ifni) {
            set_ev_prop!(IfniIncr, 1. / (vprot * rates.vprot_inh_ifni + 1.));
        }
        if can_decrease!(Ifni) && ifni_secretion {
            set_ev_prop!(IfniDecr);
        }

        // Pstat
        if can_increase!(Pstat) {
            let ifne_lo = &lattice.cytokines[cell_i][0];
            // [CONSISTENCY: 0x19cfa3]
            set_ev_prop!(
                PstatIncr,
                ifne_lo / (rates.mm_pstat + ifne_lo) / (vprot * rates.vprot_inh_pstat + 1.)
            );
        }
        if can_decrease!(Pstat) {
            set_ev_prop!(PstatDecr);
        }

        // Isg
        if can_increase!(Isg) {
            set_ev_prop!(IsgIncr, mol_count!(Pstat), rates.k_isg0);
        }
        if can_decrease!(Isg) {
            set_ev_prop!(IsgDecr);
        }

        // cell death
        if is_active!(Vprot) {
            set_ev_prop!(Die);
        }
        //------------------------------------------------------------------------------------------
    }

    fn reset_cells_ifn_events_props(propens: &mut Propensities, lattice: &Lattice, rates: &Rates) {
        let r = Event::PstatIncr;
        let (event_i, rate_k) = (r.to_index(), r.rate_coef(rates));
        for cell_i in 0..Lattice::CAPACITY {
            let &cell = &lattice.cells[cell_i];
            if !cell.alive {
                continue;
            }
            Simulation::unset_cell_event_prop(propens, cell_i, event_i);
            if Cell::can_increase(Pstat, &cell.molecules) {
                let vprot = cell.molecules[Vprot as usize] as f64;
                let ifne_lo = lattice.cytokines[cell_i][0];
                Simulation::set_event_propensity(
                    propens,
                    cell_i,
                    event_i,
                    rate_k * ifne_lo // [CONSISTENCY:0x19cfa3]
                    / (rates.mm_pstat + ifne_lo)
                    / (vprot * rates.vprot_inh_pstat + 1.),
                )
            }
        }
    }

    fn compute_propensities(
        lattice: &Lattice,
        rates: &Rates,
        ifni_secretion: bool,
    ) -> Propensities {
        let mut propens: Propensities = [[0.; PROPENS_EVENTS_SIZE]; PROPENS_TREE_SIZE];
        for cell_i in 0..lattice.cells.len() {
            Simulation::set_cell_events_props(&mut propens, lattice, rates, cell_i, ifni_secretion)
        }
        propens
    }

    fn find_event(propens: &Propensities, rho: f64) -> (usize, usize) {
        // select event class
        let mut acc = 0.;
        let mut event_i = 0;
        for ei in 0..PROPENS_EVENTS_SIZE {
            acc += propens[0][ei];
            if acc > rho {
                break;
            } else {
                event_i += 1
            }
        }

        // reuse random number
        let mut rho2 = rho - (acc - propens[0][event_i]);
        debug_assert!(rho2 < propens[0][event_i]);

        // select cell
        let mut cell_i = 0; // in-tree
        while cell_i < PROPENS_TREE_CELL_INDEX_BASE {
            let next_left = 2 * cell_i + 1;
            let next_left_psum = propens[next_left][event_i];
            if rho2 < next_left_psum {
                cell_i = next_left
            } else {
                rho2 -= next_left_psum;
                cell_i = next_left + 1
            }
        }

        debug_assert!(propens[cell_i][event_i] > 0.);
        (cell_i - PROPENS_TREE_CELL_INDEX_BASE, event_i)
    }

    fn ifn_transport_step(lattice: &mut Lattice, rates: &Rates, ifni_secretion: bool) {
        let prev: CytokineArray = lattice.cytokines.clone();
        let q_ifne_dt = rates.q_ifne * TIMESTEP;
        let k_ifn_sec_dt = rates.k_ifn_sec * TIMESTEP;
        for (cell_i, neighs) in lattice.neighborhoods.iter().enumerate() {
            let (prev_lo, prev_hi) = (prev[cell_i][0], prev[cell_i][1]);
            let (mut lo, mut hi) = (prev_lo, prev_hi);

            // decay
            lo -= q_ifne_dt * prev_lo;
            hi -= q_ifne_dt * prev_hi;

            // secretion
            if ifni_secretion && Cell::is_active(Ifni, &lattice.cells[cell_i].molecules) {
                debug_assert!(lattice.cells[cell_i].alive);
                lo += k_ifn_sec_dt;
            }

            // transport: lower -> upper
            let l2g = K_IFNE_LU_DT * prev_lo;
            lo -= l2g;
            hi += l2g;

            // transport: upper -> lower
            let g2l = K_IFNE_UL_DT * prev_hi;
            hi -= g2l;
            lo += g2l;

            // transport: inter-neighbor exchange in both lower and upper
            lo -= (Lattice::N_NEIGHBORS as f64) * K_IFNE_LL_DT * prev_lo;
            hi -= (Lattice::N_NEIGHBORS as f64) * K_IFNE_UU_DT * prev_hi;
            for neigh_i in neighs.iter() {
                let prev_neigh_lo_hi = prev[*neigh_i];
                lo += K_IFNE_LL_DT * prev_neigh_lo_hi[0];
                hi += K_IFNE_UU_DT * prev_neigh_lo_hi[1];
            }
            lattice.cytokines[cell_i] = [lo, hi]
        }
    }

    pub fn simulate(
        lattice: &mut Lattice,
        rates: &Rates,
        rng: &mut StdRng,
        tspan: (f64, f64),
        files_out: bool,
        images_out: bool,
        files_out_interval: f64,
        ifni_secretion: bool,
        in_sep_thread: bool,
        init_frame_out: bool,
        workers: &Option<ThreadPool>,
    ) {
        // (currently, these 3 parameters are redundant)
        debug_assert!(in_sep_thread == workers.is_none());
        debug_assert!(in_sep_thread == !ifni_secretion);

        let mut propens = Simulation::compute_propensities(lattice, rates, ifni_secretion);
        let (mut t, mut t_next_ifn, mut t_next_files_out) = (
            tspan.0,
            tspan.0 + TIMESTEP,
            tspan.0 + (if init_frame_out { 0. } else { files_out_interval }),
        );
        if !in_sep_thread {
            if OUT_FILE_NAME_TIME_IN_MIN {
                print!("{:.0}m:", t / MIN);
            } else {
                print!("{:.0}h:", t / HOUR);
            }
            std::io::stdout().flush().unwrap()
        }
        loop {
            // if t >= t_next_print_out && !in_sep_thread { t_next_print_out += 1.*HOUR }
            if files_out && t >= t_next_files_out {
                if !in_sep_thread {
                    // spawn in a separate thread
                    print!(".");
                    std::io::stdout().flush().unwrap();
                    let (la, rr) = (lattice.clone(), rates.clone());
                    workers.as_ref().unwrap().execute(move || { la.out(t, &rr, images_out)});
                }
                t_next_files_out += files_out_interval;
            }
            if t >= tspan.1 {
                if OUT_FILE_NAME_TIME_IN_MIN {
                    print!(":{:.0}m ", tspan.1 / MIN);
                } else {
                    print!(":{:.0}h ", tspan.1 / HOUR);
                }
                std::io::stdout().flush().unwrap();
                break;
            }
            let sum_propens: f64 = propens[0].iter().sum();
            t += -(rng.gen_range(0.0..1.0) as f64).ln() / sum_propens; // exponential variate
            if t > t_next_ifn {
                t = t_next_ifn;
                t_next_ifn += TIMESTEP;
                Simulation::ifn_transport_step(lattice, rates, ifni_secretion);
                Simulation::reset_cells_ifn_events_props(&mut propens, lattice, rates);
            } else {
                let (cell_i, event_i) =
                    Simulation::find_event(&propens, rng.gen_range(0.0..sum_propens));
                for cell_j in Event::occur(event_i, lattice, cell_i).iter() {
                    Simulation::unset_cell_events_props(&mut propens, *cell_j);
                    Simulation::set_cell_events_props(
                        &mut propens,
                        lattice,
                        rates,
                        *cell_j,
                        ifni_secretion,
                    );
                }
            }
        } // loop
    } // simulate()
}
