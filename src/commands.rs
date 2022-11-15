// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use crate::config::THREAD_STACK_SIZE;
use crate::lattice::Lattice;
use crate::molecule::Mol;
use crate::rates::Rates;
use crate::simulation::Simulation;

use rand::rngs::StdRng;
use rand_distr::{Distribution, Poisson};

pub fn add_virus(lattice: &mut Lattice, rng: &mut StdRng, moi: f64) {
    let poisson = Poisson::new(moi).unwrap();
    lattice.cells.iter_mut().for_each(|c| {
        if c.alive {
            let vi = poisson.sample(rng) as u64;
            debug_assert!(vi < u8::max_value() as u64);
            let has_virus = vi > 0 || c.molecules[Mol::Vinf as usize] > 0;
            c.molecules[Mol::Vinf as usize] = if has_virus { 1 } else { 0 } as u8
        }
    });
}

pub fn set_upper_ifne(lattice: &mut Lattice, upper_ifne_mlcs: f64) {
    lattice.cytokines = [[0., upper_ifne_mlcs]; Lattice::CAPACITY]
}

pub fn remove_ifne(lattice: &mut Lattice) {
    lattice.cytokines = [[0., 0.]; Lattice::CAPACITY]
}

pub fn run_simulation_quietly(
    lattice: &mut Lattice,
    rates: &Rates,
    rng: &mut StdRng,
    tspan: (f64, f64),
    images_out: bool,
    init_frame_out: bool,
) {
    run_simulation_(
        lattice,
        rates,
        rng,
        tspan,
        /*files_out:*/ false,
        images_out,
        /*files_out_interval*/ -1.,
        init_frame_out,
    )
}

pub fn run_simulation(
    lattice: &mut Lattice,
    rates: &Rates,
    rng: &mut StdRng,
    tspan: (f64, f64),
    images_out: bool,
    files_out_interval: f64,
    init_frame_out: bool,
) {
    run_simulation_(
        lattice,
        rates,
        rng,
        tspan,
        /*files_out:*/ true,
        images_out,
        files_out_interval,
        init_frame_out,
    )
}

fn run_simulation_(
    lattice: &mut Lattice,
    rates: &Rates,
    rng: &mut StdRng,
    tspan: (f64, f64),
    files_out: bool,
    images_out: bool,
    files_out_interval: f64,
    init_frame_out: bool,
) {
    let workers = Some(
        threadpool::Builder::new()
            .num_threads(num_cpus::get())
            .thread_stack_size(THREAD_STACK_SIZE)
            .build(),
    );
    Simulation::simulate(
        lattice,
        rates,
        rng,
        tspan,
        files_out,
        images_out,
        files_out_interval,
        /*ifni_secretion:*/ true,
        /*in_sep_thread:*/ false,
        init_frame_out,
        &workers,
    );
    workers.unwrap().join()
}
