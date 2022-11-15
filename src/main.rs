// -------------------------------------------------------------------------------------------------
// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// This code features the research article:
//
//       "Antagonism between viral infection and innate immunity at the single-cell level"
//
//                               by Frederic Grabowski et al.
//                                [TODO:JOURNAL-NAME], 202X
//
// The simulation mimicks the innate immune response to an infection with an RNA virus.
// The hard-coded and externally parametrized interactions between host cell and virus
// are specific to the respiratory syncytial virus (RSV). Infected cells attempt to produce
// and secrete interferon, which alerts the non-infected bystander cells about the nearby
// threat. The simulator executes alternating phases of (deterministic) interferon diffusion
// and (stochastic) chemical kinetics.
//
// For more info, see file ReadMe.md.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).
// -------------------------------------------------------------------------------------------------

mod cell;
mod commands;
mod config;
mod event;
mod lattice;
mod molecule;
mod protocol;
mod randomness;
mod rates;
mod simulation;
mod units;

use config::THREAD_STACK_SIZE;
use lattice::Lattice;
use protocol::Protocol;
use randomness::initialize_generator;
use rates::Rates;

use std::env;

fn print_usage_info() -> bool {
    if env::args().len() == 1 || env::args().any(|x| x == "-h" || x == "--help") {
        println!("Usage:");
        let exe_path = &env::args().collect::<Vec<_>>()[0];
        for invocation in [
            [ exe_path, "[parameters JSON file] [protocol file] <-i|--images>"],
            [ exe_path, "[-h|--help]"],
            [ exe_path, "[-v|--version]"],
        ] {
            println!(" {}",
                     invocation.into_iter().map(|s|s.to_string()).collect::<Vec<_>>().join(" "));
        }
        return true;
    }
    false
}

fn print_version_info() -> bool {
    if env::args().any(|x| x == "-v" || x == "--version") {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return true;
    }
    false
}

fn execute_protocol() -> bool {
    let argv = env::args().collect::<Vec<String>>();
    let rates = Rates::from_json_file(&argv[1]);
    let protocol = Protocol::from_text_file(&argv[2]);
    let images_out = env::args().any(|x| x == "-i" || x == "--images");

    std::thread::Builder::new()
        .name("protocol_execution".into())
        .stack_size(THREAD_STACK_SIZE)
        .spawn(move || {
            let mut generator = initialize_generator();
            let mut lattice = Lattice::new(&mut generator);
            protocol.execute(&mut lattice, &rates, &mut generator, images_out);
        })
        .expect("☠ @ protocol_execution thread")
        .join()
        .expect("☠ @ threads join");
    true
}

fn main() {
    let _ = print_usage_info() || print_version_info() || execute_protocol();
}
