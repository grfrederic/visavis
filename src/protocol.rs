// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use std::fs::File;
use std::io::{self, BufRead};
use std::num::ParseFloatError;
use std::path::Path;
use std::str::FromStr;

use rand::rngs::StdRng;

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{char, digit1, multispace1},
    combinator::{map_res, opt, recognize},
    error::ErrorKind,
    number::complete::double,
    sequence::{delimited, pair, separated_pair, terminated, tuple},
};

use crate::commands::{add_virus, remove_ifne, run_simulation, run_simulation_quietly, set_upper_ifne};
use crate::lattice::Lattice;
use crate::rates::Rates;
use crate::units::{MIN, HOUR, DAY, conversion};

pub struct Protocol {
    pub commands: Vec<String>,
}

impl Protocol {
    pub fn from_text_file(protocol_file_path: &String) -> Self {
        let mut lines = Vec::<String>::new();
        let protocol_file_path = Path::new(protocol_file_path);
        let file = File::open(protocol_file_path).expect("☠ 🕮 Protocol");
        let reader = io::BufReader::new(file);
        for line in reader.lines().flatten() {
            lines.push(line)
        }
        Protocol { commands: lines }
    }

    pub fn execute(
        &self,
        lattice: &mut Lattice,
        rates: &Rates,
        rng: &mut StdRng,
        out_images: bool,
    ) {
        let factor = || double::<&str, (_, ErrorKind)>;
        let number = || pair::<_, _, _, (_, ErrorKind), _, _>(opt(char('-')), digit1);

        let in_unit_of = |u| move |s| -> Result<f64, ParseFloatError> { Ok(u * f64::from_str(s)?) };
        let minutes = || map_res( terminated(recognize(number()), char('m')), in_unit_of(MIN));
        let hours = || map_res( terminated(recognize(number()), char('h')), in_unit_of(HOUR));
        let days = || map_res( terminated(recognize(number()), char('d')), in_unit_of(DAY));

        let time = || alt((minutes(), hours(), days()));
        let timespan = || separated_pair(time(), tag("..."), time());
        let every = || delimited(char('['), time(), char(']'));
        let never = || tag("[]");

        let cmd_run = || tuple((tag("run"), multispace1, timespan(), multispace1, every()));
        let cmd_run_quiet = || tuple((tag("run"), multispace1, timespan(), multispace1, never()));
        let cmd_set_ifne = || tuple((tag("=IFN"), multispace1, factor(), multispace1, tag("U/ml")));
        let cmd_del_ifn = || tag::<_, &str, (_, ErrorKind)>("!IFN");
        let cmd_add_rsv = || tuple((tag("+RSV"), multispace1, factor(), multispace1, tag("MOI")));

        let mut out_init_frame = false; // whether initial frame in output
        for command in self.commands.iter() {
            if let Ok((_, (_, _, tspan, _, dt))) = cmd_run()(&command) {
                run_simulation(lattice, rates, rng, tspan, out_images, dt, out_init_frame);
                out_init_frame = false;
            } else if let Ok((_, (_, _, tspan, _, _))) = cmd_run_quiet()(&command) {
                run_simulation_quietly(lattice, rates, rng, tspan, out_images, out_init_frame);
                out_init_frame = false;
            } else if let Ok((_, (_, _, ifne_uml, _, _))) = cmd_set_ifne()(&command) {
                set_upper_ifne(lattice, ifne_uml * conversion::IFNE_U_PER_ML_TO_MOLECULE_COUNT);
                out_init_frame = true;
            } else if let Ok(_) = cmd_del_ifn()(&command) {
                remove_ifne(lattice);
                out_init_frame = true;
            } else if let Ok((_, (_, _, moi, _, _))) = cmd_add_rsv()(&command) {
                add_virus(lattice, rng, moi);
                out_init_frame = true;
            } else {
                panic!("☠ @ command: {:?}", command);
            }
        }
        println!();
    }
}
