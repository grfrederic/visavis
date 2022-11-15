// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use crate::units::MIN;

use std::fs;

use serde::{Deserialize, Serialize};
use serde_json::from_str;

pub const TIMESTEP: f64 = 0.1 * MIN;

// chemical reaction rates (stochastic)
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Rates {
    pub vinf_incr: f64,
    pub vinf_decr: f64,
    pub vrna_incr: f64,
    pub vrna_decr: f64,
    pub vprot_incr: f64,
    pub vprot_decr: f64,
    pub pirf3_incr: f64,
    pub pirf3_decr: f64,
    pub ifni_incr: f64,
    pub ifni_decr: f64,
    pub pstat_incr: f64,
    pub pstat_decr: f64,
    pub isg_incr: f64,
    pub isg_decr: f64,
    pub k_isg0: f64,
    pub mm_pstat: f64,
    pub die: f64,
    pub k_ifn_sec: f64,
    pub q_ifne: f64,
    pub vprot_inh_pirf3: f64,
    pub vprot_inh_ifni: f64,
    pub vprot_inh_pstat: f64,
    pub isg_inh_vrna: f64,
    pub isg_inh_vprot: f64,
    pub isg_pro_pirf3: f64,
}

impl Rates {
    pub fn from_json_file(params_filename: &String) -> Self {
        let contents = fs::read_to_string(params_filename).expect("☠ 🕮 JSON");
        from_str(&contents).unwrap()
    }
}

// transport of extracelluar interferon-beta (deterministic)
pub mod transport {
    use crate::lattice::Lattice;
    use crate::rates::TIMESTEP;
    use crate::units::MIN;

    pub const K_IFNE_LL_DT: f64 = 0.5 / (Lattice::N_NEIGHBORS as f64) / MIN * TIMESTEP;
    pub const K_IFNE_UU_DT: f64 = K_IFNE_LL_DT;
    pub const K_IFNE_LU_DT: f64 = 5. * 0.1 / MIN * TIMESTEP;
    pub const K_IFNE_UL_DT: f64 = 5. * 0.001 / MIN * TIMESTEP;
}
