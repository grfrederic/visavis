// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

pub const MIN: f64 = 1.; // time unit is 1 minute
pub const HOUR: f64 = 60. * MIN;
pub const DAY: f64 = 24. * HOUR;

pub mod conversion {
    pub const IFNE_U_PER_ML_TO_MOLECULE_COUNT: f64 = 300.;
}
