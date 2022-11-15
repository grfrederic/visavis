// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

// memory
pub const THREAD_STACK_SIZE: usize = 8 * 1024 * 1024; // larger lattice would require larger stack!

// output files
pub const OUT_FILE_NAME_TIME_IN_MIN: bool = true; // if false, then hours are used
