// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

pub enum Mol {
    Vinf,  // viral infective particles (presence of the just-entered virus)
    Vrna,  // viral RNA
    Vprot, // viral proteins
    Pirf3, // phospho-IRF3
    Ifni,  // intracellular interferon (predominantly beta but also lambda to a some extent)
    Pstat, // phospho-STAT1/2
    Isg,   // (generalized) proteins of interferon-stimulated genes (accounts for RIG-I, PKR, OAS1)
}

pub const N_MOLECULE_SPECIES: usize = 7;
