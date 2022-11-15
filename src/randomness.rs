// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use rand::{rngs::StdRng, SeedableRng};

fn gen_seed_from_time() -> [u8; 32] {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap();
    macro_rules! wrap {
        ($a: expr) => {
            (($a as u128) % (u8::max_value() as u128)) as u8
        };
    }
    let mut s = [0u8; 32];
    s[0] = wrap!(now.as_millis());
    s[15] = wrap!(now.as_secs());
    s[31] = wrap!(now.as_nanos());
    s
}

pub fn initialize_generator() -> StdRng {
    SeedableRng::from_seed(gen_seed_from_time())
}
