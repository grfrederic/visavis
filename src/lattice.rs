// VIS-A-VIS, a simulator of Viral Infection Spread And Viral Infection Self-containment.
//
// Copyright (2022) Marek Kochanczyk & Frederic Grabowski (IPPT PAN, Warsaw).
// Licensed under the 3-Clause BSD license (https://opensource.org/licenses/BSD-3-Clause).

use crate::cell::Cell;
use crate::config::OUT_FILE_NAME_TIME_IN_MIN;
use crate::molecule::{Mol, Mol::{Vinf, Vrna, Vprot, Pirf3, Pstat}, N_MOLECULE_SPECIES};
use crate::rates::Rates;
use crate::units::{MIN, HOUR};

use cairo::{Context, Format, ImageSurface};
use rand::{rngs::StdRng, seq::SliceRandom};
use std::f64::consts::PI;
use std::fs::{File, OpenOptions};
use std::io::{prelude::*, LineWriter};

type CellArray = [Cell; Lattice::CAPACITY];
pub type CytokineArray = [[f64; 2]; Lattice::CAPACITY]; // IFNe: lo,hi
type Neighborhoods = [[usize; Lattice::N_NEIGHBORS]; Lattice::CAPACITY];

#[derive(Clone)]
pub struct Lattice {
    pub neighborhoods: Neighborhoods,
    pub cells: CellArray,
    pub cytokines: CytokineArray,
}

impl Lattice {
    pub const N_NEIGHBORS: usize = 6; // fixed "kissing number" of the lattice, do not change

    pub const WIDTH: usize = 100; // if you increase this, increase also config::THREAD_STACK_SIZE!
    pub const HEIGHT: usize = Lattice::WIDTH; // (non-square lattice shapes are also supported)
    pub const CAPACITY: usize = Lattice::WIDTH * Lattice::HEIGHT;

    pub const OCCUPANCY: f64 = 1.0; // used as ceil(WIDTH * HEIGHT * the given fraction)

    // lattice output
    pub const NEIGHS_TO_FILE: bool = true; // whether lattice neighbor indices are to be dumped
    pub const IMAGE_RESOLUTION: u16 = 100; // default: 100
    pub const IMAGE_RECTANGULAR: bool = true; // if true, the parallelogram-shaped lattice is
                                              // right-to-left wrapped to form a rectangle

    pub fn new(rng: &mut StdRng) -> Self {
        Lattice {
            neighborhoods: Lattice::generate_neighborhods(),
            cells: Lattice::populate_cells(rng),
            cytokines: [[0., 0.]; Lattice::CAPACITY],
        }
    }

    fn generate_neighborhods() -> Neighborhoods {
        let mut nbhoods = [[usize::max_value(); Lattice::N_NEIGHBORS]; Lattice::CAPACITY];
        fn as_index(x: usize, y: usize) -> usize {
            x + y * Lattice::WIDTH
        }
        for (i, nbs) in nbhoods.iter_mut().enumerate() {
            let (x, y) = ((i % Lattice::WIDTH) as usize, (i / Lattice::WIDTH) as usize);
            let (east, west) = (
                (x + 1) % Lattice::WIDTH,
                (x + Lattice::WIDTH - 1) % Lattice::WIDTH,
            );
            let (south, north) = (
                (y + 1) % Lattice::HEIGHT,
                (y + Lattice::HEIGHT - 1) % Lattice::HEIGHT,
            );
            *nbs = [
                as_index(east, y),
                as_index(west, y),
                as_index(x, south),
                as_index(x, north),
                as_index(west, south),
                as_index(east, north),
            ];
        }
        if Lattice::NEIGHS_TO_FILE {
            let nbsf = File::create("neighbors.csv").expect("☠ ☆ neighs");
            let mut nbsf = LineWriter::new(nbsf);
            nbsf.write_all(b"left,right\n").expect("☠ ✏ neighs");
            for (i, nbs) in nbhoods.iter_mut().enumerate() {
                for nbi in nbs.iter() {
                    if nbi > &i {
                        nbsf.write_fmt(format_args!("{:},{:}\n", i, nbi))
                            .expect("☠ ✏ neighs");
                    }
                }
            }
        }
        nbhoods
    }

    fn populate_cells(rng: &mut StdRng) -> CellArray {
        let mut cells = [Cell {
            alive: true,
            molecules: [0; N_MOLECULE_SPECIES],
        }; Lattice::CAPACITY];
        let n_free_nodes = ((1.0 - Lattice::OCCUPANCY) * (cells.len() as f64)) as usize;
        (0..cells.len())
            .collect::<Vec<_>>()
            .choose_multiple(rng, n_free_nodes)
            .for_each(|i| cells[*i].alive = false);
        cells
    }

    fn save_png(&self, time: f64, rates: &Rates) {
        const IMG_SCALING: f64 = 20. * ((Lattice::IMAGE_RESOLUTION as f64) / 100.);
        const R: f64 = IMG_SCALING;
        const H: f64 = IMG_SCALING * 1.732_050 / 2.;
        const X0: f64 = 2. * H;
        const Y0: f64 = 1.5 * R;
        const HEIGHT: f64 = (1.5 * (Lattice::HEIGHT as f64) + 1.5) * R;
        const WIDTH: f64 = (2. * (Lattice::WIDTH as f64)
            + 1.
            + (if Lattice::IMAGE_RECTANGULAR { 2 } else { Lattice::HEIGHT }) as f64)
            * H;

        let sf = ImageSurface::create(Format::Rgb24, WIDTH as i32, HEIGHT as i32).unwrap();
        let cx = Context::new(&sf).unwrap();
        cx.set_source_rgb(0., 0., 0.);
        cx.paint().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));
        cx.set_line_width(0.02 * IMG_SCALING);

        for cell_i in 0..Lattice::CAPACITY {
            // cell index --> its (x, y) coordinates
            let (mut i, j) = (cell_i % Lattice::WIDTH, cell_i / Lattice::WIDTH);
            if Lattice::IMAGE_RECTANGULAR {
                i = (i + j / 2) % Lattice::WIDTH
            }
            let (x, y) = (
                X0 + (2. * (i as f64) + (if Lattice::IMAGE_RECTANGULAR {j % 2} else {j} as f64)) * H,
                Y0 + 1.5 * (j as f64) * R,
            );

            // -- hexagon

            // contour
            cx.move_to(x, y + R * 0.99);
            for a in 2..=6 {
                let z = f64::from(a) * PI / 3.;
                cx.rel_line_to(R * 0.99 * z.sin(), R * 0.99 * z.cos())
            }
            cx.close_path();
            cx.set_source_rgb(0.1, 0.1, 0.1);
            cx.stroke_preserve().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));

            // fill (according to IFNe in lower subcompartment)
            let ifne_lo = self.cytokines[cell_i][0];
            let sat = ifne_lo / (ifne_lo + rates.mm_pstat);
            cx.set_source_rgb(0.15 + 0.85*sat, 0.15 + 0.85*sat, 0.15);
            cx.fill().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));

            // if cell not alive, do not draw rings in hexagon
            if !self.cells[cell_i].alive {
                continue;
            }

            // -- hexagon interior: ring and circle

            // a molecule count
            macro_rules! mlf {
                ($m:ident) => {
                    self.cells[cell_i].molecules[$m as usize] as f64
                };
            }
            // max count of a molecule
            macro_rules! mxf {
                ($m:ident) => {
                    Cell::MAX.molecules[$m as usize] as f64
                };
            }
            // fraction of molecules
            macro_rules! fxn {
                ($m:ident) => {
                    mlf!($m) / mxf!($m)
                };
            }

            // outer ring (according to viral infection progression)
            let infxn = ((if mlf!(Vinf) > 0. { 1. } else { 0. }) + mlf!(Vrna) + mlf!(Vprot))
                / (1. + mxf!(Vrna) + mxf!(Vprot));
            cx.set_source_rgb(0.15 + 0.85*infxn, 0.15, 0.15 + 0.85*infxn);
            cx.arc(x, y, 0.72 * R, 0., 2. * PI);
            cx.fill_preserve().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));
            cx.set_source_rgb(0., 0., 0.);
            cx.stroke().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));

            // inner circle
            cx.set_source_rgb(0.15 + 0.85*fxn!(Pirf3), 0.15 + 0.85*fxn!(Pstat), 0.15);
            cx.arc(x, y, 0.40 * R, 0., 2. * PI);
            cx.fill_preserve().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));
            cx.set_source_rgb(0.15, 0.15, 0.15);
            cx.stroke().unwrap_or_else(|err| println!("☠ ✏ lattice: {:?}", err));
        } // for each cell (lattice node)

        // write out image to a PNG file
        let png_fn = if OUT_FILE_NAME_TIME_IN_MIN {
            ["t_", &format!("{:0>4.0}", time / MIN), "m.png"].concat()
        } else {
            ["t_", &format!("{:0>4.0}", time / HOUR), "h.png"].concat()
        };
        let mut png = File::create(png_fn).expect("☠ ☆ PNG.");
        sf.write_to_png(&mut png).expect("☠ ✏ PNG.");
    }

    fn save_csv(&self, time: f64) {
        // create and open CSV file for writing
        let csv_fn = if OUT_FILE_NAME_TIME_IN_MIN {
            ["t_", &format!("{:0>4.0}", time / MIN), "m.csv"].concat()
        } else {
            ["t_", &format!("{:0>4.0}", time / HOUR), "h.csv"].concat()
        };
        let mut csv = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(csv_fn)
            .expect("☠ ☆ CSV");

        // write out header
        let hdr = "id,alive,Vinf,VRNA,Vprot,pIRF3,IFNi,pSTAT,ISG,IFNeL,IFNeU\n";
        csv.write_all(hdr.as_bytes()).expect("☠ ✏ CSV");

        // write out the state of each cell and the amount of IFNe above the cell
        for cell_i in 0..Lattice::CAPACITY {
            let mut line: Vec<String> = vec![
                cell_i.to_string(),
                (if self.cells[cell_i].alive { "1" } else { "0" }).to_string(),
            ];
            macro_rules! count_s {
                ($m:ident) => {
                    self.cells[cell_i].molecules[$m as usize].to_string()
                };
            }
            for m in vec![
                Mol::Vinf,
                Mol::Vrna,
                Mol::Vprot,
                Mol::Pirf3,
                Mol::Ifni,
                Mol::Pstat,
                Mol::Isg,
            ] {
                line.push(count_s!(m))
            }
            for j in 0..=1 {
                line.push(format!("{:.3e}", self.cytokines[cell_i][j]))
            }
            let mut line_s = line.join(",");
            line_s.push('\n');
            csv.write_all(line_s.as_bytes()).expect("☠ ✏ CSV");
        } // for each cell/lattice node
    }

    // save output file(s)
    pub fn out(&self, time: f64, rates: &Rates, dump_image: bool) {
        if dump_image {
            self.save_png(time, rates);
        }
        self.save_csv(time);
    }
}

#[test]
fn test_lattice_neighborhood_reflectivity() {
    use rand::SeedableRng;
    let mut rng: StdRng = SeedableRng::from_seed([123; 32]);
    let nbhoods = &Lattice::new(&mut rng).neighborhoods;
    for i in 0..nbhoods.len() {
        assert_eq!(nbhoods[i].len(), Lattice::N_NEIGHBORS);
        assert_eq!(nbhoods[ nbhoods[i][0/*E */] ][1/*W */], i);
        assert_eq!(nbhoods[ nbhoods[i][2/*S */] ][3/*N */], i);
        assert_eq!(nbhoods[ nbhoods[i][4/*SW*/] ][5/*NE*/], i);
    }
}
