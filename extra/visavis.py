# python interface for running simulations

import re
import string
import random
import shutil
import subprocess

from importlib import import_module

from pathlib import Path

from typing import Any, Optional

import termcolor
import json

from .simulation_result import SimulationResult
from .annotate import annotate_all



def _random_name(k: int) -> str:
    return ''.join(random.choices(string.ascii_uppercase, k=k))


def read_parameters_from_py(module_path: str) -> Any:
    return import_module(module_path).parameters


class VisAVisClient:
    """
    Wraps running simulations into python class.
    """
    def __init__(
        self,
        infection_bin: Path = Path('infection'),
        sim_root: Path = Path('/tmp'),  # where simulation dirs go
    ):
        if isinstance(infection_bin, str):
            infection_bin = Path(infection_bin)

        if isinstance(sim_root, str):
            sim_root = Path(sim_root)

        self._infection_bin = infection_bin
        if not self._infection_bin.is_file():
            self._infection_bin = self._infection_bin.with_suffix('.exe')

        assert self._infection_bin.is_file()

        self._sim_root = sim_root
        self._sim_root.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        parameters_json: Any,            # ...thing that serializes to json
        protocol_file_path: Path,
        dir_name: Optional[str] = None,  # name the dir with results
        clean_up: bool = True,           # remove files?
        verbose: bool = True,            # print information about progress
        images:   bool = False,          # save output images
        annotate: bool = False,          # annotate output images?
    ) -> Optional[SimulationResult]:

        if isinstance(protocol_file_path, str):
            protocol_file_path = Path(protocol_file_path)

        if annotate:
            assert images

        if dir_name is None:
            simulation_dir = self._sim_root / f'sim_{_random_name(16)}'
        else:
            simulation_dir = self._sim_root / dir_name
        simulation_dir.mkdir()

        parameters = simulation_dir / 'parameters.json'
        with open(parameters, 'w') as f:
            json.dump(parameters_json, f)

        assert protocol_file_path.exists()
        protocol_file_src_path = protocol_file_path
        protocol_file_dst_path = simulation_dir / protocol_file_src_path.name
        shutil.copy(protocol_file_src_path.absolute(),
                    protocol_file_dst_path.absolute())

        if verbose:
            termcolor.cprint(f'Starting simulation {simulation_dir}\n', end='', color='green')

        ret = subprocess.call([
            self._infection_bin.absolute(),
            parameters.absolute(),
            protocol_file_dst_path.absolute(),
            *(['--images'] if images else []),
        ], cwd=simulation_dir, stdout=subprocess.DEVNULL)


        if ret:
            if verbose:
                termcolor.cprint(f"Simulation {simulation_dir} failed\n", end='', color='red')
            if clean_up:
                shutil.rmtree(simulation_dir)
            return None
        else:
            if verbose:
                termcolor.cprint(f'Finished simulation {simulation_dir}\n', end='', color='green')

        res = SimulationResult(simulation_dir)

        if clean_up:
            shutil.rmtree(simulation_dir)
        elif annotate:
            if verbose:
                termcolor.cprint(f'Starting annotation {simulation_dir}\n', end='', color='green')
            annotate_all(res)
            if verbose:
                termcolor.cprint(f'Finished annotation {simulation_dir}\n', end='', color='green')

        return res
