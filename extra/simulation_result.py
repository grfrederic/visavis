import re

from pathlib import Path

import numpy as np
import pandas as pd


class SimulationResult:
    """
    Wraps directory with simulation results into python class.
    """
    def __init__(self, simulation_dir: Path):
        simulation_dir = Path(simulation_dir)

        self._simulation_dir = simulation_dir

        assert self._simulation_dir.is_dir()

        # load all the csvs
        self._neighbors = pd.read_csv(self._simulation_dir / 'neighbors.csv')

        states_parts = []
        for states_part_path in self._simulation_dir.glob('t_*.csv'):
            (n, unit), = re.findall(r'(-?\d+)(h|m)\.csv', states_part_path.name)
            hour = int(n) / (1 if unit == 'h' else 60)

            states_part = pd.read_csv(states_part_path)
            states_part['hour'] = hour

            states_parts.append(states_part)

        self._states = pd.concat(states_parts)
        self._calculate_acts()

    def _calculate_acts(self):
        self._states['Vinf_act'] = self._states['Vinf'] >= 1
        self._states['VRNA_act'] = self._states['VRNA'] >= 3
        self._states['Vprot_act'] = self._states['Vprot'] >= 1  # 3
        self._states['pIRF3_act'] = self._states['pIRF3'] >= 3
        self._states['IFNi_act'] = self._states['IFNi'] >= 3
        self._states['pSTAT_act'] = self._states['pSTAT'] >= 1
        self._states['ISG_act'] = self._states['ISG'] >= 3

    @property
    def simulation_dir(self):
        return self._simulation_dir

    @property
    def neighbors(self):
        return self._neighbors

    @property
    def states(self):
        return self._states

    @property
    def hours(self):
        return sorted(set(self._states['hour']))

    def at_hour(self, hour: int):
        return self.states[self.states['hour'] == hour]

    def ngh_at_hour(self, hour: int):
        ah = self.at_hour(hour)
        return (
            self.neighbors
            .merge(ah, left_on='left', right_on='id')
            .merge(ah, left_on='right', right_on='id', suffixes=('_left', '_right'))
            .drop(['left', 'right'], axis=1)
        )

    def ngh_ks_at_hour(self, col1: str, col2: str, hour: int):
        sn = self.ngh_at_hour(hour)

        c1l, c1r = f"{col1}_left", f"{col1}_right"
        c2l, c2r = f"{col2}_left", f"{col2}_right"

        qaa = ( sn[c1l] &  sn[c2r]).sum() + ( sn[c1r] &  sn[c2l]).sum()
        qan = ( sn[c1l] & ~sn[c2r]).sum() + ( sn[c1r] & ~sn[c2l]).sum()
        qna = (~sn[c1l] &  sn[c2r]).sum() + (~sn[c1r] &  sn[c2l]).sum()
        qnn = (~sn[c1l] & ~sn[c2r]).sum() + (~sn[c1r] & ~sn[c2l]).sum()

        qa = qaa + qan
        qn = qna + qnn

        return qaa / qa + qnn / qn - 1 if qa * qn != 0 else np.nan

    def ks_at_hour(self, col1: str, col2: str, hour: int):
        state = self.at_hour(hour)
        state_col1 = state[col1]
        state_col2 = state[col2]

        qaa = ( state_col1 &  state_col2).sum()
        qan = ( state_col1 & ~state_col2).sum()
        qna = (~state_col1 &  state_col2).sum()
        qnn = (~state_col1 & ~state_col2).sum()

        qa = qaa + qan
        qn = qna + qnn

        return qaa / qa + qnn / qn - 1 if qa * qn != 0 else np.nan
