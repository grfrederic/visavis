from pathlib import Path
from typing import Union
import re

import numpy as np
import pandas as pd

import multiprocessing as mp

from matplotlib.font_manager import FontProperties
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import seaborn as sns

from .simulation_result import SimulationResult



def pairwise_at_hour(result: SimulationResult, hour: Union[int, float], active_x=True, active_y=True):
    ah = result.at_hour(hour)
    ah = ah[ah['alive'] > 0]
    ah = ah[[col for col in ah if "_act" in col]]
    ah = ah.rename(columns={col: col.replace("_act", "") for col in ah})
    ah_x = (ah == active_x).astype(int)
    ah_y = (ah == active_y).astype(int)
    return ah_y.transpose() @ ah_x



def plot_triangle_table(ax, table, bold_x=True, bold_y=True):
    table = table.copy()

    table = table.loc[reversed(table.columns)]
    table.drop(columns=[table.columns[-1]], inplace=True)
    table.drop(index=[table.columns[0]], inplace=True)

    mask = 1 - np.tril(np.ones_like(table))
    mask = mask[::-1]

    #sns.set(font_scale=10)
    mul = 1000  # => permille
    with sns.axes_style("white"):
        ax = sns.heatmap(data=table, mask=mask, annot=True, square=True, cbar=False,
                         linewidths=1, fmt='.0f', cmap=sns.color_palette("BrBG", 22),
                         vmin=-1.2*mul, vmax=1.*mul, ax=ax)
        ax.xaxis.set_ticks_position('top')

    # bold fonts
    fm = FontProperties(weight='bold')
    if bold_x:
        ax.set_xticklabels(ax.get_xticklabels(), fontproperties=fm, rotation=0)
    else:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    if bold_y:
        ax.set_yticklabels(ax.get_yticklabels(), fontproperties=fm, rotation=0)
    else:
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)


def plot_summary(ax, cells, hour: Union[int, float]):
    cells_alive = cells[cells['alive'] > 0]

    n_fields = len(cells)
    n_alive = len(cells_alive)
    n_dead = n_fields - n_alive

    def print_count(ax, what, count, y, kws):
        mul = 1000  # => permille
        ax.annotate(f"{what}:", (0.0, y), **kws)
        ax.annotate(f"{mul * count / n_fields:.0f}", (0.60, y), **kws, ha='right')
        ax.annotate(f"({count})", (0.65, y), **kws, color='gray')

    ax.annotate(f"t = {hour}h", (0.0, 1.0), xycoords='axes fraction', fontsize=18)

    mol_act_max = [
        ('Vinf', 1, 1),
        ('VRNA', 3, 3),
        ('Vprot', 3, 3),
        ('pIRF3', 3, 3),
        ('IFNi', 3, 3),
        ('pSTAT', 1, 3),
        ('ISG', 1, 3),
    ]

    y = 0.95
    y_step = 0.025
    mol_kws = dict(xycoords='axes fraction', fontsize=14)
    for i, (mol, act, max_val) in enumerate(mol_act_max):
        for level in range(max_val + 1):
            kws = {**mol_kws, **(dict(weight='bold') if level >= act else {})}
            mol_on_level = np.sum(cells_alive[mol] == level)
            print_count(ax, f"{mol} {level}", mol_on_level, y, kws)
            y -= y_step
        y -= y_step

    print_count(ax, "fields", n_fields, y, mol_kws)
    y -= y_step
    print_count(ax, "dead", n_dead, y, mol_kws)
    y -= 2*y_step

    ax.annotate("IFNeU/node:", (0.0, y), **mol_kws)
    ax.annotate(f"{int(np.mean(cells['IFNeU']))}", (0.9, y), **mol_kws, ha='right')
    y -= y_step
    ax.annotate("IFNeL/node:", (0.0, y), **mol_kws)
    ax.annotate(f"{int(np.mean(cells['IFNeL']))}", (0.9, y), **mol_kws, ha='right')


_hour_to_image_path_cache = {}
def find_image_path_by_hour(result: SimulationResult, hour: Union[int, float]):
    if hour not in _hour_to_image_path_cache:
        for img_path in result.simulation_dir.glob("*.png"):
            (n, unit), = re.findall(r'(-?\d+)(h|m)\.png', img_path.name)
            h = int(n) / (1 if unit == 'h' else 60)
            _hour_to_image_path_cache[h] = img_path
            
    return _hour_to_image_path_cache[hour]
 

def annotate(result: SimulationResult, hour: Union[int, float], verbose=False):
    cell_img_path = find_image_path_by_hour(result, hour)
    assert cell_img_path.exists(), "no image for this hour"
    cell_img = mpimg.imread(cell_img_path)

    cells = result.at_hour(hour)
    assert len(cells), "no csv for this hour"

    cell_img_height, cell_img_width, _ = cell_img.shape

    if cell_img_width / cell_img_height > 2:
        print("WARNING! image probably already annotated. quitting")
        return

    # 1.1 because merging gives a bit more than 2 lengths
    widths = [1.1 * 2 * cell_img_width / cell_img_height, 0.8, 1, 1]
    heights = [1, 1]

    figsize = (4 * sum(widths), 4 * sum(heights))

    fig = plt.figure(constrained_layout=False, figsize=figsize)


    spec = fig.add_gridspec(ncols=4, nrows=2,
                            width_ratios=widths,
                            height_ratios=heights,
                            left=-0.02, right=0.99,
                            bottom=0.01, top=0.96)

    # cell image
    ax_cell_img = fig.add_subplot(spec[:, 0])
    ax_cell_img.axis('off')
    # offset image bc other plots have labels sticking out over them
    ax_cell_img.imshow(cell_img, aspect='equal', clip_on=False,
                       extent=(0, cell_img_width, 0, cell_img_height))
    ax_cell_img.set_xlim((-0.07 * cell_img_width, 0.93 * cell_img_width))
    ax_cell_img.set_ylim((-0.015 * cell_img_height, 0.985 * cell_img_height))

    # summary columns
    ax_summ = fig.add_subplot(spec[:, 1])
    ax_summ.axis('off')
    plot_summary(ax_summ, cells, hour)

    # triangle tables
    axs_tables = [[fig.add_subplot(spec[0, 2]), fig.add_subplot(spec[0, 3])],
                  [fig.add_subplot(spec[1, 2]), fig.add_subplot(spec[1, 3])]]

    for i, row_tables in enumerate(axs_tables):
        for j, ax in enumerate(row_tables):
            active_x = 1 - i
            active_y = 1 - j

            mul = 1000  # => permille
            pr = mul * pairwise_at_hour(result, hour, active_x=active_x,
                                                      active_y=active_y) \
               / len(cells)  # number of fields

            plot_triangle_table(ax, pr, bold_x=active_x, bold_y=active_y)

    plt.savefig(cell_img_path, dpi=200)

    if verbose:
        print(f"finished annotating: {cell_img_path}")


def annotate_all(result):
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.starmap(
            annotate,
            [(result, hour) for hour in result.hours]
        )
