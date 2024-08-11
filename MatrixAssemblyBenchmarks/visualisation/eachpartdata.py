from collections import defaultdict
from json import load
import matplotlib.pyplot as plt
from numpy import array, prod, zeros
from os import listdir
from os.path import isdir, join
import seaborn as sns
from typing import Union

from common import change_func_name, draw_legend
from constants import FOLDER


class EachPartData:
    def __init__(self):
        data = self.__load_each_part_data()
        one_case = next(iter(next(iter(data.values()))))
        self.data = data
        self.nps = tuple(data.keys())
        self.fs = tuple(k for k, v in one_case.items() if isinstance(v, dict))
        self.times = tuple(next(iter(one_case.values())).keys())
        self.stages = {
            key: tuple(value) for key, value in next(iter(one_case.values())).items()
        }
        self.__info_key = tuple(
            k for k, v in one_case.items() if not isinstance(v, dict)
        )

    def __load_each_part_data(self):
        """__load_each_part_data()"""

        def set_stage_name(data):
            build_mappings = {
                0: "find_owner_I",
                1: "union_ghost_rows",
                2: "assembly_neighbors",
                3: "partition_and_setup_cache_snd!",
                4: "ExchangeGraph",
                5: "exchange_I",
                6: "exchange_J",
                7: "exchange_V",
                8: "fetch_I",
                9: "fetch_J",
                10: "fetch_V",
                11: "store_recv_data!",
                12: "set_rows_fa",
                13: "find_owner_J",
                14: "union_ghost_cols",
                15: "split_and_compress!",
                16: "pack_cache",
                17: "set_flag_assembled",
                18: "PSparseMatrix",
            }
            rebuild_mappings = {
                0: "unpack_cache",
                1: "partition_and_setup_cache_snd!",
                2: "exchange!",
                3: "fetch",
                4: "store_recv_data!",
                5: "split_and_compress!",
            }
            return {
                f: {
                    time: {
                        (build_mappings if time == "build_time" else rebuild_mappings)[
                            i
                        ]: v_3
                        for i, v_3 in enumerate(v_2)
                    }
                    for time, v_2 in v_1.items()
                }
                for f, v_1 in data.items()
            }

        def convert_to_dict(d):
            if isinstance(d, defaultdict):
                d = {k: convert_to_dict(v) for k, v in sorted(d.items())}
            return d

        folder = join(FOLDER, "each_part")
        file_name = "summary.json"
        dicts = defaultdict(list)
        for nc in listdir(folder):
            path = join(folder, nc)
            if not isdir(path):
                continue
            node_core = tuple(map(int, nc.strip("()").split(",")))
            np = int(prod(node_core))
            for subpath in listdir(path):
                full_path = join(path, subpath)
                if not isdir(full_path):
                    continue
                size_str, part_str, _ = subpath.split("_")
                size = tuple(map(int, size_str.strip("()").split(",")))
                part = tuple(map(int, part_str.strip("()").split(",")))
                with open(join(full_path, file_name)) as f:
                    dicts[np].append(
                        {
                            **set_stage_name(change_func_name(load(f))),
                            "node_core": node_core,
                            "part": part,
                            "data_size": size,
                        }
                    )
        data = convert_to_dict(dicts)
        for k, v in data.items():
            data[k] = tuple(v)
        return data

    def get_data(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_data(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
        if isinstance(nps, int):
            nps = (nps,)
        elif nps is None or not nps:
            nps = self.nps
        elif isinstance(nps, list):
            nps = tuple(nps)
        if isinstance(fs, str):
            fs = (fs,)
        elif fs is None or not fs:
            fs = self.fs
        elif isinstance(fs, list):
            fs = tuple(fs)
        if isinstance(times, str):
            times = (times,)
        elif times is None or not times:
            times = self.times
        elif isinstance(times, list):
            times = tuple(times)
        if nps == self.nps and fs == self.fs and times == self.times:
            return self.data
        fs = (*fs, *self.__info_key)
        return {
            np: tuple(
                {
                    f: (
                        {time: case[f][time] for time in times}
                        if f in self.fs
                        else case[f]
                    )
                    for f in fs
                }
                for case in self.data[np]
            )
            for np in nps
        }

    def get_plot_data_stages_colors(self, mode: str = "each_part"):
        """get_plot_data_stages_colors(mode: str)"""
        if mode == "each_part":
            build_colors = sns.color_palette("Set1", 9, 0.9) + sns.color_palette(
                "Set3", 10, 0.9
            )
            rebuild_colors = [
                build_colors[16],
                build_colors[3],
                build_colors[7],
                build_colors[10],
                build_colors[11],
                build_colors[15],
            ]
            return (
                self.data,
                self.stages,
                {
                    self.times[0]: build_colors,
                    self.times[1]: rebuild_colors,
                },
            )
        elif mode == "each_stage":
            mappings = {
                self.times[0]: {
                    "find_owner_I": "auxiliary",
                    "union_ghost_rows": "auxiliary",
                    "assembly_neighbors": "exchange",
                    "partition_and_setup_cache_snd!": "send_data_preparation",
                    "ExchangeGraph": "auxiliary",
                    "exchange_I": "exchange",
                    "exchange_J": "exchange",
                    "exchange_V": "exchange",
                    "fetch_I": "exchange",
                    "fetch_J": "exchange",
                    "fetch_V": "exchange",
                    "store_recv_data!": "data_consolidation",
                    "set_rows_fa": "auxiliary",
                    "find_owner_J": "auxiliary",
                    "union_ghost_cols": "auxiliary",
                    "split_and_compress!": "data_filing",
                    "pack_cache": "auxiliary",
                    "set_flag_assembled": "auxiliary",
                    "PSparseMatrix": "auxiliary",
                },
                self.times[1]: {
                    "unpack_cache": "auxiliary",
                    "partition_and_setup_cache_snd!": "send_data_preparation",
                    "exchange!": "exchange",
                    "fetch": "exchange",
                    "store_recv_data!": "data_consolidation",
                    "split_and_compress!": "data_filing",
                },
            }
            stages = {
                time: (
                    "auxiliary",
                    "send_data_preparation",
                    "exchange",
                    "data_consolidation",
                    "data_filing",
                )
                for time in self.times
            }
            colors = {time: sns.color_palette("Set1", 5, 0.9) for time in self.times}
            data = {
                np: tuple(
                    {
                        f: (
                            {
                                time: {stage: 0.0 for stage in stages[time]}
                                for time in v_3.keys()
                            }
                            if isinstance(v_3, dict)
                            else v_3
                        )
                        for f, v_3 in v_2.items()
                    }
                    for v_2 in v_1
                )
                for np, v_1 in self.data.items()
            }
            for new_v_1, v_1 in zip(data.values(), self.data.values()):
                for new_v_2, v_2 in zip(new_v_1, v_1):
                    for new_v_3, v_3 in zip(new_v_2.values(), v_2.values()):
                        if not isinstance(new_v_3, dict):
                            continue
                        for new_v_4, (time, v_4) in zip(new_v_3.values(), v_3.items()):
                            for stage, value in v_4.items():
                                new_v_4[mappings[time][stage]] += value
            return data, stages, colors
        else:
            mappings = {
                self.times[0]: {
                    "find_owner_I": "computation",
                    "union_ghost_rows": "computation",
                    "assembly_neighbors": "communication",
                    "partition_and_setup_cache_snd!": "computation",
                    "ExchangeGraph": "computation",
                    "exchange_I": "communication",
                    "exchange_J": "communication",
                    "exchange_V": "communication",
                    "fetch_I": "communication",
                    "fetch_J": "communication",
                    "fetch_V": "communication",
                    "store_recv_data!": "computation",
                    "set_rows_fa": "computation",
                    "find_owner_J": "computation",
                    "union_ghost_cols": "computation",
                    "split_and_compress!": "computation",
                    "pack_cache": "computation",
                    "set_flag_assembled": "computation",
                    "PSparseMatrix": "computation",
                },
                self.times[1]: {
                    "unpack_cache": "computation",
                    "partition_and_setup_cache_snd!": "computation",
                    "exchange!": "communication",
                    "fetch": "communication",
                    "store_recv_data!": "computation",
                    "split_and_compress!": "computation",
                },
            }
            stages = {time: ("computation", "communication") for time in self.times}
            colors = {time: sns.color_palette("Set1", 2, 0.9) for time in self.times}
            data = {
                np: tuple(
                    {
                        f: (
                            {
                                time: {stage: 0.0 for stage in stages[time]}
                                for time in v_3.keys()
                            }
                            if isinstance(v_3, dict)
                            else v_3
                        )
                        for f, v_3 in v_2.items()
                    }
                    for v_2 in v_1
                )
                for np, v_1 in self.data.items()
            }
            for new_v_1, v_1 in zip(data.values(), self.data.values()):
                for new_v_2, v_2 in zip(new_v_1, v_1):
                    for new_v_3, v_3 in zip(new_v_2.values(), v_2.values()):
                        if not isinstance(new_v_3, dict):
                            continue
                        for new_v_4, (time, v_4) in zip(new_v_3.values(), v_3.items()):
                            for stage, value in v_4.items():
                                new_v_4[mappings[time][stage]] += value
            return data, stages, colors

    def draw_data(
        self,
        mode: str = "each_part",
        draw_type: str = "lines",
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (10, 10),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_data(mode: str, draw_type: str, nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plots(ax, time, show_ylabel=True):
            MARKER = "x"
            data, stages, colors = self.get_plot_data_stages_colors(mode)
            for f in fs:
                y = [data[np][0][f][time] for np in nps]
                for i, np in enumerate(nps):
                    for point in data[np]:
                        if sum(point[f][time].values()) < sum(y[i].values()):
                            y[i] = point[f][time]
                y_sum = zeros(len(nps))
                # total_times = array([sum(y_np.values()) for y_np in y])
                for i, stage in enumerate(stages[time]):
                    y_stage = array([y_np[stage] for y_np in y])
                    # y_stage = array([y_np[stage] for y_np in y]) / total_times
                    prev_y_sum = y_sum.copy()
                    y_sum += y_stage
                    ax.plot(x, y_sum, color=colors[time][i], marker="", label=stage)
                    ax.fill_between(
                        x, y_sum, prev_y_sum, color=colors[time][i], alpha=0.5
                    )
            ax.set_xscale("log")
            # ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Wall-clock time (s)")

        def scatters(ax, time, show_ylabel=True):
            data, stages, colors = self.get_plot_data_stages_colors(mode)
            for f in fs:
                y = [data[np][0][f][time] for np in nps]
                for i, np in enumerate(nps):
                    for point in data[np]:
                        if sum(point[f][time].values()) < sum(y[i].values()):
                            y[i] = point[f][time]
                log_width = 0.03
                left_edges = array(x) * 10 ** (-log_width / 2)
                right_edges = array(x) * 10 ** (log_width / 2)
                widths = right_edges - left_edges
                prev_y_sum = zeros(len(nps))
                total_times = array([sum(y_np.values()) for y_np in y])
                for i, stage in enumerate(stages[time]):
                    y_stage = array([y_np[stage] for y_np in y]) * 100.0 / total_times
                    ax.bar(
                        x,
                        y_stage,
                        color=colors[time][i],
                        label=stage,
                        bottom=prev_y_sum,
                        width=widths,
                    )
                    prev_y_sum += y_stage
            ax.set_xscale("log")
            # ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel(r"Percentage of total wall-clock time (\%)")

        if isinstance(nps, int):
            nps = (nps,)
        elif nps is None or not nps:
            nps = self.nps
        if isinstance(fs, str):
            fs = (fs,)
        elif fs is None or not fs:
            fs = self.fs
        if times is None or not times:
            times = self.times
        elif not isinstance(times, str) and len(times) == 1:
            times = times[0]
        elif isinstance(times, list):
            times = tuple(times)
        sns.set_theme(
            context="paper",
            style="whitegrid",
            rc={
                "text.usetex": plt.rcParams["text.usetex"],
                "font.family": plt.rcParams["font.family"],
                "font.serif": plt.rcParams["font.serif"],
            },
        )
        x = nps
        func = plots if draw_type == "lines" else scatters
        if isinstance(times, tuple):
            fig, axs = plt.subplots(
                1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=False
            )
            if save_file_name is None or not separate_legend:
                for i, (ax, time) in enumerate(zip(axs, times)):
                    func(ax, time, i == 0)
                    ax.legend()
            else:
                for i, (ax, time) in enumerate(zip(axs, times)):
                    func(ax, time, i == 0)
            if save_file_name is None:
                plt.show()
            else:
                if separate_legend:
                    draw_legend(axs[0], figsize, save_file_name)
                plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
        else:
            fig, ax = plt.subplots(figsize=figsize)
            func(ax, times)
            if save_file_name is None:
                ax.legend()
                plt.show()
            else:
                if separate_legend:
                    draw_legend(ax, figsize, save_file_name)
                else:
                    ax.legend()
                plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
        plt.close(fig)
