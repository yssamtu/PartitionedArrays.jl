from collections import defaultdict
from json import load
from numpy import prod
from os import listdir
from os.path import isdir, join
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union

from common import change_func_name, draw_legend
from constants import (
    SPECIAL_LINE_COLOUR,
    SPECIAL_LINE_WIDTH,
    SPECIAL_LINE_STYLE,
    FOLDER,
)


class StrongData:
    def __init__(self, computer: str = "das5"):
        data = self.__load_strong_scaling_data(computer)
        one_case = next(iter(next(iter(next(iter(data.values())).values()))))
        self.data = data
        self.data_sizes = tuple(data.keys())
        self.nps = tuple(next(iter(data.values())).keys())
        self.fs = tuple(k for k, v in one_case.items() if isinstance(v, dict))
        self.times = tuple(next(iter(one_case.values())).keys())
        self.__info_key = tuple(
            k for k, v in one_case.items() if not isinstance(v, dict)
        )
        self.__speedup = self.get_speedup()
        self.__efficiency = self.get_efficiency()

    def __load_strong_scaling_data(self, computer: str):
        """__load_strong_scaling_data(computer: str)"""

        def convert_to_dict(d):
            def f(d):
                if isinstance(d, defaultdict):
                    d = {k: convert_to_dict(v) for k, v in sorted(d.items())}
                return d

            if isinstance(d, defaultdict):
                d = {k: f(v) for k, v in reversed(sorted(d.items()))}
            return d

        if computer == "das5":
            folder = join(FOLDER, "strong_scaling")
        else:
            folder = join(FOLDER, "strong_scaling_old_snellius")
        file_name = "summary.json"
        dicts = defaultdict(lambda: defaultdict(list))
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
                size = int(size_str.strip("()").split(",")[0])
                part = tuple(map(int, part_str.strip("()").split(",")))
                with open(join(full_path, file_name)) as f:
                    dicts[size][np].append(
                        {**change_func_name(load(f)), "node_core": node_core, "part": part}
                    )
        data = convert_to_dict(dicts)
        for v in data.values():
            for k, v_1 in v.items():
                v[k] = tuple(v_1)
        return data

    def get_data(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_data(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
        if isinstance(data_sizes, int):
            data_sizes = (data_sizes,)
        elif data_sizes is None or not data_sizes:
            data_sizes = self.data_sizes
        elif isinstance(data_sizes, list):
            data_sizes = tuple(data_sizes)
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
        if (
            data_sizes == self.data_sizes
            and nps == self.nps
            and fs == self.fs
            and times == self.times
        ):
            return self.data
        fs = (*fs, *self.__info_key)
        return {
            data_size: {
                np: tuple(
                    {
                        f: (
                            {time: case[f][time] for time in times}
                            if f in self.fs
                            else case[f]
                        )
                        for f in fs
                    }
                    for case in self.data[data_size][np]
                )
                for np in nps
            }
            for data_size in data_sizes
        }

    def get_gathered_data(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_gathered_data(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
        if isinstance(nps, int):
            nps = (nps,)
        elif nps is None or not nps:
            nps = self.nps
        elif isinstance(nps, list):
            nps = tuple(nps)
        data = self.get_data(data_sizes=data_sizes, nps=nps, fs=fs, times=times)
        return self.__to_StrongGatheredData(nps=nps, data=data)

    def get_speedup(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_speedup(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
        if isinstance(data_sizes, int):
            data_sizes = (data_sizes,)
        elif data_sizes is None or not data_sizes:
            data_sizes = self.data_sizes
        elif isinstance(data_sizes, list):
            data_sizes = tuple(data_sizes)
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
        if (
            data_sizes == self.data_sizes
            and nps == self.nps
            and fs == self.fs
            and times == self.times
            and hasattr(self, "_StrongData__speedup")
        ):
            return self.__speedup
        data = self.get_gathered_data(
            data_sizes=data_sizes, nps=nps, fs=fs, times=times
        )
        base_data = (
            data
            if 1 in nps
            else self.get_gathered_data(
                data_sizes=data_sizes, nps=1, fs=fs, times=times
            )
        )
        for case, base_case in zip(data.values(), base_data.values()):
            for v, base_v in zip(case.values(), base_case.values()):
                for (k_1, v_1), base_v_1 in zip(v.items(), base_v.values()):
                    base = base_v_1[0]
                    v[k_1] = tuple(base / value if value else 0.0 for value in v_1)
        return data

    def get_efficiency(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_efficiency(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
        if isinstance(data_sizes, int):
            data_sizes = (data_sizes,)
        elif data_sizes is None or not data_sizes:
            data_sizes = self.data_sizes
        elif isinstance(data_sizes, list):
            data_sizes = tuple(data_sizes)
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
        if (
            data_sizes == self.data_sizes
            and nps == self.nps
            and fs == self.fs
            and times == self.times
            and hasattr(self, "_StrongData__efficiency")
        ):
            return self.__efficiency
        data = self.get_gathered_data(
            data_sizes=data_sizes, nps=nps, fs=fs, times=times
        )
        base_data = (
            data
            if 1 in nps
            else self.get_gathered_data(
                data_sizes=data_sizes, nps=1, fs=fs, times=times
            )
        )
        for case, base_case in zip(data.values(), base_data.values()):
            for v, base_v in zip(case.values(), base_case.values()):
                for (k_1, v_1), base_v_1 in zip(v.items(), base_v.values()):
                    base = base_v_1[0]
                    v[k_1] = tuple(
                        base / value / np if value else 0.0
                        for value, np in zip(v_1, nps)
                    )
        return data

    def __to_StrongGatheredData(self, nps: tuple[int, ...] = None, data=None):
        if data is None:
            data = self.data
        result = {}
        for data_size, v_1 in data.items():
            one_case = next(iter(v_1.values()))[0]
            result[data_size] = {
                f: {time: [float("inf")] * len(nps) for time in one_case[f].keys()}
                for f in one_case.keys()
                if isinstance(one_case[f], dict)
            }
            for i, np in enumerate(nps):
                for d in v_1[np]:
                    for f, time_data in result[data_size].items():
                        for time, v_3 in d[f].items():
                            if v_3 < time_data[time][i]:
                                time_data[time][i] = v_3
        for v in result.values():
            for v_1 in v.values():
                for k, v_3 in v_1.items():
                    v_1[k] = tuple(v_3)
        return result

    def draw_data(
        self,
        draw_type: str = "lines",
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        draw_limit: bool = True,
        figsize: tuple[int, int] = (10, 10),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_data(draw_type: str, data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], draw_limit: bool, figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plots(ax, time, show_ylabel=True):
            for data_size, marker in zip(data_sizes, MARKERS):
                for i, f in enumerate(fs):
                    y = [float("inf")] * len(nps)
                    for j, np in enumerate(nps):
                        y[j] = min(
                            y[j],
                            *(point[f][time] for point in self.data[data_size][np]),
                        )
                    ax.plot(
                        x, y, color=colors[i], marker=marker, label=f"{f}_{data_size}"#, alpha=0.5
                    )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Wall-clock time (s)")

        def scatters(ax, time, show_ylabel=True):
            for data_size, marker in zip(data_sizes, MARKERS):
                for i, f in enumerate(fs):
                    y = []
                    for j, np in enumerate(nps):
                        y.extend((point[f][time] for point in self.data[data_size][np]))
                    ax.scatter(
                        x, y, color=colors[i], marker=marker, label=f"{f}_{data_size}"
                    )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Wall-clock time (s)")

        colors = sns.color_palette("Set1", 9, 0.9)
        MARKERS = ("x", "+", "v", ".", "*")
        if isinstance(data_sizes, int):
            data_sizes = (data_sizes,)
        elif data_sizes is None or not data_sizes:
            data_sizes = self.data_sizes
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
            palette=colors,
            rc={
                "text.usetex": plt.rcParams["text.usetex"],
                "font.family": plt.rcParams["font.family"],
                "font.serif": plt.rcParams["font.serif"],
            },
        )
        if draw_type == "lines":
            x = nps
            func = plots
        else:
            x = tuple(np for np in nps for _ in self.data[data_sizes[0]][np])
            func = scatters
        if isinstance(times, tuple):
            fig, axs = plt.subplots(
                1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=True
            )
            ylabel_mask = (True, False)
            if draw_limit and len(data_sizes) == 1:
                if save_file_name is None or not separate_legend:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        ax.axvline(
                            x=pow(data_sizes[0] - 1, 3) / 25000,
                            color=SPECIAL_LINE_COLOUR,
                            linewidth=SPECIAL_LINE_WIDTH,
                            linestyle=SPECIAL_LINE_STYLE,
                            label="25K rows/rank",
                        )
                        func(ax, time, show_ylabel)
                        ax.legend()
                else:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        ax.axvline(
                            x=pow(data_sizes[0] - 1, 3) / 25000,
                            color=SPECIAL_LINE_COLOUR,
                            linewidth=SPECIAL_LINE_WIDTH,
                            linestyle=SPECIAL_LINE_STYLE,
                            label="25K rows/rank",
                        )
                        func(ax, time, show_ylabel)
            else:
                if save_file_name is None or not separate_legend:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        func(ax, time, show_ylabel)
                        ax.legend()
                else:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        func(ax, time, show_ylabel)
            if save_file_name is None:
                plt.show()
            else:
                if separate_legend:
                    draw_legend(axs[0], figsize, save_file_name)
                plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
        else:
            fig, ax = plt.subplots(figsize=figsize)
            if draw_limit and len(data_sizes) == 1:
                ax.axvline(
                    x=pow(data_sizes[0] - 1, 3) / 25000,
                    color=SPECIAL_LINE_COLOUR,
                    linewidth=SPECIAL_LINE_WIDTH,
                    linestyle=SPECIAL_LINE_STYLE,
                    label="25K rows/rank",
                )
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

    def draw_complexity(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (6, 6),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_complexity(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""
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
            style="ticks",
            palette=sns.color_palette("Set1", 9, 0.9),
            rc={
                "text.usetex": plt.rcParams["text.usetex"],
                "font.family": plt.rcParams["font.family"],
                "font.serif": plt.rcParams["font.serif"],
            },
        )
        MARKER = "x"
        x = tuple(k**3 for k in self.data.keys())

        def plot(ax, f, time, show_y_label=True):
            for np in nps:
                y = []
                for data_size in self.data.keys():
                    y.append(min(case[f][time] for case in self.data[data_size][np]))
                ax.plot(
                    x,
                    tuple(i / y[-1] * y[-1] for i in y),
                    marker=MARKER,
                    label=f"{f}_{np}",
                )
            ax.plot(
                x, tuple((i / x[-1]) * y[-1] for i in x), marker=MARKER, label=r"$O(n)$"
            )
            ax.plot(
                x,
                tuple((i / x[-1]) ** 0.5 * y[-1] for i in x),
                marker=MARKER,
                label=r"$O(\sqrt n)$",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("Problem size")
            if show_y_label:
                ax.set_ylabel("Wall-clock time (s)")

        if isinstance(times, tuple):
            figsize = (figsize[0] * 2 * 1.1, figsize[1])
            if save_file_name is None or not separate_legend:
                for f in fs:
                    fig, axs = plt.subplots(1, 2, figsize=figsize)
                    for i, (ax, time) in enumerate(zip(axs, times)):
                        plot(ax, f, time, i == 0)
                        ax.legend()
                    if save_file_name is None:
                        plt.show()
                    else:
                        plt.savefig(f"{save_file_name}_{f}.pdf", bbox_inches="tight")
                    plt.close(fig)
            else:
                for f in fs:
                    fig, axs = plt.subplots(1, 2, figsize=figsize)
                    for i, (ax, time) in enumerate(zip(axs, times)):
                        plot(ax, f, time, i == 0)
                    draw_legend(axs[0], figsize, f"{save_file_name}_{f}")
                    plt.savefig(f"{save_file_name}_{f}.pdf", bbox_inches="tight")
                    plt.close(fig)
        else:
            for f in fs:
                fig, ax = plt.subplots(figsize=figsize)
                plot(ax, f, times)
                if save_file_name is None:
                    ax.legend()
                    plt.show()
                else:
                    if separate_legend:
                        draw_legend(ax, figsize, f"{save_file_name}_{f}")
                    plt.savefig(f"{save_file_name}_{f}.pdf", bbox_inches="tight")
                plt.close(fig)

    def draw_speedup(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (8, 8),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_speedup(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plot(data_size, ax, time, show_ylabel=True):
            y = nps
            ax.plot(
                nps,
                y,
                color=SPECIAL_LINE_COLOUR,
                linewidth=SPECIAL_LINE_WIDTH,
                label="theoretical",
                linestyle=SPECIAL_LINE_STYLE,
            )
            for f in fs:
                y = speedup[data_size][f][time]
                ax.plot(nps, y, marker=MARKER, label=f"{f}_{data_size}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"CPU cores")
            if show_ylabel:
                ax.set_ylabel(r"Speed-up")

        if isinstance(data_sizes, int):
            data_sizes = (data_sizes,)
        elif data_sizes is None or not data_sizes:
            data_sizes = self.data_sizes
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
        speedup = self.get_speedup(data_sizes=data_sizes, nps=nps, fs=fs, times=times)
        sns.set_theme(
            context="paper",
            style="whitegrid",
            palette=sns.color_palette("Set1", 9, 0.9),
            rc={
                "text.usetex": plt.rcParams["text.usetex"],
                "font.family": plt.rcParams["font.family"],
                "font.serif": plt.rcParams["font.serif"],
            },
        )
        MARKER = "x"
        if isinstance(times, tuple):
            for data_size in data_sizes:
                fig, axs = plt.subplots(
                    1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=True
                )
                if save_file_name is None or not separate_legend:
                    for i, (ax, time) in enumerate(zip(axs, times)):
                        plot(data_size, ax, time, i == 0)
                        ax.legend()
                else:
                    for i, (ax, time) in enumerate(zip(axs, times)):
                        plot(data_size, ax, time, i == 0)
                if save_file_name is None:
                    plt.show()
                else:
                    if separate_legend:
                        draw_legend(axs[0], figsize, save_file_name)
                    plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
                plt.close(fig)
        else:
            for data_size in data_sizes:
                fig, ax = plt.subplots(figsize=figsize)
                plot(data_size, ax, times)
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

    def draw_efficiency(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (8, 8),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_efficiency(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plot(data_size, ax, time, show_ylabel=True):
            y = [1.0] * len(nps)
            ax.plot(nps, y, label="theoretical", linestyle=SPECIAL_LINE_STYLE)
            for f in fs:
                y = efficiency[data_size][f][time]
                ax.plot(nps, y, marker=MARKER, label=f"{f}_{data_size}")
            ax.set_xscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Efficiency")

        if isinstance(data_sizes, int):
            data_sizes = (data_sizes,)
        elif data_sizes is None or not data_sizes:
            data_sizes = self.data_sizes
        if isinstance(nps, int):
            nps = (nps,)
        elif nps is None or not data_sizes:
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
        efficiency = self.get_efficiency(
            data_sizes=data_sizes, nps=nps, fs=fs, times=times
        )
        sns.set_theme(
            context="paper",
            style="whitegrid",
            palette=sns.color_palette("Set1", 9, 0.9),
            rc={
                "text.usetex": plt.rcParams["text.usetex"],
                "font.family": plt.rcParams["font.family"],
                "font.serif": plt.rcParams["font.serif"],
            },
        )
        MARKER = "x"
        if isinstance(times, tuple):
            for data_size in data_sizes:
                fig, axs = plt.subplots(
                    1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=True
                )
                if save_file_name is None or not separate_legend:
                    for i, (ax, time) in enumerate(zip(axs, times)):
                        plot(data_size, ax, time, i == 0)
                        ax.legend()
                else:
                    for i, (ax, time) in enumerate(zip(axs, times)):
                        plot(data_size, ax, time, i == 0)
                if save_file_name is None:
                    plt.show()
                else:
                    if separate_legend:
                        draw_legend(axs[0], figsize, save_file_name)
                    plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
                plt.close(fig)
        else:
            for data_size in data_sizes:
                fig, ax = plt.subplots(figsize=figsize)
                plot(data_size, ax, times)
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
