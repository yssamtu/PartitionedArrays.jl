from collections import defaultdict
from json import load
import matplotlib.pyplot as plt
from numpy import prod
from os import listdir
from os.path import isdir, join
import seaborn as sns
from typing import Union

from common import change_func_name, draw_legend
from constants import (
    FOLDER,
    SPECIAL_LINE_COLOUR,
    SPECIAL_LINE_STYLE,
    SPECIAL_LINE_WIDTH,
)


class WeakData:
    def __init__(self):
        data = self.__load_weak_scaling_data()
        one_case = next(iter(next(iter(data.values()))))
        self.data = data
        self.nps = tuple(data.keys())
        self.fs = tuple(k for k, v in one_case.items() if isinstance(v, dict))
        self.times = tuple(next(iter(one_case.values())).keys())
        self.__info_key = tuple(
            k for k, v in one_case.items() if not isinstance(v, dict)
        )
        self.__speedup = self.get_speedup()
        self.__efficiency = self.get_efficiency()

    def __load_weak_scaling_data(self):
        """load_weak_scaling_data()"""

        def convert_to_dict(d):
            if isinstance(d, defaultdict):
                d = {k: convert_to_dict(v) for k, v in sorted(d.items())}
            return d

        folder = join(FOLDER, "weak_scaling")
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
                            **change_func_name(load(f)),
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

    def get_gathered_data(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_gathered_data(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
        if isinstance(nps, int):
            nps = (nps,)
        elif nps is None or not nps:
            nps = self.nps
        elif isinstance(nps, list):
            nps = tuple(nps)
        data = self.get_data(nps=nps, fs=fs, times=times)
        return self.__to_WeakGatheredData(nps=nps, data=data)

    def get_speedup(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_speedup(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
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
            nps == self.nps
            and fs == self.fs
            and times == self.times
            and hasattr(self, "_WeakData__speedup")
        ):
            return self.__speedup
        data = self.get_gathered_data(nps=nps, fs=fs, times=times)
        base_data = (
            data if 1 in nps else self.get_gathered_data(nps=1, fs=fs, times=times)
        )
        for case, base_case in zip(data.values(), base_data.values()):
            for (k, v), base_v in zip(case.items(), base_case.values()):
                base = base_v[0]
                case[k] = tuple(
                    base / value * np if value else 0.0 for value, np in zip(v, nps)
                )
        return data

    def get_efficiency(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
    ):
        """get_efficiency(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]])"""
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
            nps == self.nps
            and fs == self.fs
            and times == self.times
            and hasattr(self, "_WeakData__efficiency")
        ):
            return self.__efficiency
        data = self.get_gathered_data(nps=nps, fs=fs, times=times)
        base_data = (
            data if 1 in nps else self.get_gathered_data(nps=1, fs=fs, times=times)
        )
        for case, base_case in zip(data.values(), base_data.values()):
            for (k, v), base_v in zip(case.items(), base_case.values()):
                base = base_v[0]
                case[k] = tuple(base / value if value else 0.0 for value in v)
        return data

    def __to_WeakGatheredData(self, nps: tuple[int, ...] = None, data=None):
        if nps is None:
            nps = self.nps
        if data is None:
            data = self.data
        one_case = next(iter(data.values()))[0]
        result = {
            f: {time: [float("inf")] * len(nps) for time in one_case[f].keys()}
            for f in one_case.keys()
            if isinstance(one_case[f], dict)
        }
        for i, np in enumerate(nps):
            for d in data[np]:
                for f, time_data in result.items():
                    for time, v_3 in d[f].items():
                        if v_3 < time_data[time][i]:
                            time_data[time][i] = v_3
        for v in result.values():
            for k, v_1 in v.items():
                v[k] = tuple(v_1)
        return result

    def draw_data(
        self,
        draw_type: str = "lines",
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (10, 10),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_data(draw_type: str, nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plots(ax, time, show_ylabel=True):
            for i, f in enumerate(fs):
                y = [float("inf")] * len(nps)
                for j, np in enumerate(nps):
                    y[j] = min(
                        y[j],
                        *(point[f][time] for point in self.data[np]),
                    )
                ax.plot(x, y, color=colors[i], marker=MARKER, label=f)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Wall-clock time (s)")

        def scatters(ax, time, show_ylabel=True):
            for i, f in enumerate(fs):
                y = []
                for j, np in enumerate(nps):
                    y.extend((point[f][time] for point in self.data[np]))
                ax.scatter(x, y, color=colors[i], marker=MARKER, label=f)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Wall-clock time (s)")

        colors = sns.color_palette("Set1", 9, 0.9)
        MARKER = "x"
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
            x = tuple(np for np in nps for _ in self.data[np])
            func = scatters
        if isinstance(times, tuple):
            fig, axs = plt.subplots(
                1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=True
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

    def draw_speedup(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (8, 8),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_speedup(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plot(ax, time, show_ylabel=True):
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
                y = speedup[f][time]
                ax.plot(nps, y, marker=MARKER, label=f"{f}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Scaled speed-up")

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
        speedup = self.get_speedup(nps=nps, fs=fs, times=times)
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
            fig, axs = plt.subplots(
                1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=True
            )
            if save_file_name is None or not separate_legend:
                for i, (ax, time) in enumerate(zip(axs, times)):
                    plot(ax, time, i == 0)
                    ax.legend()
            else:
                for i, (ax, time) in enumerate(zip(axs, times)):
                    plot(ax, time, i == 0)
            if save_file_name is None:
                plt.show()
            else:
                if separate_legend:
                    draw_legend(axs[0], figsize, save_file_name)
                plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
            plt.close(fig)
        else:
            fig, ax = plt.subplots(figsize=figsize)
            plot(ax, times)
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
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (8, 8),
        save_file_name: str = None,
        separate_legend: bool = True,
    ):
        """draw_efficiency(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str, separate_legend: bool)"""

        def plot(ax, time, show_ylabel=True):
            y = [1.0] * len(nps)
            ax.plot(
                nps,
                y,
                color=SPECIAL_LINE_COLOUR,
                linewidth=SPECIAL_LINE_WIDTH,
                label="theoretical",
                linestyle=SPECIAL_LINE_STYLE,
            )
            for f in fs:
                y = efficiency[f][time]
                ax.plot(nps, y, marker=MARKER, label=f)
            ax.set_xscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Efficiency")

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
        efficiency = self.get_efficiency(nps=nps, fs=fs, times=times)
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
            fig, axs = plt.subplots(
                1, 2, figsize=(figsize[0] * 2 * 1.1, figsize[1]), sharey=True
            )
            if save_file_name is None or not separate_legend:
                for i, (ax, time) in enumerate(zip(axs, times)):
                    plot(ax, time, i == 0)
                    ax.legend()
            else:
                for i, (ax, time) in enumerate(zip(axs, times)):
                    plot(ax, time, i == 0)
            if save_file_name is None:
                plt.show()
            else:
                if separate_legend:
                    draw_legend(axs[0], figsize, save_file_name)
                plt.savefig(f"{save_file_name}.pdf", bbox_inches="tight")
            plt.close(fig)
        else:
            fig, ax = plt.subplots(figsize=figsize)
            plot(ax, times)
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
