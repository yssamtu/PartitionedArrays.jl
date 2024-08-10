from collections import defaultdict
from json import load
from numpy import prod
from os import listdir
from os.path import isdir, join
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union


class StrongData:
    def __init__(self, data):
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
            for data_size, marker in zip(data_sizes, markers):
                for i, f in enumerate(fs):
                    y = [float("inf")] * len(nps)
                    for j, np in enumerate(nps):
                        y[j] = min(
                            y[j],
                            *(point[f][time] for point in self.data[data_size][np]),
                        )
                    ax.plot(
                        x, y, color=colors[i], marker=marker, label=f"{f}_{data_size}"
                    )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Wall-clock time (s)")

        def scatters(ax, time, show_ylabel=True):
            for data_size, marker in zip(data_sizes, markers):
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
        markers = ("x", "+", "v", ".", "*")
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
                1, 2, figsize=(figsize[0] * 2, figsize[1]), sharey=True
            )
            ylabel_mask = (True, False)
            if draw_limit and len(data_sizes) == 1:
                if save_file_name is None or not separate_legend:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        func(ax, time, show_ylabel)
                        ax.axvline(
                            x=pow(data_sizes[0] - 1, 3) / 25000,
                            color="#AFA938",
                            linewidth=2,
                            label="25K rows/rank",
                        )
                        ax.legend()
                else:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        func(ax, time, show_ylabel)
                        ax.axvline(
                            x=pow(data_sizes[0] - 1, 3) / 25000,
                            color="#AFA938",
                            linewidth=2,
                            label="25K rows/rank",
                        )
            else:
                if save_file_name is None or not separate_legend:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        func(ax, time, show_ylabel)
                        ax.legend()
                else:
                    for ax, time, show_ylabel in zip(axs, times, ylabel_mask):
                        func(ax, time, show_ylabel)
        else:
            fig, ax = plt.subplots(figsize=figsize)
            func(ax, times)
            if draw_limit and len(data_sizes) == 1:
                ax.axvline(
                    x=pow(data_sizes[0] - 1, 3) / 25000,
                    color="#AFA938",
                    linewidth=2,
                    label="25K rows/rank",
                )
            if save_file_name is None or not separate_legend:
                ax.legend()
        if save_file_name is None:
            plt.show()
        else:
            if separate_legend:
                fig_legend, ax_legend = plt.subplots(
                    figsize=(figsize[0] * 0.5, figsize[1] * 0.5)
                )
                handles, labels = ax.get_legend_handles_labels()
                legend = fig_legend.legend(
                    handles, labels, loc="center", ncol=2, frameon=False
                )
                ax_legend.axis("off")
                fig_legend.canvas.draw()
                bbox = legend.get_window_extent().transformed(
                    fig_legend.dpi_scale_trans.inverted()
                )
                fig_legend.set_size_inches(bbox.width, bbox.height)
                fig_legend.savefig(f"{save_file_name}_legend.pdf")
                plt.close(fig_legend)
            plt.savefig(f"{save_file_name}.pdf")
        plt.close(fig)

    def draw_complexity(
        self,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (6, 6),
        save_file_name: str = None,
    ):
        """draw_complexity(nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str)"""
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
        marker = "x"
        x = tuple(k**3 for k in self.data.keys())

        def plot(ax, f, time, show_y_label=True):
            for np in nps:
                y = []
                for data_size in self.data.keys():
                    y.append(min(case[f][time] for case in self.data[data_size][np]))
                ax.plot(
                    x,
                    tuple(i / y[-1] * y[-1] for i in y),
                    marker=marker,
                    label=f"{f}_{np}",
                )
            ax.plot(
                x, tuple((i / x[-1]) * y[-1] for i in x), marker=marker, label=r"$O(n)$"
            )
            ax.plot(
                x,
                tuple((i / x[-1]) ** 0.5 * y[-1] for i in x),
                marker=marker,
                label=r"$O(\sqrt n)$",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("Problem size")
            if show_y_label:
                ax.set_ylabel("Wall-clock time (s)")
            ax.legend()

        if isinstance(times, tuple):
            figsize = (figsize[0] * 2, figsize[1])
            for f in fs:
                fig, axs = plt.subplots(1, 2, figsize=figsize)
                plot(axs[0], f, times[0], True)
                plot(axs[1], f, times[1], False)
                if save_file_name is None:
                    plt.show()
                else:
                    plt.savefig(f"{save_file_name}_{f}.pdf")
                plt.close(fig)
        else:
            for f in fs:
                fig, ax = plt.subplots(figsize=figsize)
                plot(ax, f, times)
                if save_file_name is None:
                    plt.show()
                else:
                    plt.savefig(f"{save_file_name}_{f}.pdf")
                plt.close(fig)

    def draw_speedup(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (8, 8),
        save_file_name: str = None,
    ):
        """draw_speedup(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str)"""

        def plot(data_size, ax, time, show_ylabel=True):
            y = nps
            ax.plot(nps, y, label="theoretical", linestyle="-.")
            for f in fs:
                y = speedup[data_size][f][time]
                ax.plot(nps, y, marker=marker, label=f"{f}_{data_size}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"CPU cores")
            if show_ylabel:
                ax.set_ylabel(r"Speed-up")
            ax.legend()

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
        marker = "x"
        if isinstance(times, tuple):
            for data_size in data_sizes:
                fig, axs = plt.subplots(
                    1, 2, figsize=(figsize[0] * 2, figsize[1]), sharey=True
                )
                plot(data_size, axs[0], times[0], True)
                plot(data_size, axs[1], times[1], False)
                if save_file_name is None:
                    plt.show()
                else:
                    plt.savefig(f"{save_file_name}.pdf")
                plt.close(fig)
        else:
            for data_size in data_sizes:
                fig, ax = plt.subplots(figsize=figsize)
                plot(data_size, ax, times)
                if save_file_name is None:
                    plt.show()
                else:
                    plt.savefig(f"{save_file_name}.pdf")
                plt.close(fig)

    def draw_efficiency(
        self,
        data_sizes: Union[int, tuple[int, ...], list[int]] = None,
        nps: Union[int, tuple[int, ...], list[int]] = None,
        fs: Union[str, tuple[str, ...], list[str]] = None,
        times: Union[str, tuple[str, str], list[str]] = None,
        figsize: tuple[int, int] = (8, 8),
        save_file_name: str = None,
    ):
        """draw_efficiency(data_sizes: Union[int, tuple[int, ...], list[int]], nps: Union[int, tuple[int, ...], list[int]], fs: Union[str, tuple[str, ...], list[str]], times: Union[str, tuple[str, str], list[str]], figsize: tuple[int, int], save_file_name: str)"""

        def plot(data_size, ax, time, show_ylabel=True):
            y = [1.0] * len(nps)
            ax.plot(nps, y, label="theoretical", linestyle="-.")
            for f in fs:
                y = efficiency[data_size][f][time]
                ax.plot(nps, y, marker=marker, label=f"{f}_{data_size}")
            ax.set_xscale("log")
            ax.set_xlabel("CPU cores")
            if show_ylabel:
                ax.set_ylabel("Efficiency")
            ax.legend()

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
        marker = "x"
        if isinstance(times, tuple):
            for data_size in data_sizes:
                fig, axs = plt.subplots(
                    1, 2, figsize=(figsize[0] * 2, figsize[1]), sharey=True
                )
                plot(data_size, axs[0], times[0], True)
                plot(data_size, axs[1], times[1], False)
                if save_file_name is None:
                    plt.show()
                else:
                    plt.savefig(f"{save_file_name}.pdf")
                plt.close(fig)
        else:
            for data_size in data_sizes:
                fig, ax = plt.subplots(figsize=figsize)
                plot(data_size, ax, times)
                if save_file_name is None:
                    plt.show()
                else:
                    plt.savefig(f"{save_file_name}.pdf")
                plt.close(fig)


def load_strong_scaling_data(folder: str, computer: str = "das5"):
    """load_strong_scaling_data(folder: str)"""

    def convert_to_dict(d):
        def f(d):
            if isinstance(d, defaultdict):
                d = {k: convert_to_dict(v) for k, v in sorted(d.items())}
            return d

        if isinstance(d, defaultdict):
            d = {k: f(v) for k, v in reversed(sorted(d.items()))}
        return d

    if computer == "das5":
        folder = join(folder, "strong_scaling")
    else:
        folder = join(folder, "strong_scaling_old_snellius")
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
                    {**load(f), "node_core": node_core, "part": part}
                )
    data = convert_to_dict(dicts)
    for v in data.values():
        for k, v_1 in v.items():
            v[k] = tuple(v_1)
    return StrongData(data)
