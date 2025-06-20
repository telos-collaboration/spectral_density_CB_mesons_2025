#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import standard_plot_main, channel_color, channel_marker


def plot(data, **kwargs):
    fig, ax = plt.subplots(1, 1, figsize=(7, 2.4), layout="constrained")

    # ax.set_xlim(0, 0.16)
    ax.set_xlabel(r"ensemble")
    ax.set_ylabel(r"$aE_n$")

    channels = ["ps", "v", "t", "av", "at", "s"]
    rep = "f"

    for channel in channels:
        to_plot = []
        for datum in data:
            if f"gevp_{rep}_{channel}_E0_mass_samples" not in datum:
                continue

            for n in range(3):
                name = datum["ensemble_name"]
                X = (
                    int(name[1]) + channels.index(channel) / 6 - 0.5
                )  # assigning position for each ensamble (not a good way)
                Y = datum[f"gevp_{rep}_{channel}_E{n}_mass_samples"].samples

                to_plot.append((Y.mean(), Y.std(), X))

        y_values, y_errors, x_values = zip(*to_plot)

        # print(np.array(to_fit_x).shape)

        ax.errorbar(
            x_values,
            y_values,
            yerr=y_errors,
            ls="none",
            alpha=0.7,
            color=channel_color(channel),
            label=channel,
            marker=channel_marker(channel),
        )

    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(["M1", "M2", "M3", "M4", "M5"])

    handles, labels = fig.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(
        by_label.values(),
        by_label.keys(),
        loc="outside upper center",
        ncol=6,
        borderaxespad=0.2,
    )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
