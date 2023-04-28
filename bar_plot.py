#!/usr/bin/env python
import matplotlib.pyplot as plt
import yaml
import sys
from pandas import DataFrame


def get_bar_chart(data: dict[int, int], bins: int, title: str, output_image: str) -> None:
    # Get the frequency of each slice of bin
    bin_size = max(data.keys())//bins + 10
    # aproximate the value of the bin to the nearest integer divisible by 10
    bin_size = int(bin_size/10)*10
    # Create a list of bins
    bin_list = [bin for bin in range(0, max(data.keys()), bin_size)]
    # Create a list of frequencies
    frequency_list = [0 for _ in range(len(bin_list))]
    # Fill the frequency list
    for key, value in data.items():
        index = key//bin_size
        frequency_list[index] += value
    # Create a dataframe
    sum_freq_list = sum(frequency_list)
    frequency_list = [value/sum_freq_list
                      for value in frequency_list]  # to be in percentage
    top_10 = sorted(frequency_list, reverse=True)[0:10]
    positions_dict = {}
    for i, value in enumerate(frequency_list):
        if value in top_10:
            # print(value, top_10.index(value))
            # print(bin_list[i-1], bin_list[i])
            positions_dict[top_10.index(value) + 1] = [bin_list[i-1], bin_list[i]]

    df = DataFrame({"bins": bin_list, "frequency": frequency_list})
    # Create the bar chart
    fig, ax = plt.subplots()
    ax.bar(df["bins"], df["frequency"], width=bin_size,
           color="#3686C9", edgecolor='black', align='edge')
    ax.set_title(title)
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("Frequency of SNPs")
    # plt.text((21555+ 266) // 2, -0.025, 'ORF1AB')
    # plt.text((21563+25384) // 2 , -0.025, 'S')
    # plt.text((25393+26220) // 2 , -0.045, 'ORF3a', rotation=90)
    # plt.text((26245+26472) // 2 , -0.045, 'E', rotation=90)
    # plt.text((26523+27191) // 2 , -0.045, 'M', rotation=90)
    # plt.text((27202+27387) // 2 , -0.045, 'ORF6', rotation=90)
    # plt.text((27394+27759) // 2 , -0.045, 'ORF7a', rotation=90)
    # plt.text((27894+28259) // 2 , -0.045, 'ORF8', rotation=90)
    # plt.text((28274+29533) // 2 , -0.045, 'N', rotation=90)
    # plt.text((29558+29674) // 2 , -0.045, 'ORF10', rotation=90)
    plt.savefig(
        f"{output_image}")
    return fig, positions_dict


def scatter(x, y, title, xlabel,ylabel):
    fig, ax = plt.subplots(1, 1)
    ax.scatter(x, y, marker="o", edgecolors="face",
               alpha=0.4, color="#3686C9", s=3)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Fazer curva em que mostra que Depth * Freq > 0.05 * Depth
    depth, freq = ([], [])
    for d, f in zip(x, y):
        if d * f < 0.05 * d:
            depth.append(d)
            freq.append(f)
    ax.scatter(depth, freq, marker="o", edgecolors="face",
               alpha=0.4, color="red", s=3)
    return fig


def get_positions_dict(yaml_dict: dict[str, str]) -> dict[int, int]:
    data_dict = {}

    for key, val in yaml_dict.items():
        data_dict[int(key)] = int(val)

    return data_dict


def main():
    input_file = sys.argv[1]
    output_image = sys.argv[2]
    with open(input_file, "r") as handle:
        yaml_data = yaml.load(handle, Loader=yaml.FullLoader)
    print(yaml_data)
    info_dict = get_positions_dict(yaml_data)
    get_bar_chart(data=info_dict, bins=100,
                  title="Monkeypox SNPs", output_image=output_image)


if __name__ == "__main__":
    main()
