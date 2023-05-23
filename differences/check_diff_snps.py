import csv
import sys
from scipy.stats import chi2_contingency


def get_snps(file_name):
    with open(file_name) as handler:
        reader = list(csv.reader(handler))
    my_set = set()

    for e in reader:
        my_set.add((e[0], e[2], e[4], e[5]))
    return my_set


def main():
    # snippy = sys.argv[1]
    # ivar = sys.argv[2]
    # s_set = get_snps(snippy)
    # i_set = get_snps(ivar)
    #
    # print("Snippy Entries: ", len(s_set))
    # print("iVar Entries: ", len(i_set))
    #
    # print("Intersection: ", len(s_set.intersection(i_set)))
    # print("Symmetric Difference: ", len(s_set.symmetric_difference(i_set)))
    # print("0,0", len(s_set.intersection(i_set)))
    #
    # print("1,0", len(s_set.difference(i_set)))
    # print("0,1", len(i_set.difference(s_set)))
    # # diffs = s_set.difference(i_set)
    # # diffs_2 = i_set.difference(s_set)
    # # diff = s_set.symmetric_difference(i_set)
    # # print(diffs)
    # # print(len(diffs), len(diffs_2), len(diff))
    #
    #
    # Define the contingency table
    table = [[20867, 269], [35, 0]]

# Perform the chi-squared test
    chi2, p, dof, expected = chi2_contingency(table)

# Print the results
    print(f"Chi-squared: {chi2:.4f}")
    print(f"p-value: {p:.4f}")


if __name__ == "__main__":
    main()
