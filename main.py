from copy import deepcopy
from itertools import combinations
from functools import reduce


def gluing(SNF):
    """
    This function is performs gluing operaration for SDNF
    """
    nf = []
    for i in range(len(SNF)):  # Run through the SNF in search of implicants
        for j in range(i + 1, len(SNF)):
            summand1 = set(SNF[i])
            summand2 = set(SNF[j])
            implicant = list(summand1 & summand2)
            implicant.sort(key=lambda x: x[-1])
            if len(implicant) == 2:  # If the implicant fits, then add it to NF
                nf.append(implicant)
    return nf


def calculation_method(function, key=None):
    """
    This function is performs calculation method for DNF
    """
    mnf = deepcopy(function)
    for index, _ in enumerate(mnf):  # Checking every implicant, if it is extra
        dict_of_variables = {}
        mnf_copy = deepcopy(mnf)
        verifiable_implicant = mnf_copy.pop(index)  # Remove checked implicant from the list
        for j in verifiable_implicant:
            dict_of_variables[j] = True if key == "dnf" else False
        for i in mnf_copy:
            for index, j in enumerate(i):  # substitution of values ​​of the checked implicant
                if j in dict_of_variables.keys():
                    i[index] = dict_of_variables[j]
                elif j[-1] in dict_of_variables.keys():
                    i[index] = not dict_of_variables[j[-1]]
                elif f"!{j}" in dict_of_variables.keys():
                    i[index] = not dict_of_variables[f"!{j}"]
        for i in range(len(mnf_copy)):
            for j in range(i + 1, len(mnf_copy)):
                result = set(mnf_copy[i]).symmetric_difference(set(mnf_copy[j]))
                if len(result) == 2 and 1 not in result:
                    mnf.remove(verifiable_implicant)
                    calculation_method(mnf, key)
    return mnf


def calculation_tabular_method(nf, SNF, key=None):
    """
    This function is performs calculation-tabular method for KNF
    """
    mnf, table, filled_columns, verifiable_implicants = [], [], [], []
    for i in nf:
        table.append([len(set(i) & set(j)) == 2 for j in SNF])
    filled_columns = [False for _ in range(len(table[0]))]
    for i in range(len(table[0])):
        verifiable_column = [j[i] for j in table]
        if verifiable_column.count(True) == 1:
            implicant = nf[verifiable_column.index(True)]
            if implicant not in mnf:
                mnf.append(implicant)
            filled_columns = list(map(
                lambda x: x[0] or x[1],
                zip(filled_columns, table[verifiable_column.index(True)])
            ))
    verifiable_implicants = [i for i in nf if i not in mnf]
    if False in filled_columns:
        min_amount = 256
        for amount in range(1, len(verifiable_implicants) + 1):
            for subset in combinations(verifiable_implicants, amount):
                set_of_verifiable_implicants = [table[nf.index(i)] for i in subset]
                set_of_verifiable_implicants = reduce(
                    lambda x, y: [i or j for i, j in zip(x, y)],
                    set_of_verifiable_implicants
                )
                if False not in list(map(
                        lambda x: x[0] or x[1],
                        zip(set_of_verifiable_implicants, filled_columns)
                )) and len(set_of_verifiable_implicants) < min_amount:
                    min_amount = len(set_of_verifiable_implicants)
        mnf.extend(subset)
    print("        ", end="")
    for i in SNF:
        sign = "*" if key == "sdnf" else "+"
        print(f"|   {sign.join(i).ljust(10, ' ')}", end="")
    print()
    for index, i in enumerate(table):
        print(f" {sign.join(nf[index]).rjust(6, ' ')} ", end="")
        for j in i:
            if j:
                print("|      x      ", end="")
            else:
                print("|             ", end="")
        print()
    return mnf


def tabular_method(SNF):
    """
    This function  performs calculation-tabular method for KNF
    """
    SNF_binary_list, karno_kart, mnf = [], [], []
    meanings = [
        ["0 0 0", "0 0 1", "0 1 1", "0 1 0"],
        ["1 0 0", "1 0 1", "1 1 1", "1 1 0"],
    ]
    for i in SNF:
        SNF_binary_list.append([str(int(not j.startswith("!"))) for j in i])
    for i in meanings:
        karno_kart.append([j.split() in SNF_binary_list for j in i])
    print("      |  00  |  01  |  11  |  10  ")
    for index, i in enumerate(karno_kart):
        print(f"  {index}   ", end="")
        for j in i:
            print(f"|  {int(j)}   ", end="")
        print()
    if karno_kart[0].count(True) % 2 == 1:
        for i in range(len(karno_kart[0])):
            verifiable_column = [j[i] for j in karno_kart]
            if all(verifiable_column):
                for index, j in enumerate(SNF_binary_list):
                    if j == meanings[0][i].split():
                        mnf.append(SNF[index][1:])
                        for k in karno_kart:
                            k[i] = False
    for index, karno_kart_string in enumerate(karno_kart):
        for iterator in range(len(karno_kart_string)):
            for iterator_for_comparison in range(iterator + 1, len(karno_kart_string)):
                if karno_kart_string[iterator] and karno_kart_string[iterator_for_comparison]:
                    check = meanings[index][iterator].split()[:-1]
                    for index_, SNF_bin_item in enumerate(SNF_binary_list):
                        if SNF_bin_item[:-1] == check:
                            mnf.append(SNF[index_][:-1])
                            karno_kart_string[iterator] = karno_kart_string[iterator_for_comparison] = False
                            break
    return mnf


def check(checked_function, key=None):
    checked_function = deepcopy(checked_function)
    test_1 = checked_function.count("(") == checked_function.count(")")
    test_3 = checked_function.count("a") == checked_function.count("b") == checked_function.count("c")
    if key == "dnf":
        checked_function = [i.split("*") for i in checked_function.split(" + ")]
    elif key == "knf":
        checked_function = [i.split("+") for i in checked_function[1:-1].split(") * (")]
    for i in checked_function:
        test_2 = len(i) == 3
        if not test_2:
            break
    test_4 = type(checked_function) == str
    if key == "dnf":
        test_5 = checked_function.count("a") - checked_function.count(" + ") == 1
    elif key == "knf":
        test_5 = checked_function.count("a") - checked_function.count(" * ") == 1

    return all([test_1, test_2, test_3])


def print_dnf(dnf):
    dnf_output = ""
    for i in dnf:
        dnf_output += f"{'*'.join(i)}+"
    print(f"DNF: {dnf_output[:-1]}")


def print_knf(knf):
    knf_output = ""
    for i in knf:
        knf_output += f"({'+'.join(i)})*"
    print(f"KNF: {knf_output[:-1]}")


def print_mdnf(mdnf):
    mdnf_output = ""
    for i in mdnf:
        mdnf_output += f"{'*'.join(i)}+"
    print(f"MDNF: {mdnf_output[:-1]}")


def print_mknf(mknf):
    mknf_output = ""
    for i in mknf:
        mknf_output += f"({'+'.join(i)})*"
    print(f"MKNF: {mknf_output[:-1]}")


def main():
    try:
        SDNF = ("!a*b*c + a*!b*!c + a*!b*c + a*b*c")
        print(SDNF)
        SKNF = ("(a+b+c) * (a+b+!c) * (a+!b+c) * (!a+!b+c)")
        print(SKNF)
        if check(SDNF, "dnf") and check(SKNF, "knf"):
            SDNF = [i.split("*") for i in SDNF.split(" + ")]
            SKNF = [i.split("+") for i in SKNF[1:-1].split(") * (")]

            print("Gluing")
            dnf = gluing(SDNF)
            knf = gluing(SKNF)
            print_dnf(dnf)
            print_knf(knf)

            print("Calculation method:")
            mdnf = calculation_method(dnf, "dnf")
            mknf = calculation_method(knf, "knf")
            print_mdnf(mdnf)
            print_mknf(mknf)

            print("Calculation-tabular method:")
            mdnf = calculation_tabular_method(dnf, SDNF, "sdnf")
            print_mdnf(mdnf)
            mknf = calculation_tabular_method(knf, SKNF, "sknf")
            print_mknf(mknf)

            print("Tabular method:")
            mdnf = tabular_method(SDNF)
            print_mdnf(mdnf)
            mknf = tabular_method(SKNF)
            print_mknf(mknf)
        else:
            print("mistake(")
    except Exception:
        print("!!!something went wrong!!!")


if __name__ == '__main__':
    main()