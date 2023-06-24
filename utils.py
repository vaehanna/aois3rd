import re
import numpy
from enum import Enum
from parser import build_truth_table


class TYPE_OF_FUNC(Enum):
    CONJUNCTIVE = 0
    DISJUNCTIVE = 1

class ImplicantIsNecessary(Exception):
    pass

class ImplicantIsUnnecessary(Exception):
    pass

def translate_to_implicant(pair, entrance):
    surrounding_translator = [
        [('~a', '~b'), ('~a', '~c'), ('~b', '~c'), ('~b'), ('~c'), ('~a', '~b', '~c')],
        [('~a', 'c'), ('~a', '~b'), ('~b', 'c'), ('c'), ('~b'), ('~a', '~b', 'c')],
        [('~a', 'b'), ('~a', 'c'), ('b', 'c'), ('b'), ('c'), ('~a', 'b', 'c')],
        [('~a', '~c'), ('~a', 'b'), ('b', '~c'), ('~c'), ('b'), ('~a', 'b', '~c')],
        [('a', '~b'), ('a', '~c'), ('~b', '~c'), ('~b'), ('~c'), ('a', '~b', '~c')],
        [('a', 'c'), ('a', '~b'), ('~b', 'c'), ('c'), ('~b'), ('a', '~b', 'c')],
        [('a', 'b'), ('a', 'c'), ('b', 'c'), ('b'), ('c'), ('a', 'b', 'c')],
        [('a', '~c'), ('a', 'b'), ('b', '~c'), ('~c'), ('b'), ('a', 'b', '~c')],
        ]
    index = 0
    for row in range(2):
        for col in range(4):
            if pair == (row, col):
                return surrounding_translator[index][entrance]
            index += 1
    raise Exception("Index is out of range!")


def find_type_of_function(function):
    result_flag = re.search(r"^~?(\w|True|False)(\*~?(\w|True|False))*(\s\+\s~?(\w|True|False)(\*~?(\w|True|False))*)*$", function)
    type_of_func = TYPE_OF_FUNC.DISJUNCTIVE
    if result_flag is None:
        result_flag = re.search(r"^\(?~?\w(\+~?\w\)?)*(\s\*\s\(?~?\w(\+~?\w\)?)*)*$", function)
        type_of_func = TYPE_OF_FUNC.CONJUNCTIVE
    if result_flag is None:
        raise Exception("Invalid input!")
    return type_of_func

#gets array of implicants as an input
def convert_to_eval(implicants, function_type):
        ##WARN! Could break if you'll change the Enum "TYPE_OF_FUNC"
        #concatenate expression to one string
        if (function_type == TYPE_OF_FUNC.CONJUNCTIVE):
            implicants = [" or ".join(impl) for impl in implicants]
            eval_string = "(" + ") and (".join(implicants) + ")"
        else:
            #convert individual implicants to strings
            implicants = [" and ".join(impl) for impl in implicants]
            eval_string = "(" + ") or (".join(implicants) + ")"
        eval_string = re.sub('~', ' not ', eval_string)
        return eval_string

def convert_to_human(implicants, function_type):
        ##WARN! Could break if you'll change the Enum "TYPE_OF_FUNC"
        #concatenate expression to one string
        if (function_type == TYPE_OF_FUNC.CONJUNCTIVE):
            implicants = ["+".join(impl) for impl in implicants]
            eval_string = "(" +") * (".join(implicants) + ")"
        else:
            #convert individual implicants to strings
            implicants = ["*".join(impl) for impl in implicants]
            eval_string = "(" + ") + (".join(implicants)  + ")"
        eval_string = re.sub(' not ', '~', eval_string)
        return eval_string


def represent_in_values(function, implicant, first_value, second_value):
    char_1 = re.search(r'\w', implicant[0])
    char_1 = char_1.group()
    function = re.sub(char_1, str(first_value), function)
    char_2 = re.search(r'\w', implicant[1])
    char_2 = char_2.group()
    function = re.sub(char_2, str(second_value), function)
    function = re.sub(r'~0', 'True', function)
    function = re.sub(r'~1', 'False', function)
    function = re.sub(r'1', 'True', function)
    function = re.sub(r'0', 'False', function)
    return function


def split_function(function):
    type_of_func = find_type_of_function(function)
    function_splitted = []

    if type_of_func == TYPE_OF_FUNC.DISJUNCTIVE:
        function = function.split(' + ')
        for elem in function:
            function_splitted.append(elem.split('*'))
    else:
        function = function.replace('(', '')
        function = function.replace(')', '')
        function = function.split(' * ')
        for elem in function:
            function_splitted.append(elem.split('+'))
    return function_splitted


def joining_rule(function):
    type_of_func = find_type_of_function(function)
    function_splitted = split_function(function)

    result = set()
    first_elem_index = 0
    while first_elem_index < len(function_splitted):
        second_elem_index = 0
        while second_elem_index < len(function_splitted):
            first_el = function_splitted[first_elem_index]
            second_el = function_splitted[second_elem_index]
            if len(set(first_el) ^ set(second_el)) == 2:
                sum = set(first_el) & set(second_el)
                function_splitted.remove(first_el)
                function_splitted.remove(second_el)
                if type_of_func == TYPE_OF_FUNC.DISJUNCTIVE:
                    result.add("*".join(sum))
                else:
                    result.add("+".join(sum))
                first_elem_index = 0
                second_elem_index = 0
                continue
            second_elem_index += 1
        first_elem_index += 1
    if type_of_func == TYPE_OF_FUNC.DISJUNCTIVE:
        result = " + ".join(result)
    else:
        result = " * ".join(result)
    for remaining in function_splitted:
        if type_of_func == TYPE_OF_FUNC.DISJUNCTIVE:
            result += ' + ' + "*".join(remaining)
        else:
            result += ' * ' + "+".join(remaining)
    return result


def find_kernel(perfect_form):
    type_of_func = find_type_of_function(perfect_form)
    simple_form = joining_rule(perfect_form)
    perfect_form_splitted = split_function(perfect_form)
    simple_form_splitted = split_function(simple_form)

    table = numpy.zeros(shape=(len(simple_form_splitted),
                               len(perfect_form_splitted)))

    for constituent in range(len(perfect_form_splitted)):
        for implicant in range(len(simple_form_splitted)):
            if set(simple_form_splitted[implicant]).issubset(perfect_form_splitted[constituent]):
                table[implicant][constituent] = 1
    kernel = []
    kernel_search_result = numpy.count_nonzero(table == 1, axis=0)
    for index in range(len(kernel_search_result)):
        if kernel_search_result[index] == 1:
            for element in range(len(table[:, index])):
                if table[element][index] == 1:
                    if simple_form_splitted[element] not in kernel:
                        kernel.append(simple_form_splitted[element])
    kernel_result = ""
    if type_of_func == TYPE_OF_FUNC.DISJUNCTIVE:
        for implicant in range(len(kernel)):
            kernel[implicant] = "*".join(kernel[implicant])
            kernel_result += kernel[implicant] + ' + '
    else:
        for implicant in range(len(kernel)):
            kernel[implicant] = "+".join(kernel[implicant])
            kernel_result += '(' + kernel[implicant] + ')' + ' * '
    kernel_result = kernel_result[: len(kernel_result) - 3]
    return kernel_result


def build_KMap(perfect_form):
    truth_table = build_truth_table(perfect_form)
    KMap_template = numpy.zeros(shape=(2, 4))
    KMap_template[0][0] = truth_table[3][0]
    KMap_template[0][1] = truth_table[3][1]
    KMap_template[0][3] = truth_table[3][2]
    KMap_template[0][2] = truth_table[3][3]
    KMap_template[1][0] = truth_table[3][4]
    KMap_template[1][1] = truth_table[3][5]
    KMap_template[1][3] = truth_table[3][6]
    KMap_template[1][2] = truth_table[3][7]
    return KMap_template


def find_surrounding(KMap, position):
    surrounding = numpy.zeros(shape=(2, 3))
    surrounding[position[0]][1] = KMap[position]
    if position[0] == 0:
        if position[1] == 0:
            surrounding[0][0] = KMap[position[0]][KMap.shape[1] - 1]
            surrounding[0][2] = KMap[position[0]][position[1] + 1]
            surrounding[1][0] = KMap[position[0] + 1][KMap.shape[1] - 1]
            surrounding[1][1] = KMap[position[0] + 1][position[1]]
            surrounding[1][2] = KMap[position[0] + 1][position[1] + 1]
        elif position[1] == KMap.shape[1] - 1:
            surrounding[0][0] = KMap[position[0]][position[1] - 1]
            surrounding[0][2] = KMap[position[0]][0]
            surrounding[1][0] = KMap[position[0] + 1][position[1] - 1]
            surrounding[1][1] = KMap[position[0] + 1][position[1]]
            surrounding[1][2] = KMap[position[0] + 1][0]
        else:
            surrounding[0][0] = KMap[position[0]][position[1] - 1]
            surrounding[0][2] = KMap[position[0]][position[1] + 1]
            surrounding[1][0] = KMap[position[0] + 1][position[1] - 1]
            surrounding[1][1] = KMap[position[0] + 1][position[1]]
            surrounding[1][2] = KMap[position[0] + 1][position[1] + 1]
    elif position[0] == 1:
        if position[1] == 0:
            surrounding[0][0] = KMap[position[0] - 1][KMap.shape[1] - 1]
            surrounding[0][1] = KMap[position[0] - 1][position[1]]
            surrounding[0][2] = KMap[position[0] - 1][position[1] + 1]
            surrounding[1][0] = KMap[position[0]][KMap.shape[1] - 1]
            surrounding[1][2] = KMap[position[0]][position[1] + 1]
        elif position[1] == KMap.shape[1] - 1:
            surrounding[0][0] = KMap[position[0] - 1][position[1] - 1]
            surrounding[0][1] = KMap[position[0] - 1][position[1]]
            surrounding[0][2] = KMap[position[0] - 1][0]
            surrounding[1][0] = KMap[position[0]][position[1] - 1]
            surrounding[1][2] = KMap[position[0]][0]
        else:
            surrounding[0][0] = KMap[position[0] - 1][position[1] - 1]
            surrounding[0][1] = KMap[position[0] - 1][position[1]]
            surrounding[0][2] = KMap[position[0] - 1][position[1] + 1]
            surrounding[1][0] = KMap[position[0]][position[1] - 1]
            surrounding[1][2] = KMap[position[0]][position[1] + 1]
    return surrounding