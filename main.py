import parser
import minimizer


def test_odd(function):
    pdnf = parser.make_pdnf(parser.build_truth_table(parser.resolve_inversions(function)))
    index = parser.to_number_form(parser.build_truth_table(pdnf))
    func = minimizer.find_odd(pdnf)
    index_of_minimized = parser.to_number_form(parser.build_truth_table(func))
    if index == index_of_minimized:
        print(function, " |||| ", func, " |||| ", "test_odd passed")
    else:
        print(function, " |||| ", func, " |||| ", "test_odd failed")


def test_Quine(function):
    pdnf = parser.make_pdnf(parser.build_truth_table(parser.resolve_inversions(function)))
    index = parser.to_number_form(parser.build_truth_table(pdnf))
    func = minimizer.minimize_Quine(pdnf)
    index_of_minimized = parser.to_number_form(parser.build_truth_table(func))
    if index == index_of_minimized:
        print(function, " |||| ", func, " |||| ", "test_Quine passed")
    else:
        print(function, " |||| ", func, " |||| ", "test_Quine failed")


def test_KMap(function):
    pdnf = parser.make_pdnf(parser.build_truth_table(parser.resolve_inversions(function)))
    index = parser.to_number_form(parser.build_truth_table(pdnf))
    func = minimizer.minimize_KMap(pdnf)
    index_of_minimized = parser.to_number_form(parser.build_truth_table(func))
    if index == index_of_minimized:
        print(function, " |||| ", func, " |||| ", "test_KMap passed")
    else:
        print(function, " |||| ", func, " |||| ", "test_KMap failed")


print("~((b+c)*~(a*c))")

function = "~((b+c)*~(a*c))"

print(parser.build_truth_table(parser.resolve_inversions(function)))

pdnf = parser.make_pdnf(parser.build_truth_table(parser.resolve_inversions(function)))

print("PDNF: ", pdnf)

print("Calculation\n", minimizer.find_odd(pdnf))

print("Calculation-table\n", minimizer.minimize_Quine(pdnf))

print("Kmaps\n", minimizer.minimize_KMap(pdnf))

print("Find Odd test")
test_odd("~((~a+~b)*~(~a*~c))")
test_odd("~((~a+~b)*~(~a*c))")
test_odd("~((~a+b)*~(~a*~c))")
test_odd("~((~a+b)*~(~a*c))")
test_odd("~((a+~b)*~(a*~c))")
test_odd("~((a+~b)*~(a*c))")
test_odd("~((a+b)*~(a*~c))")
test_odd("~((a+b)*~(a*c))")
test_odd("~((~a+~b)*~(~b*~c))")
test_odd("~((~a+~b)*~(~b*c))")
test_odd("~((~a+b)*~(b*~c))")
test_odd("~((~a+b)*~(b*c))")
test_odd("~((a+~b)*~(~b*~c))")
test_odd("~((a+~b)*~(~b*c))")
test_odd("~((a+b)*~(b*~c))")
test_odd("~((a+b)*~(b*c))")
test_odd("~((~b+~c)*~(~a*~c))")
test_odd("~((~b+c)*~(~b*c))")
test_odd("~((b+~c)*~(~b*~c))")
test_odd("~((b+c)*~(~a*c))")
test_odd("~((~b+~c)*~(a*~c))")
test_odd("~((~b+c)*~(a*~c))")
test_odd("~((b+~c)*~(a*c))")
test_odd("~((b+c)*~(a*c))")
test_odd("~((~a+~c)*~(~b*~c))")
test_odd("~((~a+c)*~(~b*c))")
test_odd("~((~a+~c)*~(b*~c))")
test_odd("~((~a+c)*~(b*c))")
test_odd("~((a+~c)*~(~b*~c))")
test_odd("~((a+c)*~(~b*c))")
test_odd("~a*b+c")
print()


print("Quine test")
test_Quine("~((~a+~b)*~(~a*~c))")
test_Quine("~((~a+~b)*~(~a*c))")
test_Quine("~((~a+b)*~(~a*~c))")
test_Quine("~((~a+b)*~(~a*c))")
test_Quine("~((a+~b)*~(a*~c))")
test_Quine("~((a+~b)*~(a*c))")
test_Quine("~((a+b)*~(a*~c))")
test_Quine("~((a+b)*~(a*c))")
test_Quine("~((~a+~b)*~(~b*~c))")
test_Quine("~((~a+~b)*~(~b*c))")
test_Quine("~((~a+b)*~(b*~c))")
test_Quine("~((~a+b)*~(b*c))")
test_Quine("~((a+~b)*~(~b*~c))")
test_Quine("~((a+~b)*~(~b*c))")
test_Quine("~((a+b)*~(b*~c))")
test_Quine("~((a+b)*~(b*c))")
test_Quine("~((~b+~c)*~(~a*~c))")
test_Quine("~((~b+c)*~(~b*c))")
test_Quine("~((b+~c)*~(~b*~c))")
test_Quine("~((b+c)*~(~a*c))")
test_Quine("~((~b+~c)*~(a*~c))")
test_Quine("~((~b+c)*~(a*~c))")
test_Quine("~((b+~c)*~(a*c))")
test_Quine("~((b+c)*~(a*c))")
test_Quine("~((~a+~c)*~(~b*~c))")
test_Quine("~((~a+c)*~(~b*c))")
test_Quine("~((~a+~c)*~(b*~c))")
test_Quine("~((~a+c)*~(b*c))")
test_Quine("~((a+~c)*~(~b*~c))")
test_Quine("~((a+c)*~(~b*c))")
test_Quine("~a*b+c")
print()


print("KMap test")
test_KMap("~((~a+~b)*~(~a*~c))")
test_KMap("~((~a+~b)*~(~a*c))")
test_KMap("~((~a+b)*~(~a*~c))")
test_KMap("~((~a+b)*~(~a*c))")
test_KMap("~((a+~b)*~(a*~c))")
test_KMap("~((a+~b)*~(a*c))")
test_KMap("~((a+b)*~(a*~c))")
test_KMap("~((a+b)*~(a*c))")
test_KMap("~((~a+~b)*~(~b*~c))")
test_KMap("~((~a+~b)*~(~b*c))")
test_KMap("~((~a+b)*~(b*~c))")
test_KMap("~((~a+b)*~(b*c))")
test_KMap("~((a+~b)*~(~b*~c))")
test_KMap("~((a+~b)*~(~b*c))")
test_KMap("~((a+b)*~(b*~c))")
test_KMap("~((a+b)*~(b*c))")
test_KMap("~((~b+~c)*~(~a*~c))")
test_KMap("~((~b+c)*~(~b*c))")
test_KMap("~((b+~c)*~(~b*~c))")
test_KMap("~((b+c)*~(~a*c))")
test_KMap("~((~b+~c)*~(a*~c))")
test_KMap("~((~b+c)*~(a*~c))")
test_KMap("~((b+~c)*~(a*c))")
test_KMap("~((b+c)*~(a*c))")
test_KMap("~((~a+~c)*~(~b*~c))")
test_KMap("~((~a+c)*~(~b*c))")
test_KMap("~((~a+~c)*~(b*~c))")
test_KMap("~((~a+c)*~(b*c))")
test_KMap("~((a+~c)*~(~b*~c))")
test_KMap("~((a+c)*~(~b*c))")
test_KMap("~a*b+c")