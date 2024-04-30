from dict_smn_ii_stations import smn_ii
from dict_smn_iii_stations import smn_iii

def print_non_repeats(list1, list2):
    for str1 in list1:
        found = False
        for str2 in list2:
            if str1 in str2:
                found = True
                break
        if not found:
            print(str1)


list_smn_ii = []
for i in range(1,len(smn_ii)):
	list_smn_ii.append(smn_ii[i][0])

list_smn_iii = []
for ii in range(1,len(smn_iii)):
	list_smn_iii.append(smn_iii[ii][0])

# Check for repeats
print_non_repeats(list_smn_iii, list_smn_ii)
