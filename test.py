import sys
import copy
import timeit
from array import array
from collections import defaultdict
import os
import psutil

genome_seq = array('b', 'actgaaccgg')
print(genome_seq)

#"../data/mtb-single-mapping-report-all.sam"

# with open("../data/toy-report-all.sam") as sam_file:
#     count = 0
#     curr_line = next(sam_file)
#     while curr_line[0] == '@':
#         print(curr_line)
#         curr_line = next(sam_file)
#
#     multi = False
#     multi_count = 0
#     unique_count = 0
#     multireads = defaultdict(list)
#
#     prev_fields = curr_line.split('\t')
#     for curr_line in sam_file:
#         curr_fields = curr_line.split('\t')
#         if curr_fields[0] == prev_fields[0]:
#             multireads[prev_fields[0]].append(prev_fields)
#             multi = True
#             multi_count += 1
#         elif multi == True:
#             multireads[prev_fields[0]].append(prev_fields)
#             multi = False
#             multi_count += 1
#         else:
#             unique_count += 1
#         prev_fields = curr_fields
#
# print(unique_count, multi_count, len(multireads))
# print(sys.getsizeof(multireads) / 1024 ** 2)
# print(sys.getsizeof(prev_fields) * multi_count / 2**20)
#
# process = psutil.Process(os.getpid())
# mem = process.memory_info()[0] / float(2 ** 20)
# print(mem)


# zero_array = array('i', [0] * 1000000)
# bases = [copy.deepcopy(zero_array), copy.deepcopy(zero_array), copy.deepcopy(zero_array), copy.deepcopy(zero_array)]
#
# print(sys.getsizeof(zero_array) / 1000 ** 2)
# print(sys.getsizeof(bases))

# bases[0][0] = 10
# print(bases[0][0:5], bases[1][0:5])

# print(timeit.timeit("import copy; from array import array; a = array('i', [0] * 4000000); b = copy.deepcopy(a); c = copy.deepcopy(a); d = copy.deepcopy(a);", number=1))


# stat = """from array import array; int_array_iter = array('i')
# for i in range(4000000):
#     int_array_iter.append(0)"""
#
# print(timeit.timeit("from array import array; intarray = array('i', [0] * 4000000)", number=1))
# print(timeit.timeit("lst = [0] * 4000000", number=1))
# print(timeit.timeit(stat, number=1))
#
# intarray = array('i', [0] * 4000000)
# lst = [0] * 4000000
# int_array_iter = array('i')
# for i in range(4000000):
#     int_array_iter.append(0)
#
#
# print(sys.getsizeof(intarray) / 1024 ** 2)
# print(sys.getsizeof(lst) / 1024 ** 2)
# print(sys.getsizeof(int_array_iter) / 1024 ** 2)


