import os
import sys

fragment_file = sys.argv[1]
assembly_file = sys.argv[2]
stdout_file = sys.argv[3]

print(fragment_file)
print(assembly_file)
print(stdout_file)

for threshold_weight in range(19, 9, -1):
    os.system("sed -i ""'4s/.*/thresholdweight=" + str(threshold_weight/100.0) + "/' config.txt")
    fname_suffix = str(threshold_weight).zfill(2)
    os.system("java ScatterSearch " + fragment_file + " " + assembly_file + fname_suffix + "| tee " + stdout_file + fname_suffix)
    print("Done with " + str(threshold_weight) + " percent")
