import os
import sys

fragment_file = sys.argv[1]
assembly_file = sys.argv[2]
stdout_file = sys.argv[3]

print(fragment_file)
print(assmebly_file)
print(stdout_file)

for threshold_weight in range(19, 9, -1):
    os.system("sed -i ""'4s/.*/thresholdweight=" + str(threshold_weight/100.0) + "/' config.txt")
    fname_prefix = str(threshold_weight).zfill(2)
    os.system("java ScatterSearch " + fname_prefix + fragment_file + " " + fname_prefix + assembly_file + "| tee " + fname_prefix + stdout_file)
    print("Done with " + str(threshold_weight) + " percent")
