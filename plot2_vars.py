import pandas as pd
import matplotlib.pyplot as plt


fname = "omp.csv"
df = pd.read_csv(fname, comment="#")
print(df)

var_names = list(df.columns)

print("var names =", var_names)

# split the df into individual vars
# assumption: column order - 0=problem size, 1=blas time, 2=basic time

problem_sizes = df[var_names[0]].values.tolist()
code1_time = df[var_names[1]].values.tolist()
code2_time = df[var_names[2]].values.tolist()
code3_time = df[var_names[3]].values.tolist()
code4_time = df[var_names[4]].values.tolist()
#code5_time = df[var_names[5]].values.tolist()

plt.title("OMP Speed up")

xlocs = [i for i in range(len(problem_sizes))]

plt.xticks(xlocs, problem_sizes)

# here, we are plotting the raw values read from the input .csv file, which
# we interpret as being "time" that maps directly to the y-axis.
#
# what if we want to plot MFLOPS instead? How do we compute MFLOPS from
# time and problem size? You may need to add some code here to compute
# MFLOPS, then modify the plt.plot() lines below to plot MFLOPS rather than time.

plt.plot(code1_time, "r-o")
plt.plot(code2_time, "b-x")
plt.plot(code3_time, "g-^")
plt.plot(code4_time,"y-.")
#plt.plot(code5_time,"k-^")

#plt.xscale("log")
#plt.yscale("log")

plt.xlabel("problem_size(n)")
plt.ylabel("Speed up")

varNames = [var_names[1], var_names[2], var_names[3], var_names[4]]
plt.legend(varNames, loc="best")

plt.grid(axis='both')

plt.show()

# EOF