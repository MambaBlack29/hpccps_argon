import matplotlib.pyplot as plt


# X-axis: number of points
points = [1600, 5400, 12800]

# Y-axis: serial wall time (seconds)
without_verlet = [7.98, 77.56, 427.49]
with_verlet_10 = [2.37, 10.75, 50.64]
with_verlet_15 = [1.79, 8.34, 34.88]
with_verlet_20 = [1.60, 6.92, 28.22]


plt.figure(figsize=(8, 5))
# plt.plot(points, without_verlet, marker="o", linewidth=2, label="Without Verlet")
plt.plot(points, with_verlet_10, marker="s", linewidth=2, label="With Verlet (update every 10 steps)")
plt.plot(points, with_verlet_15, marker="^", linewidth=2, label="With Verlet (update every 15 steps)")
plt.plot(points, with_verlet_20, marker="d", linewidth=2, label="With Verlet (update every 20 steps)")

plt.xlabel("Number of points")
plt.ylabel("Serial wall time (s)")
plt.title("Serial Wall Time vs Number of Points")
plt.grid(True, linestyle="--", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()
