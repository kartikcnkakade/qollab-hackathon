
from qiskit import QuantumCircuit
from math import pi
import numpy as np
import matplotlib.pyplot as plt


# Doubly controlled Ry rotation
def ccry(circ, angle, control1, control2, target):
    circ.cry(angle / 2, control2, target)
    circ.cx(control1, control2)
    circ.cry(-angle / 2, control2, target)
    circ.cx(control1, control2)
    circ.cry(angle / 2, control1, target)


# R01(theta) acting on qubits 1 and 2
def R01(circ, theta):
    circ.x(2)
    circ.cry(2 * theta, 2, 1)
    circ.x(2)


# Controlled R01(theta), control qubit 0
def controlled_R01(circ, theta):
    circ.x(2)
    ccry(circ, 2 * theta, 0, 2, 1)
    circ.x(2)


# R12(theta) acting on qubits 1 and 2
def R12(circ, theta):
    circ.cx(2, 1)
    circ.cry(2 * theta, 1, 2)
    circ.cx(2, 1)


# Controlled R12(theta), control qubit 0
def controlled_R12(circ, theta):
    circ.cx(2, 1)
    ccry(circ, 2 * theta, 0, 1, 2)
    circ.cx(2, 1)


# Convert Qiskit bitstrings to labels 0, 1, 2, 3
def rename(counts):
    new_counts = {"0": 0, "1": 0, "2": 0, "3": 0}
    trans = {"00": "0", "10": "2",
             "01": "1", "11": "3"}

    for key, value in counts.items():
        new_counts[trans[key]] = value

    return new_counts


def run_circuit(alpha, beta, theta0, theta1, chi, runs):

    circ = QuantumCircuit(3, 2)

    # Initialize qubit 0 in alpha|0> + beta|1>
    if not np.isclose(abs(alpha)**2 + abs(beta)**2, 1):
        raise ValueError("Initial state must be normalized")

    theta = 2 * np.arctan2(abs(beta), abs(alpha))
    phase = np.angle(beta) - np.angle(alpha)

    circ.ry(theta, 0)
    circ.p(phase, 0)

    # First pair of rotations
    R01(circ, theta0)
    controlled_R01(circ, theta1 - theta0)

    # Hadamard on qubit 0
    circ.h(0)

    # Second pair of rotations
    R12(circ, chi)
    controlled_R12(circ, pi / 2 - 2 * chi)

    # Hadamard on qubit 0
    circ.h(0)

    # Final measurements
    circ.measure(1, 0)
    circ.measure(2, 1)

    print(circ.draw("text"))

    # Qollab's selected built-in simulator
    job = backend.run(circ, shots=runs)

    out = job.result().get_counts()
    out = rename(out)

    for key in out:
        out[key] /= runs

    return out


statistics = run_circuit(
    1 / np.sqrt(2),
    -1 / np.sqrt(2),
    pi / 4,
    pi / 2,
    pi / 3,
    10000
)

print("Measurement probabilities:")
print(statistics)

# Histogram with numerical labels
fig, ax = plt.subplots()

bars = ax.bar(
    list(statistics.keys()),
    list(statistics.values())
)

ax.bar_label(bars, fmt="%.4f", padding=3)

ax.set_xlabel("Measurement outcome")
ax.set_ylabel("Probability")
ax.set_title("Measurement statistics")
ax.set_ylim(0, max(statistics.values()) * 1.15)

plt.show()
