
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


# Controlled Hadamard using basic gates
def controlled_h(circ, control, target):
    circ.ry(-pi / 4, target)
    circ.h(target)
    circ.cx(control, target)
    circ.h(target)
    circ.ry(pi / 4, target)


# Swap Qiskit's bit order
def rename(counts):
    new_counts = {"00": 0, "01": 0, "10": 0, "11": 0}

    for key, value in counts.items():
        new_counts[key[::-1]] = value

    return new_counts


def run_circuit(alpha, beta, theta0, theta1, chi, runs):

    # q0, q1 = original qubits
    # q2 = ancilla recording first measurement
    circ = QuantumCircuit(3, 2)

    # Initialize alpha|0> + beta|1> on q0
    if not np.isclose(abs(alpha)**2 + abs(beta)**2, 1):
        raise ValueError("Initial state must be normalized")

    theta = 2 * np.arctan2(abs(beta), abs(alpha))
    phase = np.angle(beta) - np.angle(alpha)

    circ.ry(theta, 0)
    circ.p(phase, 0)

    # First part of original circuit
    circ.ry(2 * theta0, 1)
    circ.cry(2 * theta1 - 2 * theta0, 0, 1)

    # Store first measurement result in q2
    circ.cx(1, 2)

    # Operations conditioned on first outcome = 1

    # Conditional X on q1
    circ.cx(2, 1)

    # Conditional H on q0
    controlled_h(circ, 2, 0)

    # Conditional Ry on q1
    circ.cry(2 * chi, 2, 1)

    # Doubly controlled Ry on q1
    ccry(circ, pi - 4 * chi, 0, 2, 1)

    # Conditional H on q0
    controlled_h(circ, 2, 0)

    # Final measurements
    circ.measure(2, 0)
    circ.measure(1, 1)

    print(circ.draw("text"))

    # Qollab built-in simulator
    job = backend.run(circ, shots=runs)

    out = job.result().get_counts()
    out = rename(out)

    for key in out:
        out[key] /= runs

    return out


def main():

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

    ax.set_xlabel("First measurement, second measurement")
    ax.set_ylabel("Probability")
    ax.set_title("Measurement statistics")
    ax.set_ylim(0, max(statistics.values()) * 1.15)

    plt.show()
