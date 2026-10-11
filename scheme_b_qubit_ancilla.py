from qiskit import QuantumCircuit, transpile
from math import pi
import numpy as np
import matplotlib.pyplot as plt


def rename(counts):
    new_counts = {"00": 0, "01": 0, "10": 0, "11": 0}
    trans = {"00": "00", "10": "01",
             "01": "10", "11": "11"}

    for key, value in counts.items():
        new_counts[trans[key]] = value

    return new_counts


def run_circuit(alpha, beta, theta0, theta1, chi, runs):

    circ = QuantumCircuit(2, 2)

    # Initialize alpha|0> + beta|1>
    if not np.isclose(abs(alpha)**2 + abs(beta)**2, 1):
        raise ValueError("Initial state must be normalized")

    theta = 2 * np.arctan2(abs(beta), abs(alpha))
    phase = np.angle(beta) - np.angle(alpha)

    circ.ry(theta, 0)
    circ.p(phase, 0)

    # First part of the circuit
    circ.ry(2 * theta0, 1)
    circ.cry(2 * theta1 - 2 * theta0, 0, 1)

    # Mid-circuit measurement
    circ.measure(1, 0)

    # Classical feed-forward
    with circ.if_test((circ.clbits[0], 1)):
        circ.x(1)
        circ.h(0)
        circ.ry(2 * chi, 1)
        circ.cry(pi - 4 * chi, 0, 1)
        circ.h(0)

    # Final measurement
    circ.measure(1, 1)

    # Display graphical circuit
    circuit_fig = circ.draw(output="mpl", fold=-1)
    circuit_fig.set_size_inches(15, 5)
    plt.show()

    # Compile for the selected IBM backend
    compiled_circ = transpile(circ, backend)

    # Run on Qollab's IBM simulator
    job = backend.run(compiled_circ, shots=runs)

    out = job.result().get_counts()
    out = rename(out)

    for key in out:
        out[key] /= runs

    return out


def main():

    statistics = run_circuit(
        1 / np.sqrt(2),
        1 / np.sqrt(2),
        pi / 4,
        pi / 2,
        pi / 3,
        1000
    )

    print("Measurement probabilities:")
    print(statistics)

    # Histogram with numerical labels
    fig, ax = plt.subplots(figsize=(8, 5))

    labels = list(statistics.keys())
    probabilities = list(statistics.values())

    bars = ax.bar(labels, probabilities)

    for bar, probability in zip(bars, probabilities):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            probability + 0.01,
            f"{probability:.4f}",
            ha="center",
            va="bottom"
        )

    ax.set_xlabel("First measurement, second measurement")
    ax.set_ylabel("Probability")
    ax.set_title("Measurement statistics")
    ax.set_ylim(0, max(probabilities) * 1.2 + 0.02)

    plt.tight_layout()
    plt.show()
