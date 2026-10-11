from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate
from math import pi
import numpy as np
import matplotlib.pyplot as plt


def rename(counts):
    new_counts = {"0": 0, "1": 0, "2": 0, "3": 0}

    trans = {
        "00": "0",
        "10": "2",
        "01": "1",
        "11": "3"
    }

    for key, value in counts.items():
        new_counts[trans[key]] = value

    return new_counts


# Rotation between states |0> and |1>
def R01(theta):
    return np.array([
        [np.cos(theta), -np.sin(theta), 0, 0],
        [np.sin(theta),  np.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ], dtype=complex)


# Rotation between states |1> and |2>
def R12(theta):
    return np.array([
        [1, 0, 0, 0],
        [0, np.cos(theta), -np.sin(theta), 0],
        [0, np.sin(theta),  np.cos(theta), 0],
        [0, 0, 0, 1]
    ], dtype=complex)


def run_circuit(alpha, beta, theta0, theta1, chi, runs):

    circ = QuantumCircuit(3, 2)

    # Initialize alpha|0> + beta|1> on qubit 0
    if not np.isclose(abs(alpha)**2 + abs(beta)**2, 1):
        raise ValueError("Initial state must be normalized")

    theta = 2 * np.arctan2(abs(beta), abs(alpha))
    phase = np.angle(beta) - np.angle(alpha)

    circ.ry(theta, 0)
    circ.p(phase, 0)

    # First rotation R01
    circ.unitary(R01(theta0), [1, 2], label="R01")

    # Controlled R01
    controlled_R01 = UnitaryGate(
        R01(theta1 - theta0),
        label="R01"
    ).control(1)

    circ.append(controlled_R01, [0, 1, 2])

    # Hadamard on qubit 0
    circ.h(0)

    # Second rotation R12
    circ.unitary(R12(chi), [1, 2], label="R12")

    # Controlled R12
    controlled_R12 = UnitaryGate(
        R12(pi / 2 - 2 * chi),
        label="R12"
    ).control(1)

    circ.append(controlled_R12, [0, 1, 2])

    # Hadamard on qubit 0
    circ.h(0)

    # Final measurements
    circ.measure(1, 0)
    circ.measure(2, 1)

    # Display graphical quantum circuit
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
        -1 / np.sqrt(2),
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

    ax.set_xlabel("Measurement outcome")
    ax.set_ylabel("Probability")
    ax.set_title("Measurement statistics")
    ax.set_ylim(0, max(probabilities) * 1.2 + 0.02)

    plt.tight_layout()
    plt.show()
