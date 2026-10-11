from qiskit import QuantumCircuit
from math import pi
import numpy as np
import matplotlib.pyplot as plt

def rename(counts):
    new_counts = {"00": 0, "01": 0, "10": 0, "11": 0}
    trans = {"00": "00", "10": "01", "01": "10", "11": "11"}

    for key, value in counts.items():
        new_counts[trans[key]] = value

    return new_counts

def run_circuit(alpha, beta, xi, phi_p, phi_m,
                kappa_2, zeta, omega_p, omega_m,
                kappa_m2, runs):

    circ = QuantumCircuit(2, 2)

    # Prepare alpha|0> + beta|1> on qubit 0
    if not np.isclose(abs(alpha)**2 + abs(beta)**2, 1):
        raise ValueError("Initial state must be normalized")

    theta = 2 * np.arctan2(abs(beta), abs(alpha))
    phase = np.angle(beta) - np.angle(alpha)

    circ.ry(theta, 0)
    circ.p(phase, 0)

    # First part of the circuit
    circ.ry(-2 * xi, 0)
    circ.ry(2 * phi_p, 1)
    circ.cry(2 * phi_m - 2 * phi_p, 0, 1)
    circ.ry(2 * xi, 0)

    # Mid-circuit measurement
    circ.measure(1, 0)

    # Classical feed-forward after first measurement
    with circ.if_test((circ.clbits[0], 1)):
        circ.ry(2 * kappa_2, 0)
        circ.ry(-2 * zeta, 0)
        circ.ry(2 * omega_p, 1)
        circ.cry(2 * omega_m - 2 * omega_p, 0, 1)
        circ.ry(2 * zeta, 0)

    # Second measurement
    circ.measure(1, 1)

    # Classical feed-forward after second measurement
    with circ.if_test((circ.clbits[1], 1)):
        circ.ry(2 * kappa_m2, 0)

    # Show graphical circuit in Data Results
    circuit_fig = circ.draw(output="mpl", fold=-1)
    circuit_fig.set_size_inches(15, 5)
    plt.show()

    # Run on Qollab's selected IBM backend
    job = backend.run(circ, shots=runs)

    out = job.result().get_counts()
    out = rename(out)

    for key in out:
        out[key] /= runs

    return out

def main():
    statistics = run_circuit(
        1 / np.sqrt(2), -1 / np.sqrt(2),
        pi / 4, pi / 2, pi / 3,
        pi / 4, pi / 2, pi / 3,
        pi / 4, pi / 3,
        10000
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
