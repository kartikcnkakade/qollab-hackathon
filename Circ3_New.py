
from qiskit import QuantumCircuit
from math import pi
import numpy as np
import matplotlib.pyplot as plt


# Doubly controlled Ry using only CRY and CX gates
def ccry(circ, angle, control1, control2, target):
    circ.cry(angle / 2, control2, target)
    circ.cx(control1, control2)
    circ.cry(-angle / 2, control2, target)
    circ.cx(control1, control2)
    circ.cry(angle / 2, control1, target)


# Swap the bit order in Qiskit's output
def rename(counts):
    new_counts = {"00": 0, "01": 0, "10": 0, "11": 0}
    for key, value in counts.items():
        new_counts[key[::-1]] = value
    return new_counts


def run_circuit(alpha, beta, xi, phi_p, phi_m,
                kappa_2, zeta, omega_p, omega_m,
                kappa_m2, runs):

    # q0, q1 = original qubits
    # q2 = ancilla recording the first measurement
    circ = QuantumCircuit(3, 2)

    # Prepare alpha|0> + beta|1> on q0
    if not np.isclose(abs(alpha)**2 + abs(beta)**2, 1):
        raise ValueError("Initial state must be normalized")

    theta = 2 * np.arctan2(abs(beta), abs(alpha))
    phase = np.angle(beta) - np.angle(alpha)

    circ.ry(theta, 0)
    circ.p(phase, 0)

    # Original first part
    circ.ry(-2 * xi, 0)
    circ.ry(2 * phi_p, 1)
    circ.cry(2 * phi_m - 2 * phi_p, 0, 1)
    circ.ry(2 * xi, 0)

    # Record the first measurement coherently in q2
    circ.cx(1, 2)

    # Operations conditioned on first outcome = 1
    circ.cry(2 * kappa_2, 2, 0)
    circ.cry(-2 * zeta, 2, 0)
    circ.cry(2 * omega_p, 2, 1)

    ccry(circ, 2 * omega_m - 2 * omega_p, 0, 2, 1)

    circ.cry(2 * zeta, 2, 0)

    # Final conditional rotation, before deferred measurement
    circ.cry(2 * kappa_m2, 1, 0)

    # Both measurements at the end
    circ.measure(2, 0)
    circ.measure(1, 1)

    print(circ.draw("text"))

    # Qollab automatically provides backend
    job = backend.run(circ, shots=runs)
    counts = job.result().get_counts()

    out = rename(counts)

    for key in out:
        out[key] /= runs

    return out


statistics = run_circuit(
    1 / np.sqrt(2), -1 / np.sqrt(2),
    pi / 4, pi / 2, pi / 3,
    pi / 4, pi / 2, pi / 3,
    pi / 4, pi / 3,
    10000
)

print("Measurement probabilities:")
print(statistics)

fig, ax = plt.subplots()

bars = ax.bar(list(statistics.keys()), list(statistics.values()))

ax.bar_label(bars, fmt="%.4f", padding=3)

ax.set_xlabel("First measurement, second measurement")
ax.set_ylabel("Probability")
ax.set_title("Measurement statistics")
ax.set_ylim(0, max(statistics.values()) * 1.15)

plt.show()
