# Implementing a Quantum Instrument with Adaptive Circuits

### Adaptivity reduces the ancilla dimension needed for a three-outcome qubit instrument

---

## 1. Quantum instruments and mid-circuit measurements on NISQ devices

### 1.1 Quantum instruments

A measurement on a quantum computer gives two things: a classical outcome and a post-measurement quantum state. A **quantum instrument** describes both at once.

An instrument with outcome set $\Omega$, from a system with Hilbert space $\mathcal H$ to $\mathcal K$, is a family of completely positive maps $\{\mathcal T_k\}_{k\in\Omega}$ whose sum is trace preserving:

$$
\mathcal T_k(\rho)=\sum_m K_{k,m}\,\rho\,K_{k,m}^\dagger,\qquad \sum_{k}\sum_m K_{k,m}^\dagger K_{k,m}=\mathbb 1 .
$$

- Outcome $k$ occurs with probability $\Pr(k)=\operatorname{Tr}[\mathcal T_k(\rho)]$.
- After outcome $k$, the system is in the state $\mathcal T_k(\rho)/\Pr(k)$.
- The **induced POVM** $E_k=\sum_m K_{k,m}^\dagger K_{k,m}$ gives the outcome statistics alone: $\Pr(k)=\operatorname{Tr}[E_k\rho]$.
- The **induced channel** $\sum_k\mathcal T_k$ describes the state change when the outcome is ignored.
- The **Lüders instrument** of a POVM $\{B_j\}$ is $\mathcal J_j(\rho)=\sqrt{B_j}\,\rho\,\sqrt{B_j}$. It is the least disturbing way to measure $\{B_j\}$.

Different instruments can share the same POVM but act differently on the state. To say two implementations realise "the same instrument", every map $\mathcal T_k$ must agree, not just the probabilities.

Special cases: a single outcome gives a channel, and a one-dimensional output space gives a POVM.

### 1.2 Mid-circuit measurements as quantum instruments

On current NISQ processors (superconducting qubits, trapped ions), a **mid-circuit measurement** reads out a subset of qubits partway through a circuit. The unmeasured qubits keep evolving coherently, and the measured qubit can be **reset** and reused. These devices have tens to a few hundred noisy qubits, limited connectivity, two-qubit gate errors around $10^{-3}$ to $10^{-2}$, and readout errors around $10^{-2}$. Mid-circuit measurement and reset are therefore mainly used to save qubits: a few ancillas are measured and reused many times.

The instrument formalism describes exactly what such an operation does.

**(a) Measuring one qubit.** A projective measurement of an ancilla qubit $A$ in the computational basis is an instrument on the joint system $SA$ with outcomes $r\in\{0,1\}$:

$$
\mathcal M_r(\sigma)=(\mathbb 1_S\otimes\ket r\!\bra r_A)\,\sigma\,(\mathbb 1_S\otimes\ket r\!\bra r_A).
$$

This is the Lüders instrument of the projective measurement $\{\mathbb 1\otimes\ket r\!\bra r\}$. If the measured qubit is discarded, the destructive version $\mathcal M_r(\sigma)=\operatorname{Tr}_A[(\mathbb 1\otimes\ket r\!\bra r)\sigma]$ maps $SA$ to $S$. A reset afterwards turns this into $\sigma\mapsto\operatorname{Tr}_A[(\mathbb 1\otimes\ket r\!\bra r)\sigma]\otimes\ket0\!\bra0_A$, so the ancilla returns to $\ket0$ and the record $r$ is kept.

**(b) What the system sees.** Prepare the ancilla in $\ket0$, apply a unitary $U$ to $SA$, then measure $A$. The operation induced on the system alone is

$$
\mathcal T_r(\rho)=\operatorname{Tr}_A\!\big[(\mathbb 1\otimes\ket r\!\bra r)\,U(\rho\otimes\ket0\!\bra0)U^\dagger\big]=K_r\,\rho\,K_r^\dagger,\qquad K_r=(\mathbb 1\otimes\bra r)\,U\,(\mathbb 1\otimes\ket0).
$$

This is a quantum instrument on $S$, with Kraus operators given by blocks of $U$. Conversely, every instrument arises this way (Ozawa's dilation theorem). The ancilla readout must have at least $r_{\mathcal T}=\sum_k\operatorname{rank}\mathcal T_k$ orthogonal outcomes, so with $n_A$ ancilla qubits measured once, $2^{n_A}\ge r_{\mathcal T}$. All schemes below are built and checked by computing these blocks $\bra r\,U\ket0$.

**(c) Several readouts.** Measuring several ancilla qubits gives a tuple of outcomes. The instrument is the product of single-qubit measurements, and any relabelling of the tuple to the desired outcome $k$ is a classical post-processing of instruments. Scheme A uses $k=2a_1+a_2$.

**(d) Feed-forward.** If the next gates depend on the recorded outcome, the overall process is a composition of instruments chosen by the previous outcome. This is an adaptive sequence of instruments, defined in Section 2.

**(e) Noise.** On hardware the formalism still applies, with noisy instruments:

- A readout error with assignment probabilities $P(b\,|\,r)$ turns $\mathcal M_r$ into $\mathcal M'_b=\sum_r P(b|r)\,\mathcal M_r$. The wrong recorded bit also steers the feed-forward wrongly.
- Unmeasured qubits decohere while the readout and the classical decision take place, typically about 1 µs.
- Non-ideal (non-QND) measurements disturb the measured qubit beyond the ideal projection.

Comparing implementations "as instruments" therefore means comparing the full set of operations $\{\mathcal T_k\}$ (their Choi matrices), not only the outcome statistics.

---

## 2. Adaptive circuits and what this project does

In an **adaptive (dynamic) circuit**, the result of a mid-circuit measurement decides, in real time, which gates come next. This is called classical feed-forward.

Mathematically, this is an **adaptive sequence of instruments (ASI)** (Def. 4 of Ref. [1]): at step $n$ an instrument $\mathcal I^{\,n,a_{n-1}}$ is chosen according to the previous outcome $a_{n-1}$. The total instrument is

$$
\mathcal T_{a_N}=\sum_{a_1,\dots,a_{N-1}}\mathcal I^{\,N,a_{N-1}}_{a_N}\circ\cdots\circ\mathcal I^{\,1}_{a_1}.
$$

Each step needs only enough ancilla to resolve its own outcomes. This trades **space** (ancilla dimension, number of qubits) for **time** (number of measurement rounds and feed-forward steps).

**This project**

- Takes one concrete, genuinely non-classical three-outcome qubit instrument $\mathcal T$ with Kraus rank 3.
- Implements it in three ways:
  - **Scheme A:** one-shot dilation with a qutrit ancilla (6 × 6 unitary). Also built on qubit hardware with a two-qubit ancilla for Qiskit.
  - **Scheme B:** adaptive circuit with one ancilla qubit and two binary measurements.
  - **Scheme C:** adaptive circuit built from **Theorem 1** of Ref. [1], where the first step is a Lüders instrument.
- Verifies all of them:
  - exactly, with SymPy;
  - as Qiskit dynamic circuits on the Aer simulator, with $10^6$ shots;
  - by instrument tomography.
- Compares their resources. Adaptivity lowers the ancilla dimension from 3 (4 on qubit hardware) to 2. Theorem 1 is one systematic way to get such a circuit, but it is not always the cheapest one.

---

## 3. The instrument $\mathcal T$

### 3.1 Parameters and notation

The system is one qubit, $\mathcal H=\mathcal K=\mathbb C^2$. There are three parameters $\theta_0,\theta_1,\chi\in(0,\tfrac\pi2)$:

$$
c_k=\cos\theta_k,\quad s_k=\sin\theta_k,\qquad p=\tfrac12(\cos\chi+\sin\chi),\quad m=\tfrac12(\cos\chi-\sin\chi),
$$

with the identities $p^2+m^2=\tfrac12$ and $2pm=\tfrac12\cos2\chi$.

Further notation:

- $P_\pm=\ket\pm\!\bra\pm$.
- $R(\alpha)=\begin{pmatrix}\cos\alpha&-\sin\alpha\\ \sin\alpha&\cos\alpha\end{pmatrix}=R_Y(2\alpha)$, with $R_Y(\varphi)=e^{-i\varphi Y/2}$, which is Qiskit's `ry`.

### 3.2 Kraus operators

$\mathcal T_k(\rho)=K_k\rho K_k^\dagger$ for $k=0,1,2$, with

$$
K_0=\begin{pmatrix}c_0&0\\0&c_1\end{pmatrix},\qquad
K_1=\begin{pmatrix}p\,s_0&m\,s_1\\ m\,s_0&p\,s_1\end{pmatrix},\qquad
K_2=\begin{pmatrix}p\,s_0&-m\,s_1\\ -m\,s_0&p\,s_1\end{pmatrix}.
$$

**Structure.** Write $S=\operatorname{diag}(s_0,s_1)$, $F_1=\cos\chi\,P_++\sin\chi\,P_-$ and $F_2=\sin\chi\,P_++\cos\chi\,P_-$. In matrix form

$$
F_1=\begin{pmatrix}p&m\\m&p\end{pmatrix},\quad F_2=\begin{pmatrix}p&-m\\-m&p\end{pmatrix},\qquad K_1=F_1S,\qquad K_2=F_2S .
$$

So $\mathcal T$ works as follows:

1. A weak $Z$ measurement with Kraus operators $\{K_0,S\}$, where $K_0^\dagger K_0+S^\dagger S=\mathbb 1$.
2. On the "not 0" branch only, a weak $X$ measurement with Kraus operators $\{F_1,F_2\}$, where $F_1^\dagger F_1+F_2^\dagger F_2=\mathbb 1$.

Each operation has one Kraus operator, so the **Kraus rank is $r_{\mathcal T}=3$**. $K_1$ and $K_2$ do not commute for generic parameters, so $\mathcal T$ is not a relabelled projective measurement.

### 3.3 Induced POVM

$$
E_0=K_0^\dagger K_0=\begin{pmatrix}c_0^2&0\\0&c_1^2\end{pmatrix},\qquad
E_1=\frac12\begin{pmatrix}s_0^2&s_0s_1\cos2\chi\\ s_0s_1\cos2\chi&s_1^2\end{pmatrix},\qquad
E_2=\frac12\begin{pmatrix}s_0^2&-s_0s_1\cos2\chi\\ -s_0s_1\cos2\chi&s_1^2\end{pmatrix}.
$$

For $E_1$: $K_1^\dagger K_1$ has diagonal $(p^2+m^2)s_k^2$ and off-diagonal $2pm\,s_0s_1$. $E_2$ is the same with $-2pm$.

**Completeness:** $E_0+E_1+E_2=\operatorname{diag}(c_0^2+s_0^2,\ c_1^2+s_1^2)=\mathbb 1$.

For $\ket\psi=\psi_0\ket0+\psi_1\ket1$:

$$
\Pr(0)=c_0^2|\psi_0|^2+c_1^2|\psi_1|^2,\qquad
\Pr(1)\text{ and }\Pr(2)=\tfrac12\big(s_0^2|\psi_0|^2+s_1^2|\psi_1|^2\big)\pm s_0s_1\cos2\chi\,\operatorname{Re}(\psi_0^*\psi_1),
$$

with the upper sign for $k=1$.

**Example.** At $\theta_0=0.37$, $\theta_1=1.11$, $\chi=0.52$:

$$
K_0=\begin{pmatrix}0.9323&0\\0&0.4447\end{pmatrix},\ 
K_1=\begin{pmatrix}0.2467&0.1661\\0.0671&0.6112\end{pmatrix},\ 
K_2=\begin{pmatrix}0.2467&-0.1661\\-0.0671&0.6112\end{pmatrix},
$$

$$
E_0=\begin{pmatrix}0.8692&0\\0&0.1977\end{pmatrix},\ 
E_1=\begin{pmatrix}0.0654&0.0820\\0.0820&0.4011\end{pmatrix},\ 
E_2=\begin{pmatrix}0.0654&-0.0820\\-0.0820&0.4011\end{pmatrix}.
$$

---

## 4. The circuits

Throughout, the system qubit is `s`. A unitary acting on "system ⊗ ancilla" is written in the basis $\ket{s\,a}$ with the system first.

### 4.1 Scheme A: qutrit ancilla, one 6 × 6 unitary

**Goal:** $U\ket\psi\ket0_a=K_0\ket\psi\ket0_a+K_1\ket\psi\ket1_a+K_2\ket\psi\ket2_a$, followed by one measurement of the qutrit.

**Construction.** $U=V_2V_1$, where $R^{jk}(\alpha)$ is a real rotation between qutrit levels $j$ and $k$:

$$
V_1=\ket0\!\bra0\otimes R^{01}(\theta_0)+\ket1\!\bra1\otimes R^{01}(\theta_1),\qquad
V_2=(H\otimes\mathbb 1)\big[\ket0\!\bra0\otimes R^{12}(\chi)+\ket1\!\bra1\otimes R^{12}(\tfrac\pi2-\chi)\big](H\otimes\mathbb 1).
$$

- $V_1$ produces $K_0\ket\psi\ket0+S\ket\psi\ket1$: the weak $Z$ step.
- $V_2$ splits the $\ket1$ branch into $F_1S\ket\psi\ket1+F_2S\ket\psi\ket2$: the weak $X$ step.

```
                 ┌──────────────┐       ┌──────────────────┐
system |ψ⟩ ──────┤      ◐       ├── H ──┤        ◐         ├── H ──────────────  K_k|ψ⟩
                 │      │       │       │        │         │
                 │ ┌────┴─────┐ │       │ ┌──────┴───────┐ │      ┌───┐
qutrit |0⟩ ══════╪═┤ R⁰¹(θ₀|θ₁)├═╪═══════╪═┤R¹²(χ | π/2−χ)├═╪══════┤ M ├═══  k ∈ {0,1,2}
                 │ └──────────┘ │       │ └──────────────┘ │      └───┘
                 └──── V₁ ──────┘       └─────── V₂ ───────┘
   ◐ = multiplexed control: left angle if the system is |0⟩, right angle if |1⟩
```

**The unitary.** Basis $\ket{s\,a}$ in the order $\ket{00},\ket{01},\ket{02},\ket{10},\ket{11},\ket{12}$:

$$
U=\begin{pmatrix}
c_0 & -s_0 & 0 & 0 & 0 & 0\\
p\,s_0 & p\,c_0 & -p & m\,s_1 & m\,c_1 & m\\
p\,s_0 & p\,c_0 & p & -m\,s_1 & -m\,c_1 & m\\
0 & 0 & 0 & c_1 & -s_1 & 0\\
m\,s_0 & m\,c_0 & m & p\,s_1 & p\,c_1 & -p\\
-m\,s_0 & -m\,c_0 & m & p\,s_1 & p\,c_1 & p
\end{pmatrix},\qquad \bra k_aU\ket0_a=K_k .
$$

$U^TU=\mathbb 1_6$ holds identically, as checked with SymPy.

### 4.2 Scheme A for Qiskit: the qutrit encoded in two qubits

Qubit hardware has no qutrits. The ancilla is two qubits $a_1a_2$ (dimension 4), with the encoding

$$
\ket0\mapsto\ket{00},\qquad \ket1\mapsto\ket{01},\qquad \ket2\mapsto\ket{10},\qquad \ket{11}\ \text{unused}.
$$

**Goal:** $U\ket\psi\ket{00}=K_0\ket\psi\ket{00}+K_1\ket\psi\ket{01}+K_2\ket\psi\ket{10}$.

**Construction.** $U=V_2V_1$, where:

- $V_1=\sum_s\ket s\!\bra s\otimes\mathbb 1_{a_1}\otimes R(\theta_s)_{a_2}$.
- $V_2=P_+\otimes G(\chi)+P_-\otimes G(\tfrac\pi2-\chi)$.
- $G(\alpha)$ is a rotation by $\alpha$ inside $\operatorname{span}\{\ket{01},\ket{10}\}$ that fixes $\ket{00}$ and $\ket{11}$.

**Kraus operators of the circuit.** These are $M_{a_1a_2}=\bra{a_1a_2}U\ket{00}$:

$$
M_{00}=K_0,\quad M_{01}=K_1,\quad M_{10}=K_2,\quad M_{11}=0 .
$$

With the relabelling $k=2a_1+a_2$, the total instrument is exactly $\mathcal T$. Its induced POVM is $\{E_0,E_1,E_2\}$ plus $A_{11}=0$, so readout $11$ never occurs in the ideal circuit and works as an error flag.

**Gate-level circuit** (5 CNOTs). The angles are $a=\theta_0+\theta_1-\tfrac\pi2$, $b=\theta_0-\theta_1+\tfrac\pi2$ and $\tau=\tfrac\pi4-\chi$:

```
                   ┌──────── V₁ ────────┐ ┌──────────────────────── V₂ ─────────────────────────┐
  s : |ψ⟩ ─────────────────●─────────────── H ───────────●──────────────────────────●────────── H ───────────────── K_k|ψ⟩
                           │                             │                          │
  a1: |0⟩ ───────────────────────────────────── RY(π/4)─⊕─ RY(τ) ─⊕─ RY(−τ) ───────⊕── RY(−π/4) ──●──── M ═══ a1
                           │                                      │                              │
  a2: |0⟩ ── RY(b) ────────⊕─── RY(a) ────────────────────────────●──────────────────────────────⊕──── M ═══ a2

                                                                                     outcome k = 2·a1 + a2
```

- **$V_1$** uses one CNOT: $R_Y(a)X^sR_Y(b)\ket0=R_Y(2\theta_s)\ket0$.
- **$V_2$:**
  1. The block $R_Y(\tfrac\pi4)\,\text{CX}\,R_Y(\tau)\,\text{CX}\,R_Y(-\tau)\,\text{CX}\,R_Y(-\tfrac\pi4)$ rotates $a_1$ only if $a_2=1$. The rotation is by $R(\chi)$ when the system is in $\ket+$ and by $R(\tfrac\pi2-\chi)$ when it is in $\ket-$.
  2. The final $\text{CNOT}_{a_1\to a_2}$ maps $\ket{11}\to\ket{10}$.
- This last CNOT can be replaced by classical post-processing: drop it and compute $a_2\leftarrow a_2\oplus a_1$ after readout. That gives 4 CNOTs.

Qiskit drawing at the example parameters:

```text
                              ┌───┐                                                     ┌───┐
  s: ────────────────■────────┤ H ├────────■──────────────────────────────────────■─────┤ H ├───────────────
      ┌─────────┐    │        └───┘      ┌─┴─┐┌────────────┐┌───┐┌─────────────┐┌─┴─┐┌──┴───┴───┐        ┌─┐
 a1: ─┤ Ry(π/4) ├────┼───────────────────┤ X ├┤ Ry(0.2654) ├┤ X ├┤ Ry(-0.2654) ├┤ X ├┤ Ry(-π/4) ├──■─────┤M├
     ┌┴─────────┴─┐┌─┴─┐┌───────────────┐└───┘└────────────┘└─┬─┘└─────────────┘└───┘└──────────┘┌─┴─┐┌─┐└╥┘
 a2: ┤ Ry(0.8308) ├┤ X ├┤ Ry(-0.090796) ├─────────────────────■──────────────────────────────────┤ X ├┤M├─╫─
     └────────────┘└───┘└───────────────┘                                                        └───┘└╥┘ ║
c: 2/══════════════════════════════════════════════════════════════════════════════════════════════════╩══╩═
                                                                                                       0  1
```

**Equivalence with the qutrit scheme.** The qutrit Kraus operators $\bra i U_6\ket0$ equal the circuit's $M_{\mathrm{enc}(i)}$, signs included. The POVMs and the Choi matrices of every operation also coincide. The dilations themselves are identical too: $(\mathbb 1\otimes\mathrm{enc})\,U_6\ket{\psi,0}=U_{\rm circ}\ket{\psi,00}$. This was checked exactly on 2000 random parameter sets, and confirmed in Aer with $10^6$ shots: probabilities within sampling error, post-measurement fidelity $1-10^{-12}$, and readout 11 never observed.

### 4.3 Scheme B: one ancilla qubit, two binary measurements, feed-forward

**Idea.** Ask two yes/no questions instead of one three-way question. First "is the outcome 0?", and only if not, "is it 1 or 2?".

**Instruments.**

- Step 1: $\{K_0,\ S\}$.
- Step 2, run only if the first answer was "not 0": $\{F_1,\ F_2\}$.
- Paths: $K_0$, $F_1S=K_1$, $F_2S=K_2$.

**Unitaries.** Basis $\ket{s\,a}=\ket{00},\ket{01},\ket{10},\ket{11}$:

$$
W_1=\ket0\!\bra0\otimes R(\theta_0)+\ket1\!\bra1\otimes R(\theta_1)=\begin{pmatrix}c_0&-s_0&0&0\\ s_0&c_0&0&0\\ 0&0&c_1&-s_1\\ 0&0&s_1&c_1\end{pmatrix},
$$

$$
W_2=P_+\otimes R(\chi)+P_-\otimes R(\tfrac\pi2-\chi)=\begin{pmatrix}p&-p&m&m\\ p&p&-m&m\\ m&m&p&-p\\ -m&m&p&p\end{pmatrix}.
$$

Here $\bra0_aW_1\ket0_a=K_0$, $\bra1_aW_1\ket0_a=S$, $\bra0_aW_2\ket0_a=F_1$ and $\bra1_aW_2\ket0_a=F_2$.

**Circuit.**

```
  s : |ψ⟩ ──────────●──────────────────────────── H ─────────────●──── H ──────────────  K_k|ψ⟩
                    │                             ┊              │
  a : |0⟩ ── RY(b) ─⊕─ RY(a) ── M ════╗           ┊ RY(2χ−π) ────⊕──────── M ═══ m₂
                                      ║           ┊                ┊
  m₁  ════════════════════════════════╩═══════ [ if m₁ = 1 : run this block ] ══════

  outcome:  k = 0 if m₁ = 0;   k = 1 + m₂ if m₁ = 1
```

- Step 1 is the same one-CNOT multiplexer as in Scheme A: $a=\theta_0+\theta_1-\tfrac\pi2$, $b=\theta_0-\theta_1+\tfrac\pi2$.
- **Reset is absorbed.** When $m_1=1$ the ancilla is in $\ket1$, and $R_Y(2\chi-\pi)\ket1=R_Y(2\chi)\ket0$. So the conditional block $\{H,\ R_Y(2\chi-\pi),\ \text{CX},\ H\}$ implements $W_2$ with no explicit reset.

Qiskit drawing:

```text
                                             ┌──────      ┌───┐          ┌───┐ ───────┐
   s: ────────────────■──────────────────────┤       ─────┤ H ├───────■──┤ H ├        ├─
      ┌────────────┐┌─┴─┐┌───────────────┐┌─┐│ If-0  ┌────┴───┴────┐┌─┴─┐└┬─┬┘  End-0 │
   a: ┤ Ry(0.8308) ├┤ X ├┤ Ry(-0.090796) ├┤M├┤       ┤ Ry(-2.1016) ├┤ X ├─┤M├─        ├─
      └────────────┘└───┘└───────────────┘└╥┘└──╥─── └─────────────┘└───┘ └╥┘  ───────┘
                                           ║ ┌──╨──┐                       ║
m1: 1/═════════════════════════════════════╩═╡ 0x1 ╞═══════════════════════╬════════════
                                           0 └─────┘                       ║
m2: 1/═════════════════════════════════════════════════════════════════════╩════════════
```

### 4.4 Scheme C: Theorem 1 of "Sequential realization of quantum instruments"

Scheme C is built directly from Theorem 1 of Ref. [1], together with the notions and corollaries around it. We first state what we use, in the notation of Ref. [1], and then apply it to $\mathcal T$.

#### 4.4.1 Results of Ref. [1] used here

**Notation (Defs. 1, 2 and 4–6 of Ref. [1]).**

- $\mathrm{Ins}(\Omega,\mathcal H,\mathcal K)$ is the set of instruments with outcome set $\Omega$ from $\mathcal L(\mathcal H)$ to $\mathcal L(\mathcal K)$. $A^{\mathcal T}$ is the induced POVM of $\mathcal T$.
- **Post-processing (Def. 1).** A POVM $\mathsf B\in\mathcal O(\Omega_{\mathsf B},\mathcal H)$ is a post-processing of $A^{\mathcal T}$, written $A^{\mathcal T}\to\mathsf B$, if $\mathsf B_j=\sum_{k\in\Omega}\nu_{kj}A^{\mathcal T}_k$ for a stochastic matrix $\nu_{kj}\ge0$ with $\sum_j\nu_{kj}=1$. If all $\nu_{kj}\in\{0,1\}$ it is a **coarse-graining**: outcomes are merged into groups.
- **Lüders instrument (Def. 2):** $\mathcal I^{\mathsf B}_j(\rho)=\sqrt{\mathsf B_j}\,\rho\,\sqrt{\mathsf B_j}$.
- **Adaptive sequence (Defs. 4–6).** A 2-step adaptive sequence consists of an **initial instrument** $\mathcal J\in\mathrm{Ins}(\Omega_{\mathsf B},\mathcal H,\mathcal H)$ and **residual instruments** $\mathcal R^j\in\mathrm{Ins}(\Omega,\mathcal H,\mathcal K)$, one per intermediate outcome $j$. Its **total instrument** is $\mathcal T_k=\sum_j\mathcal R^j_k\circ\mathcal J_j$.

> **Theorem 1 (Ref. [1]).** Given a total instrument $\mathcal T\in\mathrm{Ins}(\Omega,\mathcal H,\mathcal K)$ and any POVM $\mathsf B\in\mathcal O(\Omega_{\mathsf B},\mathcal H)$ such that the induced POVM of $\mathcal T$ can be post-processed to $\mathsf B$ ($A^{\mathcal T}\to\mathsf B$), there exist an initial instrument $\mathcal J\in\mathrm{Ins}(\Omega_{\mathsf B},\mathcal H,\mathcal H)$ and residual instruments $\mathcal R^j\in\mathrm{Ins}(\Omega,\mathcal H,\mathcal K)$, $j\in\Omega_{\mathsf B}$, such that
> $$\mathcal T_k=\sum_{j\in\Omega_{\mathsf B}}\mathcal R^j_k\circ\mathcal J_j\qquad\forall k\in\Omega,$$
> where the initial instrument $\mathcal J_j(\rho)=\sqrt{\mathsf B_j}\,\rho\,\sqrt{\mathsf B_j}$ is the Lüders instrument of $\mathsf B$.

**The construction in the proof** is what we implement.

- Let $\{T'_{k,m}\}$ be Kraus operators of $\mathcal T_k$. Define $T_{k,jm}=\sqrt{\nu_{kj}}\,T'_{k,m}$ (Eq. (10)). Then $\mathsf B_j=\sum_{k,m}T_{k,jm}^\dagger T_{k,jm}$, and the input range of each $T_{k,jm}$ lies in the support of $\mathsf B_j$, with projector $\Pi_j$.
- For $\dim\mathcal H\le\dim\mathcal K$, the residual Kraus operators are (Eq. (18))
  $$R^j_{k,m}=T_{k,jm}\,\mathsf B_j^{-1/2}+c_{k,jm}V_{k,jm}(\mathbb 1-\Pi_j),\qquad \sum_{k,m}|c_{k,jm}|^2=1,$$
  where $\mathsf B_j^{-1/2}$ is the Moore–Penrose pseudo-inverse of $\sqrt{\mathsf B_j}$ and $V_{k,jm}$ is an isometry.
- Normalisation holds by $\sum_{k,m}R^{j\dagger}_{k,m}R^j_{k,m}=\Pi_j+(\mathbb 1-\Pi_j)=\mathbb 1$ (Eq. (19)).
- Correctness holds by $R^j_{k,m}\sqrt{\mathsf B_j}=T_{k,jm}\Pi_j=T_{k,jm}$ (Eq. (20)).

**Other results of Ref. [1] used in Scheme C.**

- **Corollary 1.1.** If $\dim\mathcal H\le\dim\mathcal K$, the Kraus rank of $\mathcal R^j$ is $m_j=\sum_{k:\,\nu_{kj}>0}\operatorname{rank}\mathcal T_k$. This fixes how large the ancilla must be in step 2.
- **Remark 3.** Paths $(j,k)$ with $\nu_{kj}=0$ never occur. After outcome $j$, the residual instrument can only produce outcomes $k$ in group $j$.
- **Remark 2.** The initial instrument depends only on $\mathsf B$, not on how $\mathcal T$ changes the state.
- **Ancilla dimension, Eq. (28).** If one $d_A$-dimensional ancilla is reset and reused for $\mathcal J$ and all $\mathcal R^j$, it needs
  $$d_A=\max\Big(g\Big\lceil\frac{|\Omega_{\mathsf B}|}{g}\Big\rceil,\ \{g\,m_j\}_{j\in\Omega_{\mathsf B}}\Big),\qquad g=\Big\lceil\frac{\dim\mathcal K}{\dim\mathcal H}\Big\rceil .$$
- **Theorem 2.** For $\mathcal T$ with Kraus rank $r_{\mathcal T}$ and $\dim\mathcal H\le\dim\mathcal K$, the minimal ancilla dimension of a 2-step adaptive sequence obtained from Theorem 1 (with classical post-processing of outcomes) is
  $$d_A=g\Big\lceil\sqrt{r_{\mathcal T}/g}\,\Big\rceil .$$
  It is reached by a coarse-graining into groups of at most $d_A/g$ rank-one operations.

#### 4.4.2 How Scheme C uses these results

For $\mathcal T$: $\mathcal H=\mathcal K=\mathbb C^2$, so $g=1$. Each $\mathcal T_k$ has one Kraus operator $T'_{k,1}=K_k$, so $r_{\mathcal T}=3$. All $\mathsf B_j$ below are full rank, so $\Pi_j=\mathbb 1$ and the kernel term in Eq. (18) vanishes.

1. **Target ancilla size (Theorem 2).** $d_A=\lceil\sqrt3\,\rceil=2$: one reusable ancilla qubit.
2. **Choose $\mathsf B$ (Def. 1).** Use the coarse-graining that merges outcomes $\{0,1\}$ and keeps $\{2\}$:
   $$\nu=\begin{pmatrix}1&0\\1&0\\0&1\end{pmatrix},\qquad \mathsf B_0=E_0+E_1=\mathbb 1-E_2,\qquad \mathsf B_1=E_2 .$$
3. **Check the ancilla budget (Corollary 1.1 and Eq. (28)).** We have $m_0=\operatorname{rank}\mathcal T_0+\operatorname{rank}\mathcal T_1=2$ and $m_1=\operatorname{rank}\mathcal T_2=1$, so $d_A=\max(|\Omega_{\mathsf B}|,m_0,m_1)=\max(2,2,1)=2$. This reaches the Theorem 2 bound.
4. **Step 1 = initial instrument (Theorem 1).** The Lüders instrument of $\mathsf B$, with Kraus operators $\sqrt{\mathsf B_0}$ and $\sqrt{\mathsf B_1}$, measured with the ancilla qubit, giving $m_1=j$.
5. **Step 2 = residual instruments (Eq. (18)).** $R^j_k=\sqrt{\nu_{kj}}\,K_k\,\mathsf B_j^{-1/2}$:
   - $j=0$: $R^0_0=L_0=K_0\mathsf B_0^{-1/2}$, $R^0_1=L_1=K_1\mathsf B_0^{-1/2}$, $R^0_2=0$. Two outcomes, so the same ancilla qubit is reused.
   - $j=1$: $R^1_2=V=K_2\mathsf B_1^{-1/2}$, $R^1_0=R^1_1=0$. One Kraus operator, which must be unitary, so no ancilla is needed.
6. **Outcome logic (Remark 3).** After $m_1=1$ the outcome is 2 with certainty. After $m_1=0$ it is $m_2\in\{0,1\}$.

The rest of this section turns each of these instruments into a circuit.

**Step 1: Lüders instrument of $\mathsf B$.** $\mathsf B_0$ and $\mathsf B_1$ commute. Diagonalise them with $R(\xi)$:

$$
\mathsf B_1=R(\xi)\operatorname{diag}(\lambda_+,\lambda_-)R(\xi)^T,\quad
\xi=\tfrac12\operatorname{atan2}(-2s_0s_1\cos2\chi,\ s_0^2-s_1^2),
$$

$$
\lambda_\pm=\tfrac14\Big[s_0^2+s_1^2\pm\sqrt{(s_0^2-s_1^2)^2+4s_0^2s_1^2\cos^22\chi}\Big],\qquad \sin\varphi_\pm=\sqrt{\lambda_\pm}.
$$

With $\mathsf J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$, the step-1 unitary is

$$
U_J=\sqrt{\mathsf B_0}\otimes\mathbb 1+\sqrt{\mathsf B_1}\otimes\mathsf J
=(R(\xi)\otimes\mathbb 1)\big[\ket0\!\bra0\otimes R(\varphi_+)+\ket1\!\bra1\otimes R(\varphi_-)\big](R(\xi)^T\otimes\mathbb 1),
$$

so that $\bra j_aU_J\ket0_a=\sqrt{\mathsf B_j}$.

**Step 2, branch $m_1=1$ (outcome 2).** The single Kraus operator $V=K_2\mathsf B_1^{-1/2}$ is unitary. It is the polar part of $K_2$:

$$
V=R(\kappa_2),\qquad \kappa_2=\operatorname{atan2}\big(m(s_1-s_0),\ p(s_0+s_1)\big).
$$

This branch needs **no ancilla and no measurement**.

**Step 2, branch $m_1=0$ (outcomes 0, 1).** The Kraus operators are $L_0=K_0\mathsf B_0^{-1/2}$ and $L_1=K_1\mathsf B_0^{-1/2}$, with $L_0^TL_0+L_1^TL_1=\mathbb 1$. Write them in polar form, $L_k=R(\kappa_k)\sqrt{M_k}$, with $M_0=L_0^TL_0=R(\zeta)\operatorname{diag}(\mu_+,\mu_-)R(\zeta)^T$ and $\cos\omega_\pm=\sqrt{\mu_\pm}$. Then

$$
U_R=\big(R(\kappa_0)\otimes\ket0\!\bra0+R(\kappa_1)\otimes\ket1\!\bra1\big)\big(\sqrt{M_0}\otimes\mathbb 1+\sqrt{M_1}\otimes\mathsf J\big).
$$

The ancilla-controlled rotation is followed directly by a measurement, so it is applied classically as $R_Y(2\kappa_{m_2})$ after reading $m_2$. No reset is needed, because in this branch the ancilla was just measured as $\ket0$.

**Paths:**

- $(m_1,m_2)=(0,0)$: $L_0\sqrt{\mathsf B_0}=K_0$.
- $(0,1)$: $L_1\sqrt{\mathsf B_0}=K_1$.
- $(1,-)$: $V\sqrt{\mathsf B_1}=K_2$.

Combinations such as "outcome 2 in branch 0" have zero probability, as stated in Remark 3 of Ref. [1].

**Circuit.** The eigenbasis rotations are merged into neighbouring single-qubit gates. With $\alpha=\varphi_++\varphi_--\tfrac\pi2$, $\beta=\varphi_+-\varphi_-+\tfrac\pi2$, $\alpha'=\omega_++\omega_--\tfrac\pi2$ and $\beta'=\omega_+-\omega_-+\tfrac\pi2$:

```
             ┌──────────── step 1: U_J (Lüders) ────────────┐
  s : |ψ⟩ ── RY(−2ξ) ──────●──────────────────────────────────── [m₁=1]: RY(2ξ+2κ₂) ────────────────────────────────────  K_k|ψ⟩
                           │                                     [m₁=0]: RY(2ξ−2ζ) ──●─────────────────── [m₂]: RY(2ζ+2κ_{m₂})
                           │                                                         │
  a : |0⟩ ── RY(β) ────────⊕── RY(α) ── M ═══ m₁                [m₁=0]: RY(β′) ──────⊕── RY(α′) ── M ═══ m₂

  outcome:  k = 2 if m₁ = 1;   k = m₂ if m₁ = 0
```

Qiskit drawing at the example parameters:

```text
      ┌────────────┐                        ┌────── ┌─────────────┐┌──────── ┌─────────────┐                        ┌────── ┌──────────────┐┌──────── ┌──────────────┐ ───────┐  ───────┐
   s: ┤ Ry(2.6873) ├──■─────────────────────┤       ┤ Ry(-2.4574) ├┤         ┤ Ry(-2.4475) ├──■─────────────────────┤ If-1  ┤ Ry(-0.19975) ├┤ Else-1  ┤ Ry(-0.56375) ├  End-1 ├─        ├─
      ├────────────┤┌─┴─┐┌──────────────┐┌─┐│ If-0  └─────────────┘│ Else-0  └┬───────────┬┘┌─┴─┐┌──────────────┐┌─┐└──╥─── └──────────────┘└──────── └──────────────┘ ───────┘   End-0 │
   a: ┤ Ry(2.0587) ├┤ X ├┤ Ry(-0.64847) ├┤M├┤       ───────────────┤         ─┤ Ry(0.841) ├─┤ X ├┤ Ry(-0.37866) ├┤M├───╫────────────────────────────────────────────────────────        ├─
      └────────────┘└───┘└──────────────┘└╥┘└──╥───                └────────  └───────────┘ └───┘└──────────────┘└╥┘   ║                                                         ───────┘
                                          ║ ┌──╨──┐                                                               ║    ║
m1: 1/════════════════════════════════════╩═╡ 0x1 ╞═══════════════════════════════════════════════════════════════╬════╬══════════
                                          0 └─────┘                                                               ║ ┌──╨──┐
m2: 1/════════════════════════════════════════════════════════════════════════════════════════════════════════════╩═╡ 0x0 ╞═══════
```

### 4.5 Verification summary

| Check | A (qutrit) | A (Qiskit) | B | C |
|---|---|---|---|---|
| Path Kraus operators $=K_k$ exactly, 1000+ random parameter sets | ✓ (SymPy) | ✓ | ✓ | ✓ |
| Aer, $10^6$ shots: outcome probabilities (largest \|z\|) | — | 0.95 | 1.46 | 1.32 |
| Aer: post-measurement fidelity | — | $1-10^{-12}$ | $1-10^{-12}$ | $1-10^{-12}$ |
| Instrument tomography on the ideal simulator: distance to $\mathcal T$ | — | 0.0065 | 0.0076 | 0.0070 |

In the tomography test, the pairwise two-sample test of "same instrument" gives p = 0.50 (A–B), 0.30 (A–C) and 0.63 (B–C), with 50,000 shots × 18 settings per scheme. All three implementations realise the same instrument $\mathcal T$.

Every formula and number in Sections 3–5 was re-checked by `verify_report.py`, with 47 of 47 checks passing. This covers:

- the symbolic identities;
- the example matrices;
- $U_6$;
- each Qiskit gate list, through `Operator`, on 2000 random parameter sets;
- the Scheme C closed forms;
- the gate angles printed in the drawings;
- the resource counts and the ancilla-dimension formulas.

---

## 5. Resource comparison: adaptivity saves ancilla dimension

All counts are taken from the actual Qiskit circuits. "Expected" values are averaged over branches for a maximally mixed input; they depend on the input because the branch probabilities are $\langle E_k\rangle$.

| Resource | A (qutrit) | A (Qiskit) | **B** | **C** |
|---|---|---|---|---|
| Ancilla dimension | 3 | 4 (one level wasted) | **2** | **2** |
| Total qubits (Hilbert-space dimension) | qubit + qutrit (6) | 3 (8) | **2 (4)** | **2 (4)** |
| CNOTs, worst / expected | n/a | 5 (4 with XOR) | **2 / 1.47** | **2 / 1.77** |
| Single-qubit gates, worst / expected | n/a | 8 | 5 / 3.40 | 7 / 6.30 |
| Qubit readouts, worst / expected | 1 qutrit | 2 | 2 / 1.47 | 2 / 1.77 |
| Mid-circuit measurements | 0 | 0 | ≤ 1 | ≤ 1 |
| Feed-forward decisions, worst / expected | 0 | 0 | 1 / 0.47 | 2 / 1.77 |
| Resets | 0 | 0 | 0 (absorbed) | 0 (not needed) |
| Required couplings | — | triangle s–a1, s–a2, a1–a2 | one edge | one edge |
| Compiled CNOTs on a line | — | 6 (5 with XOR) | ≤ 2 | ≤ 2 |
| Classical precomputation of angles | closed form | closed form | closed form | matrix square roots and polar decompositions |

**Where the savings come from**

1. **Ancilla dimension.**
   - A one-shot dilation must resolve all $r_{\mathcal T}=3$ Kraus operators in one measurement, so $d_A\ge3$. On qubit hardware this rounds up to two qubits ($d_A=4$), and one level is wasted.
   - An $N$-step adaptive sequence only needs each step to resolve its own outcomes. For two steps with an unchanged system dimension, $r_{\mathcal T}\le d_A^2$ (Eq. (52) of Ref. [1] with $g=1$), giving $d_A\ge\lceil\sqrt3\,\rceil=2$. This is also the Theorem 2 value, and both B and C reach it.
   - In general, an $N$-step sequence obtained by applying Theorem 1 recursively (Corollary 3.1 of Ref. [1]) needs $d_A=g\lceil(r_{\mathcal T}/g^{N-1})^{1/N}\rceil$ (Eq. (44)). Equivalently, with $n_A$ ancilla qubits, $N\,n_A\ge\log_2(g\,r_{\mathcal T})$ (Corollary 3.2): steps can be traded against ancilla qubits.
2. **Gates and connectivity.** Scheme A needs a doubly controlled rotation (3 CNOTs) plus a relabelling CNOT on a fully connected triangle of qubits. B and C replace the extra control by a measurement and classical control. They need at most 2 CNOTs on a single coupled pair.
3. **Fewer operations on average.** The second step only runs on one branch.
4. **The price.** B and C need mid-circuit measurement, real-time feed-forward, and the system qubit idles while it waits for results. Worst-case wall-clock time per shot is longer for B and C (about 3.3 and 3.8 µs, against 2.45 µs for A with the XOR variant and 2.75 µs with 5 CNOTs, with 1 µs readout and 0.5 µs feed-forward latency).

In noisy simulations, B and C reproduce $\mathcal T$ more accurately than A. Both columns use 0.05% single-qubit and 1% CNOT depolarising error and 1.5% readout error. The TVD column comes from an exact density-matrix simulation that also includes idling of the system qubit during readout (1 µs) and feed-forward (0.5 µs). The Choi-distance column comes from instrument tomography in Aer with the same gate and readout noise, but without idling.

| | Outcome-distribution error (TVD) | Choi distance to ideal (tomography) |
|---|---|---|
| A (Qiskit, 5 CNOTs) | 0.029 | 0.073 |
| B | 0.008 | 0.043 |
| C | 0.019 | 0.047 |

Only when readout is very slow (about 4 µs) does A, which never waits for a mid-circuit result, catch up in post-measurement fidelity. B and C still give better outcome statistics.

---

## 6. Theorem 1 is not always optimal

Scheme C is what Theorem 1 produces; Scheme B is a hand-designed adaptive circuit. Comparing them shows that Theorem 1 guarantees existence, and reaches the optimal ancilla dimension within its class, but does not in general minimise the practical cost. The reasons are structural.

1. **The first step is forced to be a Lüders instrument.** Theorem 1 always starts with $\sqrt{\mathsf B_j}$. Whenever $\sqrt{\mathsf B_j}$ is not already the Kraus operator we want, the residual instrument must undo it with $K_k\mathsf B_j^{-1/2}$. This brings in:
   - basis rotations into the eigenbasis of $\mathsf B$, here $R(\pm\xi)$;
   - polar corrections $R(\kappa_k)$, which are extra feed-forward operations after the final measurement.

   A first step with Kraus operators $U_j\sqrt{\mathsf B_j}$ (any unitaries $U_j$) works equally well. For example, take the "natural" first step $\{K_0,\sqrt{\mathbb 1-K_0^\dagger K_0}\}$, as Scheme B does. Then outcome 0 needs no correction at all. Theorem 1 does not include such non-Lüders initial instruments, although Theorem 3 of Ref. [1] uses a related freedom: there the initial Kraus operators are $M_j\sqrt{\mathsf B_j}$ with an isometry $M_j$ onto a smaller intermediate space.
2. **The grouping is a free choice, and the theorem does not say which is best.** Every coarse-graining $\nu$ gives a valid sequence. They all have the same ancilla dimension here, but different expected costs. How often the second step runs depends on the probability of the group measured alone in step 1. Scheme B is in fact also a Theorem 1 sequence, for the grouping $\{0\}|\{1,2\}$: because $K_0\ge0$ here, $\sqrt{E_0}=K_0$ and $\sqrt{E_1+E_2}=S$. It is cheaper than the $\{0,1\}|\{2\}$ grouping of Scheme C:
   - fewer CNOTs and single-qubit gates on average;
   - one feed-forward decision instead of two;
   - no correction after the last measurement;
   - closed-form angles.

   The best grouping depends on the instrument and the input, so it has to be optimised on top of the theorem.
3. **The theorem optimises ancilla dimension only.** Theorem 2 gives the minimal $d_A$ among sequences built from Theorem 1. It says nothing about:
   - two-qubit gate count;
   - circuit depth;
   - number of feed-forward rounds;
   - classical latency;
   - classical precomputation, which needs matrix square roots, inverses and polar decompositions.

   On NISQ hardware these quantities often dominate the error budget.
4. **Theorem 1 covers only part of all adaptive realisations.** It requires the first-step POVM to be a post-processing of the induced POVM. As Ref. [1] notes in its Section 6 (Proposition 1), a 2-step sequence with a Lüders initial instrument of $\mathsf B$ exists if and only if $\mathsf B$ is *compatible* with $\mathcal T$, and compatible POVMs need not be post-processings of $A^{\mathcal T}$. So any compatible POVM can serve as the first step of a 2-step sequence. So some valid, and possibly cheaper, adaptive circuits are not reachable from Theorem 1 at all.
5. **Extra Kraus operators when the dimension shrinks.** If $\dim\mathcal H>\dim\mathcal K$ and some $\mathsf B_j$ is rank deficient, the kernel term of Eq. (18) does not fit into the output space. The residual instruments then need additional Kraus operators (Corollary 1.1, Eq. (23) of Ref. [1]). This can remove the ancilla advantage, and Theorem 3 of Ref. [1] is needed to recover it. That case does not arise for $\mathcal T$, but it shows that the plain construction is not uniformly optimal.

**Conclusion.** Adaptivity provably lowers the ancilla requirement: from a qutrit (two qubits on hardware) to a single reusable qubit for this rank-3 instrument. Theorem 1 gives a systematic, always-applicable recipe that achieves this minimal ancilla dimension. Among the circuits that reach that dimension, though, the cheapest one in gates, feed-forward and noise comes from also choosing the grouping, and the form of the first instrument, to suit the instrument at hand.

---

**References**

[1] S. Sau and M. Sedlák, *Sequential realization of quantum instruments*, Quantum (journal version in the project files). Results used: Defs. 1–6, Theorem 1 with its proof (Eqs. (10), (18)–(20)), Remarks 2–3, Corollary 1.1, Eq. (28), Theorem 2, Corollaries 3.1–3.2 (Eq. (44)), Eq. (52), Theorem 3, Proposition 1.

[2] M. Ozawa, *Quantum measuring processes of continuous observables*, J. Math. Phys. 25, 79 (1984): dilation of instruments.

[3] R. Iten, R. Colbeck, M. Christandl, *Quantum circuits for quantum channels*, Phys. Rev. A 95, 052316 (2017).

[4] M. DeCross et al., *Qubit-reuse compilation with mid-circuit measurement and reset*, Phys. Rev. X 13, 041057 (2023).
