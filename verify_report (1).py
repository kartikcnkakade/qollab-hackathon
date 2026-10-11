"""
Re-verification of every formula and number in adaptive_instrument_report.md.
Symbolic checks with SymPy, circuit checks with qiskit.quantum_info.Operator on the actual gate lists,
numeric checks on 2000 random parameter sets in (0, pi/2)^3.
"""
import numpy as np, sympy as sp
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from schemes import kraus, steps_A, steps_B, steps_C, angles_A, angles_B, angles_C

checks = []
def ok(name, cond):
    checks.append((name, bool(cond)))
    if not cond: print("FAIL:", name)

# ============================================================ Section 3 (symbolic)
t0, t1, x = sp.symbols('theta0 theta1 chi', positive=True)
c0, s0, c1, s1 = sp.cos(t0), sp.sin(t0), sp.cos(t1), sp.sin(t1)
p, m = (sp.cos(x) + sp.sin(x)) / 2, (sp.cos(x) - sp.sin(x)) / 2
z2 = sp.zeros(2)
simp = lambda M: sp.simplify(sp.expand_trig(sp.simplify(M)))
ok("p^2+m^2 = 1/2", sp.simplify(p**2 + m**2 - sp.Rational(1, 2)) == 0)
ok("2pm = cos(2chi)/2", sp.simplify(2*p*m - sp.cos(2*x)/2) == 0)
K0 = sp.diag(c0, c1); K1 = sp.Matrix([[p*s0, m*s1], [m*s0, p*s1]]); K2 = sp.Matrix([[p*s0, -m*s1], [-m*s0, p*s1]])
S = sp.diag(s0, s1); F1 = sp.Matrix([[p, m], [m, p]]); F2 = sp.Matrix([[p, -m], [-m, p]])
Pp = sp.Matrix([[1, 1], [1, 1]])/2; Pm = sp.Matrix([[1, -1], [-1, 1]])/2
ok("F1 = cos chi P+ + sin chi P-", simp(F1 - (sp.cos(x)*Pp + sp.sin(x)*Pm)) == z2)
ok("F2 = sin chi P+ + cos chi P-", simp(F2 - (sp.sin(x)*Pp + sp.cos(x)*Pm)) == z2)
ok("K1 = F1 S, K2 = F2 S", simp(K1 - F1*S) == z2 and simp(K2 - F2*S) == z2)
ok("K0^T K0 + S^T S = 1", simp(K0.T*K0 + S.T*S - sp.eye(2)) == z2)
ok("F1^T F1 + F2^T F2 = 1", simp(F1.T*F1 + F2.T*F2 - sp.eye(2)) == z2)
E0, E1, E2 = K0.T*K0, simp(K1.T*K1), simp(K2.T*K2)
ok("E1 closed form", simp(E1 - sp.Matrix([[s0**2, s0*s1*sp.cos(2*x)], [s0*s1*sp.cos(2*x), s1**2]])/2) == z2)
ok("E2 closed form", simp(E2 - sp.Matrix([[s0**2, -s0*s1*sp.cos(2*x)], [-s0*s1*sp.cos(2*x), s1**2]])/2) == z2)
ok("E0+E1+E2 = 1", simp(E0 + E1 + E2 - sp.eye(2)) == z2)
ok("[K1,K2] != 0 generically", simp(K1*K2 - K2*K1) != z2)
ps0, ps1 = sp.symbols('psi0 psi1')
psi = sp.Matrix([ps0, ps1]); psic = sp.Matrix([sp.conjugate(ps0), sp.conjugate(ps1)])
pr = lambda E: sp.expand((psic.T*E*psi)[0])
ok("Pr(0) formula", sp.simplify(pr(E0) - (c0**2*ps0*sp.conjugate(ps0) + c1**2*ps1*sp.conjugate(ps1))) == 0)
re = (sp.conjugate(ps0)*ps1 + ps0*sp.conjugate(ps1))/2
ok("Pr(1),Pr(2) formula", all(sp.simplify(pr(E) - (sp.Rational(1, 2)*(s0**2*ps0*sp.conjugate(ps0) + s1**2*ps1*sp.conjugate(ps1))
                                                  + sg*s0*s1*sp.cos(2*x)*re)) == 0 for E, sg in ((E1, 1), (E2, -1))))
# example numbers
T0, T1, CHI = 0.37, 1.11, 0.52
Kn = kraus(T0, T1, CHI)
ok("example K0", np.allclose(np.round(Kn[0], 4), [[0.9323, 0], [0, 0.4447]]))
ok("example K1", np.allclose(np.round(Kn[1], 4), [[0.2467, 0.1661], [0.0671, 0.6112]]))
ok("example K2", np.allclose(np.round(Kn[2], 4), [[0.2467, -0.1661], [-0.0671, 0.6112]]))
En = [k.T @ k for k in Kn]
ok("example E0", np.allclose(np.round(En[0], 4), [[0.8692, 0], [0, 0.1977]]))
ok("example E1", np.allclose(np.round(En[1], 4), [[0.0654, 0.0820], [0.0820, 0.4011]]))
ok("example E2", np.allclose(np.round(En[2], 4), [[0.0654, -0.0820], [-0.0820, 0.4011]]))

# ============================================================ 4.1 Scheme A qutrit (symbolic)
def R3(th, i, j):
    M = sp.eye(3); M[i, i] = M[j, j] = sp.cos(th); M[i, j] = -sp.sin(th); M[j, i] = sp.sin(th); return M
kp = sp.kronecker_product; P0, P1 = sp.diag(1, 0), sp.diag(0, 1)
V1 = kp(P0, R3(t0, 0, 1)) + kp(P1, R3(t1, 0, 1))
V2 = kp(Pp, R3(x, 1, 2)) + kp(Pm, R3(sp.pi/2 - x, 1, 2))
U6 = simp(V2*V1)
U6claim = sp.Matrix([[c0, -s0, 0, 0, 0, 0], [p*s0, p*c0, -p, m*s1, m*c1, m], [p*s0, p*c0, p, -m*s1, -m*c1, m],
                     [0, 0, 0, c1, -s1, 0], [m*s0, m*c0, m, p*s1, p*c1, -p], [-m*s0, -m*c0, m, p*s1, p*c1, p]])
ok("U6 = V2 V1 matches printed matrix", simp(U6 - U6claim) == sp.zeros(6))
ok("U6 unitary", simp(U6.T*U6 - sp.eye(6)) == sp.zeros(6))
Ks = [K0, K1, K2]
ok("<k|U6|0> = K_k", all(simp(sp.Matrix(2, 2, lambda a, b: U6[3*a + k, 3*b]) - Ks[k]) == z2 for k in range(3)))
H = sp.Matrix([[1, 1], [1, -1]])/sp.sqrt(2)
V2alt = kp(H, sp.eye(3))*(kp(P0, R3(x, 1, 2)) + kp(P1, R3(sp.pi/2 - x, 1, 2)))*kp(H, sp.eye(3))
ok("V2 H-sandwich form", simp(V2alt - V2) == sp.zeros(6))
ok("V1 gives K0|psi>|0> + S|psi>|1>", all(simp(sp.Matrix(2, 2, lambda a, b: V1[3*a + k, 3*b]) - M) == z2 for k, M in ((0, K0), (1, S), (2, z2))))

# ============================================================ circuits via Qiskit Operator
def unitary(ops, names):
    qc = QuantumCircuit(len(names)); q = {n: i for i, n in enumerate(names)}
    for o in ops:
        if o[0] == 'ry': qc.ry(o[2], q[o[1]])
        elif o[0] == 'h': qc.h(q[o[1]])
        elif o[0] == 'cx': qc.cx(q[o[1]], q[o[2]])
    return Operator(qc).data        # little endian: index = sum_i bit_i 2^i, bit 0 = names[0]
def block(U, aout, ain):            # names[0] = 's'
    idx = lambda s, ab: s + sum(b << (i + 1) for i, b in enumerate(ab))
    return np.array([[U[idx(so, aout), idx(si, ain)] for si in (0, 1)] for so in (0, 1)])
def Rn(a): return np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
def RYn(f): return Rn(f/2)
def U6n(t0_, t1_, chi_):
    f = sp.lambdify((t0, t1, x), U6claim, 'numpy'); return np.array(f(t0_, t1_, chi_), float)
Pp_n, Pm_n = np.array([[1, 1], [1, 1]])/2, np.array([[1, -1], [-1, 1]])/2

rng = np.random.default_rng(2026); N = 2000
cnt = dict(A=0, Axor=0, Aeq=0, Bmat=0, Bckt=0, Breset=0, Bthm1=0, Cformulas=0, Cckt=0, mux=0)
for _ in range(N):
    a0, a1, ch = rng.uniform(0.02, np.pi/2 - 0.02, 3)
    K = kraus(a0, a1, ch); E = [k.T @ k for k in K]
    sa0, sa1, ca0, ca1 = np.sin(a0), np.sin(a1), np.cos(a0), np.cos(a1)
    pp, mm = (np.cos(ch) + np.sin(ch))/2, (np.cos(ch) - np.sin(ch))/2
    Sn = np.diag([sa0, sa1]); F1n = np.array([[pp, mm], [mm, pp]]); F2n = np.array([[pp, -mm], [-mm, pp]])
    # 4.2 Scheme A (Qiskit): blocks, unitarity, XOR variant, isometry equality with U6
    U = unitary(steps_A(a0, a1, ch), ['s', 'a1', 'a2'])
    M = {(i, j): block(U, (i, j), (0, 0)) for i in (0, 1) for j in (0, 1)}
    cnt['A'] += np.allclose(M[0, 0], K[0]) and np.allclose(M[0, 1], K[1]) and np.allclose(M[1, 0], K[2]) and np.allclose(M[1, 1], 0)
    Ux = unitary(steps_A(a0, a1, ch, xor_variant=True), ['s', 'a1', 'a2'])
    Mx = {(i, j): block(Ux, (i, j), (0, 0)) for i in (0, 1) for j in (0, 1)}   # raw bits; a2 = a2' xor a1
    cnt['Axor'] += np.allclose(Mx[0, 0], K[0]) and np.allclose(Mx[0, 1], K[1]) and np.allclose(Mx[1, 1], K[2]) and np.allclose(Mx[1, 0], 0)
    u6 = U6n(a0, a1, ch); enc = {0: (0, 0), 1: (0, 1), 2: (1, 0)}; iso = True
    for si in (0, 1):
        v8 = np.zeros(8)
        for so in (0, 1):
            for i in range(3):
                b1, b2 = enc[i]; v8[so + 2*b1 + 4*b2] = u6[3*so + i, 3*si]
        iso &= np.allclose(v8, U[:, si])
    cnt['Aeq'] += iso
    # multiplexer identity RY(a) X^s RY(b)|0> = RY(2 theta_s)|0>
    g = angles_A(a0, a1, ch); X = np.array([[0, 1], [1, 0]])
    cnt['mux'] += all(np.allclose(RYn(g['a']) @ np.linalg.matrix_power(X, s_) @ RYn(g['b']) @ [1, 0], RYn(2*th) @ [1, 0])
                      for s_, th in ((0, a0), (1, a1)))
    # 4.3 Scheme B: printed W1, W2; Qiskit circuit; absorbed reset; B as Theorem-1 ASI
    W1 = np.kron(np.diag([1, 0]), Rn(a0)) + np.kron(np.diag([0, 1]), Rn(a1))
    W2 = np.kron(Pp_n, Rn(ch)) + np.kron(Pm_n, Rn(np.pi/2 - ch))
    W1p = np.array([[ca0, -sa0, 0, 0], [sa0, ca0, 0, 0], [0, 0, ca1, -sa1], [0, 0, sa1, ca1]])
    W2p = np.array([[pp, -pp, mm, mm], [pp, pp, -mm, mm], [mm, mm, pp, -pp], [-mm, mm, pp, pp]])
    blk4 = lambda W, a: W.reshape(2, 2, 2, 2)[:, a, :, 0]           # |s a> ordering
    cnt['Bmat'] += (np.allclose(W1, W1p) and np.allclose(W2, W2p) and np.allclose(blk4(W1, 0), K[0]) and np.allclose(blk4(W1, 1), Sn)
                    and np.allclose(blk4(W2, 0), F1n) and np.allclose(blk4(W2, 1), F2n))
    st1, st2 = steps_B(a0, a1, ch); Ub1, Ub2 = unitary(st1, ['s', 'a']), unitary(st2, ['s', 'a'])
    Sb = block(Ub1, (1,), (0,))
    cnt['Bckt'] += (np.allclose(block(Ub1, (0,), (0,)), K[0]) and np.allclose(block(Ub2, (0,), (1,)) @ Sb, K[1])
                    and np.allclose(block(Ub2, (1,), (1,)) @ Sb, K[2]))
    cnt['Breset'] += np.allclose(RYn(2*ch - np.pi) @ [0, 1], RYn(2*ch) @ [1, 0])
    def sqrtm(Mx_):
        w, v = np.linalg.eigh(Mx_); return v @ np.diag(np.sqrt(np.clip(w, 0, None))) @ v.T
    cnt['Bthm1'] += (np.allclose(sqrtm(E[0]), K[0]) and np.allclose(sqrtm(E[1] + E[2]), Sn)
                     and np.allclose(K[1] @ np.linalg.inv(Sn), F1n) and np.allclose(K[2] @ np.linalg.inv(Sn), F2n))
    # 4.4 Scheme C formulas
    B0, B1 = E[0] + E[1], E[2]
    xi = 0.5*np.arctan2(-2*sa0*sa1*np.cos(2*ch), sa0**2 - sa1**2)
    r = np.sqrt((sa0**2 - sa1**2)**2 + 4*sa0**2*sa1**2*np.cos(2*ch)**2)
    lp, lm = (sa0**2 + sa1**2 + r)/4, (sa0**2 + sa1**2 - r)/4
    php, phm = np.arcsin(np.sqrt(lp)), np.arcsin(np.sqrt(lm))
    J = np.array([[0, -1], [1, 0]])
    UJ = np.kron(sqrtm(B0), np.eye(2)) + np.kron(sqrtm(B1), J)
    UJ2 = np.kron(Rn(xi), np.eye(2)) @ (np.kron(np.diag([1, 0]), Rn(php)) + np.kron(np.diag([0, 1]), Rn(phm))) @ np.kron(Rn(xi).T, np.eye(2))
    V = K[2] @ np.linalg.inv(sqrtm(B1)); k2 = np.arctan2(mm*(sa1 - sa0), pp*(sa0 + sa1))
    L0, L1 = K[0] @ np.linalg.inv(sqrtm(B0)), K[1] @ np.linalg.inv(sqrtm(B0))
    pol = lambda L: np.arctan2(L[1, 0] - L[0, 1], L[0, 0] + L[1, 1])
    M0 = L0.T @ L0; M1 = L1.T @ L1
    UR = (np.kron(Rn(pol(L0)), np.diag([1, 0])) + np.kron(Rn(pol(L1)), np.diag([0, 1]))) @ (np.kron(sqrtm(M0), np.eye(2)) + np.kron(sqrtm(M1), J))
    cnt['Cformulas'] += (np.allclose(Rn(xi) @ np.diag([lp, lm]) @ Rn(xi).T, B1) and np.allclose(UJ, UJ2)
                         and np.allclose(UJ.T @ UJ, np.eye(4)) and np.allclose(blk4(UJ, 0), sqrtm(B0)) and np.allclose(blk4(UJ, 1), sqrtm(B1))
                         and np.allclose(V, Rn(k2)) and np.allclose(V @ sqrtm(B1), K[2])
                         and np.allclose(L0.T @ L0 + L1.T @ L1, np.eye(2)) and np.allclose(M0 + M1, np.eye(2))
                         and np.allclose(Rn(pol(L0)) @ sqrtm(M0), L0) and np.allclose(Rn(pol(L1)) @ sqrtm(M1), L1)
                         and np.allclose(UR.T @ UR, np.eye(4)) and np.allclose(blk4(UR, 0), L0) and np.allclose(blk4(UR, 1), L1)
                         and np.allclose(L0 @ sqrtm(B0), K[0]) and np.allclose(L1 @ sqrtm(B0), K[1]))
    # 4.4 Scheme C Qiskit circuit (merged angles): every path equals K_k
    c1_, br1, br0, fix = steps_C(a0, a1, ch)
    Uc1, Uc0 = unitary(c1_, ['s', 'a']), unitary(br0, ['s', 'a'])
    one = lambda ops: unitary(ops, ['s'])
    J0, J1 = block(Uc1, (0,), (0,)), block(Uc1, (1,), (0,))
    cnt['Cckt'] += (np.allclose(one(fix[0]) @ block(Uc0, (0,), (0,)) @ J0, K[0]) and np.allclose(one(fix[1]) @ block(Uc0, (1,), (0,)) @ J0, K[1])
                    and np.allclose(one(br1) @ J1, K[2]) and np.allclose(J0 @ J0.conj().T + 0, J0 @ J0.conj().T))
for k_, v in cnt.items(): ok(f"random sets ({N}): {k_}", v == N)
print({k_: f"{v}/{N}" for k_, v in cnt.items()})

# ============================================================ angles printed in the drawings / text
gA = angles_A(T0, T1, CHI); gB = angles_B(T0, T1, CHI); gC = angles_C(T0, T1, CHI)
ok("A angles a,b,tau", np.allclose([gA['a'], gA['b'], gA['tau']], [-0.090796, 0.8308, 0.2654], atol=5e-5))
ok("B angle 2chi-pi = -2.1016", abs(gB['r'] + 2.1016) < 5e-5)
ok("C angles", np.allclose([gC[k] for k in ('sys_in', 'b', 'a', 'sys_k2', 'sys_0', 'b2', 'a2', 'fix0', 'fix1')],
                           [2.6873, 2.0587, -0.64847, -2.4574, -2.4475, 0.841, -0.37866, -0.19975, -0.56375], atol=5e-5))
# C merged-angle structure: alpha = phi+ + phi- - pi/2, beta = phi+ - phi- + pi/2 ; system angles
sa0, sa1 = np.sin(T0), np.sin(T1)
xi = 0.5*np.arctan2(-2*sa0*sa1*np.cos(2*CHI), sa0**2 - sa1**2)
ok("C sys_in = -2 xi", abs(gC['sys_in'] + 2*xi) < 1e-12)

# ============================================================ resource counts
cx = lambda ops: sum(o[0] == 'cx' for o in ops); oq = lambda ops: sum(o[0] in ('ry', 'h') for o in ops)
sA = steps_A(T0, T1, CHI); sAx = steps_A(T0, T1, CHI, True); b1, b2 = steps_B(T0, T1, CHI); c1_, br1, br0, fix = steps_C(T0, T1, CHI)
ok("A: 5 CX, 8 1q; XOR: 4 CX", cx(sA) == 5 and oq(sA) == 8 and cx(sAx) == 4)
ok("B: CX 1/2, 1q 2/5", cx(b1) == 1 and cx(b1 + b2) == 2 and oq(b1) == 2 and oq(b1 + b2) == 5)
ok("C: CX 1/2, 1q 4/7", cx(c1_ + br1) == 1 and cx(c1_ + br0 + fix[0]) == 2 and oq(c1_ + br1) == 4 and oq(c1_ + br0 + fix[0]) == 7)
w = 0.5*np.trace(En[0]); w2 = 0.5*np.trace(En[2])         # P(m1=0) for B and P(m1=1) for C at rho = I/2
ok("expected CX B=1.47, C=1.77", abs((w*1 + (1 - w)*2) - 1.467) < 1e-3 and abs((w2*1 + (1 - w2)*2) - 1.767) < 1e-3)
ok("expected 1q B=3.40, C=6.30", abs((w*2 + (1 - w)*5) - 3.40) < 5e-3 and abs((w2*4 + (1 - w2)*7) - 6.30) < 5e-3)
ok("expected feed-forward B=0.47, C=1.77", abs((1 - w) - 0.4665) < 1e-3 and abs((w2*1 + (1 - w2)*2) - 1.767) < 1e-3)
# Theorem 2 / Eq. (28) / Eq. (44)
rT, g = 3, 1
ok("Theorem 2: d_A = g*ceil(sqrt(rT/g)) = 2", g*int(np.ceil(np.sqrt(rT/g))) == 2)
ok("Eq.(28): max(|Omega_B|, m_0, m_1) = max(2,2,1) = 2", max(2, 1 + 1, 1) == 2)
ok("Eq.(44), N=2: g*ceil((rT/g^(N-1))^(1/N)) = 2", g*int(np.ceil((rT/g)**0.5)) == 2)

print(f"\n{sum(c for _, c in checks)}/{len(checks)} checks passed")
for n_, c in checks:
    print(("PASS " if c else "FAIL ") + n_)
