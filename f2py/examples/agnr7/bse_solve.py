import numpy as np
import time


def four_polarization(
    alpha: float,
    nm_dev: int,
    nen: int,
    en: np.ndarray,
    nop: int,
    ndiag: int,
    G_lesser: np.ndarray,
    G_greater: np.ndarray,
    G_retarded: np.ndarray,
    i: int,
    j: int,
    k: int,
    l: int,
) -> complex:
    """
    Computes the P4 IPA tensor using Green's functions G_lesser, G_greater, and G_retarded.

    Args:
        alpha (float): Mixing parameter.
        nm_dev (int): Dimension of the matrix.
        nen (int): Number of energy points.
        en (np.ndarray): Energy values.
        nop (int): Offset for energy index.
        ndiag (int): Not used in this function.
        G_lesser (np.ndarray): Lesser Green's function.
        G_greater (np.ndarray): Greater Green's function.
        G_retarded (np.ndarray): Retarded Green's function.
        i, j, k, l (int): Indices.

    Returns:
        complex: P4 IPA tensor value.
    """
    dE = en[1] - en[0]
    weights = dE / (2 * np.pi)
    L0 = 0.0j

    for ie in range(nop, nen):
        L0 += (1.0 - alpha) * (
            G_lesser[j, l, ie] * np.conj(G_retarded[i, k, ie - nop])
            + G_retarded[j, l, ie] * G_lesser[k, i, ie - nop]
        ) + alpha * 0.5 * (
            G_greater[j, l, ie] * G_lesser[k, i, ie - nop]
            - G_lesser[j, l, ie] * G_greater[k, i, ie - nop]
        )

    L0 *= weights
    return L0


def four_polarization(
    alpha: float,
    nm_dev: int,
    nen: int,
    en: np.ndarray,
    nop: int,
    ndiag: int,
    G_lesser: np.ndarray,
    G_greater: np.ndarray,
    G_retarded: np.ndarray,
    i: int,
    j: int,
    k: int,
    l: int,
) -> complex:
    """
    Computes the P4 IPA tensor based on the given inputs.

    Parameters:
        alpha (float): Weighting factor.
        nm_dev (int): Size of the nm_dev dimension.
        nen (int): Size of the nen dimension.
        en (ndarray): Energy values.
        nop (int): Number of occupied states.
        ndiag (int): Number of diagonal elements.
        G_lesser (ndarray): Lesser Green's function.
        G_greater (ndarray): Greater Green's function.
        G_retarded (ndarray): Retarded Green's function.
        i, j, k, l (int): Indices.

    Returns:
        L0 (complex): Resulting P4 IPA tensor.
    """
    dE = en[1] - en[0]
    weights = dE / (2 * np.pi)

    # Calculate P4_IPA from GG
    L0 = (1.0 - alpha) * (
        np.sum(G_lesser[j, l, (nop + 1) : nen] * np.conj(G_retarded[i, k, : nen - nop]))
        + np.sum(G_retarded[j, l, (nop + 1) : nen] * G_lesser[k, i, : nen - nop])
    )
    L0 += (
        alpha
        * 0.5
        * (
            np.sum(G_greater[j, l, (nop + 1) : nen] * G_lesser[k, i, : nen - nop])
            - np.sum(G_lesser[j, l, (nop + 1) : nen] * G_greater[k, i, : nen - nop])
        )
    )
    L0 *= weights

    return L0


def bse_fullsolve(
    alpha: float,
    spindeg: float,
    nm_dev: int,
    ndiag: int,
    nen: int,
    en: np.ndarray,
    nop: int,
    G_lesser: np.ndarray,
    G_greater: np.ndarray,
    G_retarded: np.ndarray,
    W: np.ndarray,
    V: np.ndarray,
) -> np.ndarray:
    """
    Solves the full Bethe-Salpeter Equation.

    Args:
        alpha (float): Mixing parameter.
        spindeg (float): Spin degeneracy.
        nm_dev (int): Dimension of the matrix.
        ndiag (int): Not used in this function.
        nen (int): Number of energy points.
        en (np.ndarray): Energy values.
        nop (int): Offset for energy index.
        G_lesser (np.ndarray): Lesser Green's function.
        G_greater (np.ndarray): Greater Green's function.
        G_retarded (np.ndarray): Retarded Green's function.
        W (np.ndarray): Static screened Coulomb interaction.
        V (np.ndarray): Bare Coulomb interaction.

    Returns:
        np.ndarray: 2-point polarization function with interacting electron-hole at frequency [[nop]].
    """
    N = nm_dev * nm_dev
    table = np.zeros((2, N), dtype=int)

    # Construct the table of reordered indices
    it = 0
    for i in range(1, nm_dev + 1):
        it += 1
        table[:, it] = [i, i]

    for i in range(1, nm_dev + 1):
        for j in range(1, nm_dev + 1):
            if i != j:
                if abs(i - j) <= ndiag:
                    it += 1
                    table[:, it] = [i, j]

    nn = it
    N = nn  # Resize the problem

    # Initialize matrices
    Lmat = np.zeros((N, N), dtype=complex)
    Mmat = np.zeros((N, N), dtype=complex)
    Amat = np.zeros((N, N), dtype=complex)

    # Compute L0_ijkl = G_jl * G_ki
    tic = time.perf_counter()
    for row in range(N):
        for col in range(N):
            i, j = table[:, row]
            k, l = table[:, col]
            if (abs(i - k) <= ndiag) and (abs(j - l) <= ndiag):
                L0ijkl = four_polarization(
                    alpha,
                    nm_dev,
                    nen,
                    en,
                    nop,
                    ndiag,
                    G_lesser,
                    G_greater,
                    G_retarded,
                    i,
                    j,
                    k,
                    l,
                )
                Lmat[row, col] = L0ijkl * spindeg

    # Construct Mmat
    for row in range(N):
        for col in range(N):
            i, j = table[:, row]
            k, l = table[:, col]
            if i == j and k == l:
                Mmat[row, col] -= 1j * V[i, k] * spindeg
            if i == k and j == l:
                Mmat[row, col] += 1j * W[i, j]
    toc = time.perf_counter()
    print(f"L0_ijkl = G_jl * G_ki: {toc - tic} seconds", flush=True)

    # Compute (I - L0 K) -> A
    tic = time.perf_counter()
    np.fill_diagonal(Amat, 1.0 + 0.0j)
    Amat -= Lmat @ Mmat
    toc = time.perf_counter()
    print(f"-L0 K: {toc - tic} seconds", flush=True)

    # Invert (I - L0 K)
    tic = time.perf_counter()
    Amat_inv = np.linalg.inv(Amat)
    toc = time.perf_counter()
    print(f"(I - L0 K)^-1: {toc - tic} seconds", flush=True)

    # Compute L = (I - L0 K) \ L0
    tic = time.perf_counter()
    Lmat = Amat_inv @ Lmat
    toc = time.perf_counter()
    print(f"L = (I - L0 K) \\ L0: {toc - tic} seconds", flush=True)

    # Compute P_retarded
    P_retarded = np.zeros((nm_dev, nm_dev), dtype=complex)
    for row in range(nm_dev):
        for col in range(nm_dev):
            i, k = table[:, row]
            P_retarded[i, k] = -1j * Mmat[row, col]

    return P_retarded
