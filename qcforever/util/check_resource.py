import math
import sys
import time

import psutil


GIB = 1024**3
MIB = 1024**2
MIN_AUTO_MEMORY = GIB
MEMORY_SAFETY_FRACTION = 0.80


def native_thread_environment(cores):
    """Native relaxation uses the per-candidate allocation, not model threads.

    Set only on the native subprocess (xTB), or temporarily around the legacy
    Gaussian adapter (PM6), restoring the caller environment afterwards.
    """
    return {key: str(cores) for key in (
        'OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')}


def respec_cores(spec_cores):

    ava_cores = get_ava_cores()

    if ava_cores < spec_cores:
        respec_cores = ava_cores
        print(f'{spec_cores} cores you specifed are not available! \n'
                f'The number of cores for QC is reduced to {ava_cores}.')
    else:
        respec_cores = spec_cores

    return respec_cores


def get_ava_cores():

    ava_cores = int(psutil.cpu_count()*(1-psutil.cpu_percent(interval=1)/100))

    if ava_cores == 0:
        while ava_cores == 0:
            time.sleep(10)
            print('Waiting resource.....')
            ava_cores = int(psutil.cpu_count()*(1-psutil.cpu_percent(interval=None)/100))

    return ava_cores


def _atomic_number(element):
    """Return an atomic number without requiring RDKit to be imported."""
    from basis_set_exchange import lut

    return lut.element_Z_from_sym(str(element).strip())


def _element_basis_information(basis, atomic_number):
    """Return spherical AO count and ECP core electrons for one element."""
    import basis_set_exchange as bse

    basis_data = bse.get_basis(
        basis,
        elements=[atomic_number],
        uncontract_spdf=True,
    )
    element_data = basis_data['elements'][str(atomic_number)]

    n_basis = 0
    for shell in element_data.get('electron_shells', []):
        angular_momenta = shell['angular_momentum']
        coefficients = shell['coefficients']

        # uncontract_spdf normally leaves one angular momentum per shell.
        # Retain support for a combined shell in case a basis definition does
        # not permit that transformation.
        if len(angular_momenta) == 1:
            n_contractions = len(coefficients)
            n_basis += n_contractions * (2 * angular_momenta[0] + 1)
        else:
            n_per_am = max(1, len(coefficients) // len(angular_momenta))
            n_basis += sum(n_per_am * (2 * am + 1) for am in angular_momenta)

    return n_basis, int(element_data.get('ecp_electrons', 0))


def molecular_basis_information(atoms, basis, total_charge=0):
    """Count molecular basis functions and explicitly treated electrons.

    Spherical harmonic functions are used for the count.  The difference from
    Cartesian functions is covered by the safety factor in the approximate
    memory estimate.
    """
    if not atoms:
        return None

    atomic_numbers = [_atomic_number(atom) for atom in atoms]
    element_information = {}
    for atomic_number in set(atomic_numbers):
        element_information[atomic_number] = _element_basis_information(
            basis, atomic_number
        )

    n_basis = sum(element_information[z][0] for z in atomic_numbers)
    ecp_electrons = sum(element_information[z][1] for z in atomic_numbers)
    n_electrons = sum(atomic_numbers) - int(total_charge) - ecp_electrons

    if n_basis <= 0 or n_electrons <= 0:
        raise ValueError('The basis-function or electron count is not positive.')

    return {
        'n_basis': n_basis,
        'n_electrons': n_electrons,
        'ecp_electrons': ecp_electrons,
    }


def estimate_memory_bytes(
    atoms,
    basis,
    total_charge=0,
    nproc=1,
    options=None,
):
    """Return a conservative, approximate peak memory requirement.

    This is intentionally an estimate rather than a program-specific exact
    formula.  Dense AO matrices give the main N_basis**2 term, while response
    calculations also require occupied-virtual trial vectors.
    """
    option_names = {str(option).lower() for option in (options or [])}

    try:
        information = molecular_basis_information(atoms, basis, total_charge)
    except (ImportError, KeyError, TypeError, ValueError) as exc:
        print(
            f'Basis information for {basis} is unavailable ({exc}). '
            'The automatic memory fallback will be used.'
        )
        return MIN_AUTO_MEMORY, None

    n_basis = information['n_basis']
    n_electrons = information['n_electrons']
    n_occupied = min(n_basis, math.ceil(n_electrons / 2))
    n_virtual = max(0, n_basis - n_occupied)

    matrix_bytes = 8 * n_basis**2
    memory_bytes = 512 * MIB
    memory_bytes += 96 * matrix_bytes
    memory_bytes += max(1, int(nproc)) * 64 * MIB

    td_options = {'uv', 'fluor', 'tadf', 'nac'}
    response_options = td_options | {'polar', 'nmr'}
    if option_names & response_options:
        memory_bytes += 8 * n_occupied * n_virtual * 40
    if option_names & td_options:
        memory_bytes *= 1.25
    if 'freq' in option_names:
        memory_bytes *= 1.50

    # Covers Cartesian-vs-spherical differences, integral buffers and
    # variation between Gaussian and GAMESS implementations.
    memory_bytes *= 1.50

    return max(MIN_AUTO_MEMORY, math.ceil(memory_bytes)), information


def _format_memory(memory_bytes, round_up=True):
    """Format memory using whole units accepted by Gaussian and QCforever."""
    if memory_bytes >= GIB:
        rounding = math.ceil if round_up else math.floor
        return f'{max(1, rounding(memory_bytes / GIB))}GB'
    rounding = math.ceil if round_up else math.floor
    return f'{max(256, rounding(memory_bytes / (256 * MIB)) * 256)}MB'


def recommend_memory(
    atoms,
    basis,
    total_charge=0,
    nproc=1,
    options=None,
    available_memory_bytes=None,
):
    """Recommend input memory, limited to 80% of currently available RAM."""
    estimated_bytes, information = estimate_memory_bytes(
        atoms=atoms,
        basis=basis,
        total_charge=total_charge,
        nproc=nproc,
        options=options,
    )

    if available_memory_bytes is None:
        available_memory_bytes = psutil.virtual_memory().available
    safe_available_bytes = max(256 * MIB, int(
        available_memory_bytes * MEMORY_SAFETY_FRACTION
    ))

    if estimated_bytes > safe_available_bytes:
        assigned_bytes = safe_available_bytes
        print(
            f'Estimated QC memory ({_format_memory(estimated_bytes)}) exceeds '
            f'the safe available memory '
            f'({_format_memory(safe_available_bytes, round_up=False)}). '
            'The input memory is limited to the available amount; the '
            'calculation may still fail because of insufficient memory.'
        )
    else:
        assigned_bytes = estimated_bytes

    assigned_memory = _format_memory(
        assigned_bytes,
        round_up=estimated_bytes <= safe_available_bytes,
    )
    if information is None:
        print(f'Automatically assigned QC memory: {assigned_memory}.')
    else:
        print(
            f'Automatically assigned QC memory: {assigned_memory} '
            f'(basis functions: {information["n_basis"]}, '
            f'electrons: {information["n_electrons"]}).'
        )

    return assigned_memory


def respec_memory(
    spec_memory,
    atoms,
    basis,
    total_charge=0,
    nproc=1,
    options=None,
    available_memory_bytes=None,
):
    """Keep user-specified memory or provide an automatic recommendation."""
    if spec_memory is not None and str(spec_memory).strip():
        return spec_memory

    return recommend_memory(
        atoms=atoms,
        basis=basis,
        total_charge=total_charge,
        nproc=nproc,
        options=options,
        available_memory_bytes=available_memory_bytes,
    )


if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]} number_of_cores'

    try:
        spec_cores = int(sys.argv[1])
    except (IndexError, ValueError):
        print(usage)
        sys.exit(1)

    print(respec_cores(spec_cores))
