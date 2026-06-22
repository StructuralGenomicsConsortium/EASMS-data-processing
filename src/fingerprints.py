from typing import Union, List, Callable, Optional
from functools import partial
import inspect
import abc
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import numpy.typing as npt
from tqdm import tqdm

from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, rdMolDescriptors, RDKFingerprint
from rdkit.Chem.AtomPairs import Pairs
from rdkit import Chem

from utils import to_mol, catch_boost_argument_error


# --- Parallel fingerprint generation -------------------------------------
# RDKit holds the GIL during fingerprint computation, so threads give no real
# speedup; we fan out across processes instead. Tunable via the FP_N_JOBS env
# var (defaults to all CPUs). Inputs below _PARALLEL_MIN_ITEMS run serially,
# since process spin-up + pickling would cost more than it saves.
_PARALLEL_MIN_ITEMS = 512


def _default_n_jobs() -> int:
    """Worker count for parallel FP generation: ``$FP_N_JOBS`` if set, else all CPUs."""
    env = os.environ.get("FP_N_JOBS")
    if env:
        try:
            return max(1, int(env))
        except ValueError:
            pass
    return os.cpu_count() or 1


def _chunk(seq: List, n: int) -> List[List]:
    """Split ``seq`` into ``n`` contiguous, order-preserving chunks."""
    n = max(1, min(n, len(seq)))
    k, m = divmod(len(seq), n)
    return [seq[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


def _wrap_handle_none(fp_func: Callable, *args, fail_size: Optional[int] = None, **kwargs) -> List:
    """
    Wraps an FP function and handles RDKit chemical exceptions by returning a list of NaN values.

    Parameters
    ----------
    fp_func (Callable):
        The function to be wrapped.
    *args:
        Variable length argument list to be passed to the function.
    **kwargs:
        Arbitrary keyword arguments to be passed to the function.
    fail_size (int or None, optional):
        The size of the list to be returned in case of failure.
        If None, the fail_size is determined by calling the function with a default argument.
        Defaults to None.

    Returns
    -------
    list:
        A list of NaN values with length equal to fail_size.

    Notes:
    -----
    If not `fail_size` is passed, will try and assume it by calling the FP function on "CCC" to get the FP length
    This can cause major overhead if lost of SMILES fail inside RDKit (what this wrapper is built to catch),
    so if you expect to see high failure rates preseting the fail length will minimize this overhead

    The `FPFunc` Class has a `_dimension` attribute that hold the FP length, thus any defined FPFunc should not
    suffer from this overhead. Any newly added FPFunc should follow this convention as well

    Raises
    ------
    Any Exception:
        If the exception thrown is not a boost C++ exception, will still raise that exception.

    Examples
    ________
    >>> _wrap_handle_none(AllChem.GetMorganFingerprintAsBitVect, Chem.MolFromSmiles("CCC"), 2)
    [nan, nan]

    >>> _wrap_handle_none(AllChem.GetMorganFingerprintAsBitVect, Chem.MolFromSmiles("CCC"), fail_size=3)
    [nan, nan, nan]
    """
    assert isinstance(fp_func, Callable), "fp_func must be a callable"
    try:
        return list(fp_func(*args, **kwargs))
    except Exception as e:  # throws boost C++ exception, which python cannot catch
        if catch_boost_argument_error(e):
            if fail_size is None:
                # attempt to get fail_size from the func if it is not passed
                fail_size = len(list(fp_func(AllChem.MolFromSmiles("CCC"))))
            return [np.nan] * fail_size
        else:
            raise


def _generate_fps_chunk(fp_cls, items, dimension):
    """Worker for parallel FP generation.

    Rebuilds the FP function from its class (cheaper and safer to pickle than the
    underlying RDKit partial) and computes fingerprints for one chunk of inputs.
    Defined at module level so ``ProcessPoolExecutor`` can pickle it.
    """
    fp = fp_cls()
    return [
        _wrap_handle_none(fp._func, to_mol(c), fail_size=dimension)
        for c in items
    ]


class BaseFPFunc(abc.ABC):
    """
    Base class for all FP functions used in any AIRCHECK pipeline

    Parameters
    ----------
    kwargs : dict
        dictionary of keyword arguments to pass to the fingerprint function in from {`argument`: `value`}

    Attributes
    ----------
    _func : functools.partial object
        The callable FP function instance as a partial with static setting arguments (e.g., 'radius') pre-set
    _binary : bool
        Whether the FP function returns binary fingerprints
    _dimension : int
        the dimensionality of the fingerprints that will be generated

    Notes
    -----
    When declaring a child of the `BaseFPFunc` class, the `_func`, `_dimension` and `_binary` attributes must be set
    during instantiation of the child.
    FP Funcs operate on rdkit.ROMol objects, not smiles and will fail if SMILES are passed
    """
    def __init__(self, **kwargs):
        self._kwargs = kwargs
        self._func: Callable = lambda: None
        self._binary: bool = False
        self._dimension: int = -1

    def __call__(self, smis, *args, use_tqdm: bool = False,
                 n_jobs: Optional[int] = None, **kwargs) -> npt.NDArray[np.int32]:
        items = list(np.atleast_1d(smis))
        if n_jobs is None:
            n_jobs = _default_n_jobs()
        n_jobs = min(n_jobs, len(items)) if items else 1

        # Parallel path: fan the per-molecule loop out across processes.
        if n_jobs > 1 and len(items) >= _PARALLEL_MIN_ITEMS:
            try:
                worker = partial(_generate_fps_chunk, self.__class__,
                                 dimension=self._dimension)
                rows = []
                with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                    for chunk_rows in ex.map(worker, _chunk(items, n_jobs)):
                        rows.extend(chunk_rows)
                return np.array(rows)
            except Exception as e:
                # Never let parallelism break the pipeline -- fall back to serial.
                print(f"  [fingerprints] parallel FP generation failed "
                      f"({type(e).__name__}: {e}); falling back to serial.")

        # Serial path (also used for small inputs / when n_jobs == 1).
        return np.array(
            [
                _wrap_handle_none(self._func, to_mol(c), fail_size=self._dimension)
                for c in tqdm(items, disable=not use_tqdm)
            ]
        )

    def __reduce__(self):
        # Concrete FP funcs take no constructor args and rebuild their RDKit
        # callable in __init__, so they pickle as "call the class". This lets a
        # ProcessPoolExecutor ship the FP func to workers without having to
        # pickle the underlying RDKit partial (which is not reliably picklable).
        return (self.__class__, ())

    def __eq__(self, other) -> bool:
        if isinstance(other, BaseFPFunc):
            if (
                inspect.signature(self._func).parameters
                == inspect.signature(other).parameters
            ):
                return True
        return False

    def generate_fps(
        self,
        smis: Union[str, Chem.rdchem.Mol, List[Union[str, Chem.rdchem.Mol]]],
        use_tqdm: bool = False,
        n_jobs: Optional[int] = None,
    ) -> np.ndarray:
        """
        Generate Fingerprints for a set of smiles
        Parameters
        ----------
        smis : str, rdkit Mol or list of rdkit Mol or str
            the SMILES or Mol objects (or multiple SMILES/Mol objects) you want to generate a fingerprint(s) for
        use_tqdm : bool, default: False
            have a tqdm task to track progress
        n_jobs : int or None, default: None
            number of worker processes to use. None = all CPUs (or $FP_N_JOBS);
            1 forces serial. Small inputs always run serially regardless.

        Returns
        -------
        ndarray
            an array of size (M, d), where M is number of Mols passes and d is the dimension of fingerprint

        Notes
        -----
        The passed list can be a mix of SMILES and Mol objects.
        If the SMILES are invalid or the Mol object(s) are None, then that molecules row of the output fingerprint
         array will be `np.nan` (e.i., the fingerprint for that molecule will be 1-d array of `np.nan` of dimension d)
        This function just wraps the __call__ method of the class

        """
        return self.__call__(smis, use_tqdm=use_tqdm, n_jobs=n_jobs)

    def to_dict(self) -> dict:
        """
        Returns the name and settings of the FP function as a dict

        Returns
        -------
        dict
            name and settings of FP function

        """
        _signature = inspect.signature(self._func)
        args = {
            k: v.default
            for k, v in _signature.parameters.items()
            if v.default is not inspect.Parameter.empty
        }
        args['name'] = self.func_name()
        return args

    def is_binary(self) -> bool:
        """
        Determines if the FP function is binary

        Returns
        -------
        bool:
            True if the FP function is binary, otherwise False

        Notes
        -----
        This function just returns the `_binary` attribute set during instantiation
        """
        return self._binary

    def func_name(self) -> str:
        """
        Returns the name of the RDKit python function used to calculate the fingerprint

        Returns
        -------
        str
            name of python FP function

        Notes
        -----
        This function meant to support reproducibility of FP calculations
        """
        if isinstance(self._func, partial):
            return self._func.func.__name__
        else:
            return self._func.__name__


class HitGenECFP4(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating ECFP4 fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 2, "nBits": 2048, "useFeatures": False})
        self._func = partial(AllChem.GetHashedMorganFingerprint, **self._kwargs)
        self._dimension = 2048


class HitGenECFP6(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating ECFP6 fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": False})
        self._func = partial(AllChem.GetHashedMorganFingerprint, **self._kwargs)
        self._dimension = 2048


class HitGenFCFP4(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating FCFP4 fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 2, "nBits": 2048, "useFeatures": True})
        self._func = partial(AllChem.GetHashedMorganFingerprint, **self._kwargs)
        self._dimension = 2048


class HitGenFCFP6(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating FCFP6 fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": True})
        self._func = partial(AllChem.GetHashedMorganFingerprint, **self._kwargs)
        self._dimension = 2048


class HitGenBinaryECFP4(BaseFPFunc):
    """
    The FP calculation used to match the binary ECFP4 fingerprints generated from HitGen's ECFP4 fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 2, "nBits": 2048, "useFeatures": False})
        self._func = partial(AllChem.GetMorganFingerprintAsBitVect, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenBinaryECFP6(BaseFPFunc):
    """
    The FP calculation used to match the binary ECFP6 fingerprints generated from HitGen's ECFP6 fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": False})
        self._func = partial(AllChem.GetMorganFingerprintAsBitVect, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenBinaryFCFP4(BaseFPFunc):
    """
    The FP calculation used to match the binary FCFP4 fingerprints generated from HitGen's FCFP4 fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 2, "nBits": 2048, "useFeatures": True})
        self._func = partial(AllChem.GetMorganFingerprintAsBitVect, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenBinaryFCFP6(BaseFPFunc):
    """
    The FP calculation used to match the binary FCFP6 fingerprints generated from HitGen's FCFP6 fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": True})
        self._func = partial(AllChem.GetMorganFingerprintAsBitVect, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenMACCS(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating MACCS fingerprints

    Notes
    -----
    Unlike other HitGen FPs, MACCS is only generated in a binary fashion by HitGen, thus no hashed/count version exists
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__()
        self._func = partial(rdMolDescriptors.GetMACCSKeysFingerprint, **self._kwargs)
        self._binary = True
        self._dimension = 167


class HitGenRDK(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating RDK fingerprints

    Notes
    -----
    Unlike other HitGen FPs, RDK is only generated in a binary fashion by HitGen, thus no hashed/count version exists
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"fpSize": 2048})
        self._func = partial(RDKFingerprint, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenAvalon(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating Avalon fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"nBits": 2048})
        self._func = partial(pyAvalonTools.GetAvalonCountFP, **self._kwargs)
        self._dimension = 2048


class HitGenBinaryAvalon(BaseFPFunc):
    """
    The FP calculation used to match the binary Avalon fingerprints generated from HitGen's Avalon fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide.
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": True})
        self._func = partial(pyAvalonTools.GetAvalonFP, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenAtomPair(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating AtomPair fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"nBits": 2048})
        self._func = partial(rdMolDescriptors.GetHashedAtomPairFingerprint, **self._kwargs)
        self._dimension = 2048


class HitGenBinaryAtomPair(BaseFPFunc):
    """
    The FP calculation used to match the binary AtomPair fingerprints generated from HitGen's AtomPair fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide.
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": True})
        self._func = partial(Pairs.GetAtomPairFingerprintAsBitVect, **self._kwargs)
        self._binary = True
        self._dimension = 2048


class HitGenTopTor(BaseFPFunc):
    """
    The FP calculation used by HitGen when generating Topological Torsion (TopTor) fingerprints

    Notes
    -----
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"nBits": 2048})
        self._func = partial(AllChem.GetHashedTopologicalTorsionFingerprint, **self._kwargs)
        self._dimension = 2048


class HitGenBinaryTopTor(BaseFPFunc):
    """
    The FP calculation used to match the binary TopTor fingerprints generated from HitGen's TopTor fingerprints

    Notes
    -----
    This FP is not directly given by HitGen, but can be calculated by just 'binarizing' the hashed FP they provide.
    All settings and attributes are preset during the instantiation of the object.
    Tweaks to FP settings should not be made, as the FP function will not match HitGen anymore
    """
    def __init__(self):
        super().__init__(**{"radius": 3, "nBits": 2048, "useFeatures": True})
        self._func = partial(AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect, **self._kwargs)
        self._binary = True
        self._dimension = 2048
