"""Chemical structure normalization, SMILES canonicalization, and InChIKey generation."""


def canonicalize_smiles(smiles: str) -> tuple[str, str]:
    """Canonicalize a SMILES string and generate its standard InChIKey.

    Args:
        smiles: Raw SMILES string.

    Returns:
        A tuple of (canonical_smiles, inchikey).

    Raises:
        ValueError: If SMILES string cannot be parsed by RDKit or is empty.
    """
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string cannot be empty or whitespace.")

    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError(
            "RDKit is required for structure normalization. Install via `pip install rdkit`."
        ) from exc

    clean_smiles = smiles.strip()
    mol = Chem.MolFromSmiles(clean_smiles)
    if mol is None:
        raise ValueError(f"Invalid or unparseable SMILES string: '{clean_smiles}'")

    can_smiles = Chem.MolToSmiles(mol, canonical=True)
    inchikey = Chem.MolToInchiKey(mol)

    if not inchikey:
        raise ValueError(f"Failed to generate InChIKey for SMILES: '{clean_smiles}'")

    return can_smiles, inchikey


def smiles_to_inchikey(smiles: str) -> str:
    """Generate standard InChIKey from a SMILES string.

    Args:
        smiles: Input SMILES string.

    Returns:
        Standard 27-character InChIKey string.
    """
    _, inchikey = canonicalize_smiles(smiles)
    return inchikey
